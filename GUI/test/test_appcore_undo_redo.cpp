#include <gtest/gtest.h>

#include <QCoreApplication>
#include <QDateTime>
#include <QElapsedTimer>
#include <QEventLoop>
#include <QFile>
#include <QGuiApplication>
#include <QString>
#include <QTemporaryDir>
#include <QThread>

#include "viewmodels/AppCore.h"
#include "core/storage/SessionStore.h"

namespace {

int cell_count(cas::gui::NotebookVM* notebook) {
    return notebook ? notebook->cellCount() : 0;
}

void pump_events_for_ms(int duration_ms) {
    QElapsedTimer timer;
    timer.start();
    while (timer.elapsed() < duration_ms) {
        QCoreApplication::processEvents(QEventLoop::AllEvents, 20);
        QThread::msleep(10);
    }
}

} // namespace

TEST(GuiUndoRedo, MultiStepUndoRedoAcrossSessionCellAndPlot) {
    ASSERT_NE(QGuiApplication::instance(), nullptr);

    QTemporaryDir temp_dir;
    ASSERT_TRUE(temp_dir.isValid());
    const QString storage_path = temp_dir.filePath("workspace.json");
    QFile::remove(storage_path);
    qputenv("CAS_GUI_STORAGE_PATH", storage_path.toUtf8());

    AppCore* core = AppCore::instance();
    ASSERT_NE(core, nullptr);

    auto* sessions = core->sessions();
    auto* notebook = core->notebook();
    auto* plot = core->plot();
    ASSERT_NE(sessions, nullptr);
    ASSERT_NE(notebook, nullptr);
    ASSERT_NE(plot, nullptr);

    const int sessions_initial = sessions->titles().size();
    const int cells_initial = cell_count(notebook);
    const QString plot_var_initial = plot->variable();
    const double plot_xmin_initial = plot->xMin();
    const double plot_xmax_initial = plot->xMax();

    // Phase A: plot settings undo/redo (coalesced edit)
    plot->setVariable("u");
    plot->setXMin(-7.0);
    plot->setXMax(7.0);
    pump_events_for_ms(900);
    EXPECT_TRUE(core->undoActionLabel().contains("Edit plot"));
    core->undo();
    EXPECT_EQ(plot->variable(), plot_var_initial);
    EXPECT_DOUBLE_EQ(plot->xMin(), plot_xmin_initial);
    EXPECT_DOUBLE_EQ(plot->xMax(), plot_xmax_initial);
    EXPECT_TRUE(core->redoActionLabel().contains("Edit plot"));
    core->redo();
    EXPECT_EQ(plot->variable(), QString("u"));
    EXPECT_DOUBLE_EQ(plot->xMin(), -7.0);
    EXPECT_DOUBLE_EQ(plot->xMax(), 7.0);

    // Phase B: session creation undo/redo
    core->newSession();
    EXPECT_EQ(sessions->titles().size(), sessions_initial + 1);
    core->undo();
    EXPECT_EQ(sessions->titles().size(), sessions_initial);
    core->redo();
    EXPECT_EQ(sessions->titles().size(), sessions_initial + 1);

    // Phase C: cell add undo/redo
    const int cells_before_add = cell_count(notebook);
    core->addEmptyCell();
    EXPECT_EQ(cell_count(notebook), cells_before_add + 1);
    core->undo();
    EXPECT_EQ(cell_count(notebook), cells_before_add);
    core->redo();
    EXPECT_EQ(cell_count(notebook), cells_before_add + 1);

    // Redo stack must be cleared by a fresh action after an undo.
    core->undo();
    EXPECT_TRUE(core->canRedo());
    core->duplicateCurrentCell();
    EXPECT_FALSE(core->canRedo());

    // Bring cells to baseline state for deterministic continuation.
    while (cell_count(notebook) > cells_initial && core->canUndo()) {
        core->undo();
    }
    EXPECT_EQ(cell_count(notebook), cells_initial);
    EXPECT_EQ(plot->variable(), QString("u"));

    // Phase D: session delete undo/redo
    if (sessions->titles().size() < 2) {
        core->newSession();
    }
    const int sessions_before_delete = sessions->titles().size();
    core->deleteSession(sessions->activeIndex());
    EXPECT_EQ(sessions->titles().size(), sessions_before_delete - 1);

    core->undo();
    EXPECT_EQ(sessions->titles().size(), sessions_before_delete);

    core->redo();
    EXPECT_EQ(sessions->titles().size(), sessions_before_delete - 1);

    // Phase E: persisted workspace must include plot settings.
    plot->setVariable("t");
    plot->setXMin(-3.5);
    plot->setXMax(11.25);
    pump_events_for_ms(1300);

    const auto loaded = cas::gui::storage::SessionStore::load();
    ASSERT_TRUE(loaded.found);
    EXPECT_EQ(loaded.workspace.plot.variable, QString("t"));
    EXPECT_DOUBLE_EQ(loaded.workspace.plot.xMin, -3.5);
    EXPECT_DOUBLE_EQ(loaded.workspace.plot.xMax, 11.25);
}

TEST(GuiUndoRedo, EditCellCoalescedUndo) {
    AppCore* core = AppCore::instance();
    ASSERT_NE(core, nullptr);

    auto* notebook = core->notebook();
    ASSERT_NE(notebook, nullptr);
    const int cells_initial = cell_count(notebook);

    core->addEmptyCell();
    auto* cell = core->notebook()->currentCell();
    ASSERT_NE(cell, nullptr);

    // Initial state
    cell->setInputLatex("first");
    pump_events_for_ms(600); // commit pending

    // Coalesced edits
    cell->setInputLatex("first edit");
    pump_events_for_ms(100);
    cell->setInputLatex("first edit 2");
    pump_events_for_ms(600); // commit

    core->undo(); // Undo "first edit 2" -> "first"
    EXPECT_EQ(core->notebook()->currentCell()->inputLatex(), QString("first"));

    core->redo(); // Redo -> "first edit 2"
    EXPECT_EQ(core->notebook()->currentCell()->inputLatex(), QString("first edit 2"));

    // Cleanup
    core->deleteCurrentCell();
    EXPECT_EQ(cell_count(notebook), cells_initial);
}

TEST(GuiUndoRedo, UndoTransaction) {
    AppCore* core = AppCore::instance();
    ASSERT_NE(core, nullptr);

    auto* sessions = core->sessions();
    const int sessions_initial = sessions->titles().size();

    // Transaction committed
    core->beginUndoTransaction("Test Commit");
    core->newSession();
    core->endUndoTransaction(true);
    EXPECT_EQ(sessions->titles().size(), sessions_initial + 1);

    core->undo();
    EXPECT_EQ(sessions->titles().size(), sessions_initial);

    // Transaction aborted
    core->beginUndoTransaction("Test Abort");
    core->newSession();
    core->endUndoTransaction(false);
    // Even if aborted, the action itself modified the state immediately,
    // but the undo step shouldn't be pushed. Wait, endUndoTransaction(false)
    // just skips pushUndoStep. So we can't undo it.
    EXPECT_EQ(sessions->titles().size(), sessions_initial + 1);
    
    // Undo should do nothing related to "Test Abort", or it might undo something earlier if available
    // But since we just did a new session that wasn't pushed to undo stack,
    // the state is now un-undoable for that last new session.
}

TEST(GuiUndoRedo, DeleteLastCellSessionRobustness) {
    AppCore* core = AppCore::instance();
    ASSERT_NE(core, nullptr);

    auto* notebook = core->notebook();
    auto* sessions = core->sessions();

    // Ensure we have exactly 1 cell
    while (cell_count(notebook) > 1) {
        core->deleteCurrentCell();
    }
    EXPECT_EQ(cell_count(notebook), 1);

    // Ensure we have exactly 1 session
    while (sessions->titles().size() > 1) {
        core->deleteSession(sessions->titles().size() - 1);
    }
    EXPECT_EQ(sessions->titles().size(), 1);

    // Attempt delete last cell
    core->deleteCurrentCell();
    EXPECT_EQ(cell_count(notebook), 1);

    // Attempt delete last session
    core->deleteSession(0);
    EXPECT_EQ(sessions->titles().size(), 1);
}

TEST(GuiUndoRedo, PersistRestoreMultiSession) {
    AppCore* core = AppCore::instance();
    ASSERT_NE(core, nullptr);

    // Create session 1
    core->newSession();
    core->renameActiveSession("Test Session 1");
    core->addEmptyCell(); // Cell 1
    core->notebook()->currentCell()->setInputLatex("x^2");
    pump_events_for_ms(600);
    
    core->addEmptyCell(); // Cell 2
    core->notebook()->currentCell()->setInputLatex("2*x");
    pump_events_for_ms(600);

    // Reorder Cell 2 up
    core->moveCurrentCellUp();
    
    // Create session 2
    core->newSession();
    core->renameActiveSession("Test Session 2");
    core->addEmptyCell();
    core->notebook()->currentCell()->setInputLatex("y^3");
    pump_events_for_ms(600);

    // Wait for auto-save
    pump_events_for_ms(1500);

    // Verify disk state
    const auto loaded = cas::gui::storage::SessionStore::load();
    ASSERT_TRUE(loaded.found);
    
    // Find our two sessions in loaded state.
    // They are added at the end of the sessions list.
    const auto& s1 = loaded.workspace.sessions[loaded.workspace.sessions.size() - 2];
    const auto& s2 = loaded.workspace.sessions[loaded.workspace.sessions.size() - 1];

    EXPECT_EQ(s1.title, QString("Test Session 1"));
    EXPECT_GE(s1.cells.size(), 2);
    // Since Cell 2 ("2*x") moved up, it should be before Cell 1 ("x^2")
    EXPECT_EQ(s1.cells[s1.cells.size() - 2].input, QString("2*x"));
    EXPECT_EQ(s1.cells[s1.cells.size() - 1].input, QString("x^2"));

    EXPECT_EQ(s2.title, QString("Test Session 2"));
    EXPECT_EQ(s2.cells.back().input, QString("y^3"));
}

TEST(GuiUndoRedo, PersistStructuredInspectorPayload) {
    ASSERT_NE(QGuiApplication::instance(), nullptr);

    QTemporaryDir temp_dir;
    ASSERT_TRUE(temp_dir.isValid());
    const QString storage_path = temp_dir.filePath("workspace.json");
    QFile::remove(storage_path);
    qputenv("CAS_GUI_STORAGE_PATH", storage_path.toUtf8());

    AppCore* core = AppCore::instance();
    ASSERT_NE(core, nullptr);

    const QString title =
        QString("Structured Output Session %1").arg(QDateTime::currentMSecsSinceEpoch());
    core->newSession();
    core->renameActiveSession(title);

    auto* current = core->notebook()->currentCell();
    ASSERT_NE(current, nullptr);
    current->setInputLatex("x+0");
    pump_events_for_ms(600);

    core->executeCurrentCell();
    pump_events_for_ms(1500);

    const auto loaded = cas::gui::storage::SessionStore::load();
    ASSERT_TRUE(loaded.found);

    const cas::gui::storage::StoredSession* stored_session = nullptr;
    for (const auto& session : loaded.workspace.sessions) {
        if (session.title == title) {
            stored_session = &session;
            break;
        }
    }
    ASSERT_NE(stored_session, nullptr);

    const cas::gui::storage::StoredCell* executed_cell = nullptr;
    for (const auto& cell : stored_session->cells) {
        if (cell.input == "x+0" && cell.hasOutput) {
            executed_cell = &cell;
            break;
        }
    }
    ASSERT_NE(executed_cell, nullptr);
    EXPECT_EQ(executed_cell->output, QString("x"));
    ASSERT_FALSE(executed_cell->alternatives.isEmpty());
    ASSERT_FALSE(executed_cell->steps.isEmpty());

    bool saw_text_representation = false;
    for (const auto& item : executed_cell->alternatives) {
        const auto map = item.toMap();
        if (map.value("id").toString() == "text" && map.value("value").toString() == "x") {
            saw_text_representation = true;
            break;
        }
    }
    EXPECT_TRUE(saw_text_representation);

    const auto last_step = executed_cell->steps.back().toMap();
    EXPECT_EQ(last_step.value("afterLatex").toString(), QString("x"));
}

TEST(GuiUndoRedo, NotebookInteractionCommandsStayCoherent) {
    AppCore* core = AppCore::instance();
    ASSERT_NE(core, nullptr);

    auto* notebook = core->notebook();
    ASSERT_NE(notebook, nullptr);

    core->reloadWorkspaceFromDisk();
    pump_events_for_ms(400);

    while (notebook->cellCount() > 1) {
        notebook->setSelectedIndex(notebook->cellCount() - 1);
        core->deleteCurrentCell();
    }

    notebook->currentCell()->setInputLatex("alpha");
    pump_events_for_ms(600);

    core->insertCellAbove();
    EXPECT_EQ(notebook->selectedIndex(), 0);
    EXPECT_EQ(notebook->currentCell()->inputLatex(), QString());

    core->selectNextCell();
    EXPECT_EQ(notebook->selectedIndex(), 1);
    EXPECT_EQ(notebook->currentCell()->inputLatex(), QString("alpha"));

    notebook->currentCell()->setInputLatex("alpha+beta");
    pump_events_for_ms(600);
    core->splitCurrentCell(5);
    ASSERT_GE(notebook->cellCount(), 3);
    EXPECT_EQ(notebook->rawCells().at(1)->inputLatex(), QString("alpha"));
    EXPECT_EQ(notebook->rawCells().at(2)->inputLatex(), QString("+beta"));

    notebook->setSelectedIndex(2);
    core->mergeWithPreviousCell();
    ASSERT_GE(notebook->cellCount(), 2);
    EXPECT_EQ(notebook->rawCells().at(1)->inputLatex(), QString("alpha\n+beta"));

    notebook->setSelectedIndex(1);
    core->insertCellBelow();
    EXPECT_EQ(notebook->selectedIndex(), 2);
}

TEST(GuiUndoRedo, ReloadWorkspaceRestoresLiveState) {
    ASSERT_NE(QGuiApplication::instance(), nullptr);

    QTemporaryDir temp_dir;
    ASSERT_TRUE(temp_dir.isValid());
    const QString storage_path = temp_dir.filePath("workspace.json");
    QFile::remove(storage_path);
    qputenv("CAS_GUI_STORAGE_PATH", storage_path.toUtf8());

    AppCore* core = AppCore::instance();
    ASSERT_NE(core, nullptr);

    core->reloadWorkspaceFromDisk();
    pump_events_for_ms(500);

    auto* notebook = core->notebook();
    auto* sessions = core->sessions();
    auto* plot = core->plot();
    ASSERT_NE(notebook, nullptr);
    ASSERT_NE(sessions, nullptr);
    ASSERT_NE(plot, nullptr);

    const QString title =
        QString("Reload Roundtrip Session %1").arg(QDateTime::currentMSecsSinceEpoch());
    core->newSession();
    core->renameActiveSession(title);

    auto* define_cell = notebook->currentCell();
    ASSERT_NE(define_cell, nullptr);
    define_cell->setInputLatex("a = 2");
    pump_events_for_ms(600);
    core->runAndAdvance();
    pump_events_for_ms(1200);

    auto* simplify_cell = notebook->currentCell();
    ASSERT_NE(simplify_cell, nullptr);
    simplify_cell->setInputLatex("x+0");
    pump_events_for_ms(600);
    core->runAndAdvance();
    pump_events_for_ms(1200);

    ASSERT_GE(notebook->cellCount(), 3);
    notebook->setSelectedIndex(1);
    plot->setVariable("u");
    plot->setXMin(-2.5);
    plot->setXMax(4.0);
    core->setSelectedInspectorTab("Steps");
    core->setSelectedRepresentationId("text");
    pump_events_for_ms(1200);

    // Dirty the live state without waiting for autosave, then force a real reload from disk.
    notebook->currentCell()->setInputLatex("corrupt");
    plot->setVariable("corrupt");
    core->reloadWorkspaceFromDisk();
    pump_events_for_ms(300);

    EXPECT_EQ(notebook->sessionTitle(), title);
    EXPECT_EQ(sessions->activeIndex(), sessions->titles().size() - 1);
    EXPECT_EQ(notebook->selectedIndex(), 1);
    EXPECT_EQ(plot->variable(), QString("u"));
    EXPECT_DOUBLE_EQ(plot->xMin(), -2.5);
    EXPECT_DOUBLE_EQ(plot->xMax(), 4.0);
    EXPECT_EQ(core->selectedInspectorTab(), QString("Steps"));
    EXPECT_EQ(core->selectedRepresentationId(), QString("text"));

    ASSERT_GE(notebook->cellCount(), 3);
    const auto* restored_define = notebook->rawCells().at(0);
    const auto* restored_simplify = notebook->rawCells().at(1);
    ASSERT_NE(restored_define, nullptr);
    ASSERT_NE(restored_simplify, nullptr);
    EXPECT_EQ(restored_define->inputLatex(), QString("a = 2"));
    EXPECT_EQ(restored_define->outputLatex(), QString("a = 2"));
    EXPECT_EQ(restored_simplify->inputLatex(), QString("x+0"));
    EXPECT_EQ(restored_simplify->outputLatex(), QString("x"));
    ASSERT_FALSE(restored_simplify->alternatives().isEmpty());
    ASSERT_FALSE(restored_simplify->steps().isEmpty());
    EXPECT_TRUE(restored_simplify->restoredFromWorkspace());

    notebook->setSelectedIndex(2);
    auto* restored_blank = notebook->currentCell();
    ASSERT_NE(restored_blank, nullptr);
    restored_blank->setInputLatex("a + 3");
    pump_events_for_ms(600);
    core->executeCurrentCell();
    pump_events_for_ms(1200);

    ASSERT_GE(notebook->cellCount(), 3);
    const auto* restored_eval = notebook->rawCells().at(2);
    ASSERT_NE(restored_eval, nullptr);
    EXPECT_EQ(restored_eval->outputLatex(), QString("5"));
}

TEST(GuiUndoRedo, ReloadThenPlotUsesRestoredWorkspaceState) {
    ASSERT_NE(QGuiApplication::instance(), nullptr);

    QTemporaryDir temp_dir;
    ASSERT_TRUE(temp_dir.isValid());
    const QString storage_path = temp_dir.filePath("workspace.json");
    QFile::remove(storage_path);
    qputenv("CAS_GUI_STORAGE_PATH", storage_path.toUtf8());

    AppCore* core = AppCore::instance();
    ASSERT_NE(core, nullptr);

    core->reloadWorkspaceFromDisk();
    pump_events_for_ms(500);

    auto* notebook = core->notebook();
    auto* plot = core->plot();
    ASSERT_NE(notebook, nullptr);
    ASSERT_NE(plot, nullptr);

    const QString title =
        QString("Plot Reload Session %1").arg(QDateTime::currentMSecsSinceEpoch());
    core->newSession();
    core->renameActiveSession(title);

    auto* plot_cell = notebook->currentCell();
    ASSERT_NE(plot_cell, nullptr);
    plot_cell->setInputLatex("x");
    pump_events_for_ms(600);

    plot->setVariable("x");
    plot->setXMin(-1.5);
    plot->setXMax(2.0);
    pump_events_for_ms(1200);

    plot_cell->setInputLatex("corrupt");
    plot->setVariable("corrupt");
    core->invokeCommand("reload_workspace");
    pump_events_for_ms(400);

    EXPECT_EQ(notebook->sessionTitle(), title);
    EXPECT_EQ(notebook->selectedIndex(), 0);
    ASSERT_GE(notebook->cellCount(), 1);
    EXPECT_EQ(notebook->currentCell()->inputLatex(), QString("x"));
    EXPECT_EQ(plot->variable(), QString("x"));
    EXPECT_DOUBLE_EQ(plot->xMin(), -1.5);
    EXPECT_DOUBLE_EQ(plot->xMax(), 2.0);
    EXPECT_EQ(plot->state(), QString("idle"));

    core->plotSelectedCell();
    pump_events_for_ms(1200);

    EXPECT_TRUE(plot->hasPoints());
    EXPECT_EQ(plot->state(), QString("ready"));
    EXPECT_TRUE(plot->status().startsWith("Plot ready"));
    ASSERT_FALSE(plot->points().isEmpty());

    const QPointF first_point = plot->points().front().toPointF();
    const QPointF last_point = plot->points().back().toPointF();
    EXPECT_LT(first_point.x(), last_point.x());
    EXPECT_LT(first_point.y(), last_point.y());
}

TEST(GuiUndoRedo, InspectorAndPlotControlsStayCoherent) {
    ASSERT_NE(QGuiApplication::instance(), nullptr);

    QTemporaryDir temp_dir;
    ASSERT_TRUE(temp_dir.isValid());
    const QString storage_path = temp_dir.filePath("workspace.json");
    QFile::remove(storage_path);
    qputenv("CAS_GUI_STORAGE_PATH", storage_path.toUtf8());

    AppCore* core = AppCore::instance();
    ASSERT_NE(core, nullptr);
    core->reloadWorkspaceFromDisk();
    pump_events_for_ms(300);

    auto* notebook = core->notebook();
    auto* plot = core->plot();
    ASSERT_NE(notebook, nullptr);
    ASSERT_NE(plot, nullptr);

    auto* current = notebook->currentCell();
    ASSERT_NE(current, nullptr);
    current->setInputLatex("x+0");
    pump_events_for_ms(600);
    core->executeCurrentCell();
    pump_events_for_ms(1200);

    ASSERT_FALSE(core->selectedCellAlternatives().isEmpty());
    EXPECT_FALSE(core->selectedRepresentationId().isEmpty());

    core->setSelectedRepresentationId("text");
    EXPECT_EQ(core->selectedRepresentationId(), QString("text"));

    core->setSelectedInspectorTab("Representations");
    EXPECT_EQ(core->selectedInspectorTab(), QString("Representations"));

    plot->applyPreset("positive_10");
    EXPECT_DOUBLE_EQ(plot->xMin(), 0.0);
    EXPECT_DOUBLE_EQ(plot->xMax(), 10.0);

    const double preset_span = plot->xMax() - plot->xMin();
    core->zoomPlotIn();
    EXPECT_LT(plot->xMax() - plot->xMin(), preset_span);

    const double before_pan_min = plot->xMin();
    const double before_pan_max = plot->xMax();
    core->panPlotRight();
    EXPECT_GT(plot->xMin(), before_pan_min);
    EXPECT_GT(plot->xMax(), before_pan_max);

    core->resetPlotRange();
    EXPECT_DOUBLE_EQ(plot->xMin(), -10.0);
    EXPECT_DOUBLE_EQ(plot->xMax(), 10.0);

    plot->setXMin(4.0);
    plot->setXMax(4.0);
    core->plotSelectedCell();
    pump_events_for_ms(200);
    EXPECT_EQ(plot->state(), QString("error"));
    EXPECT_FALSE(plot->hasPoints());

    core->resetPlotRange();
    core->plotSelectedCell();
    pump_events_for_ms(1200);
    EXPECT_EQ(plot->state(), QString("ready"));
    EXPECT_TRUE(plot->sampleCount() > 0);

    pump_events_for_ms(1200);
    const auto loaded = cas::gui::storage::SessionStore::load();
    ASSERT_TRUE(loaded.found);
    EXPECT_EQ(loaded.workspace.selectedInspectorTab, QString("Representations"));
    EXPECT_EQ(loaded.workspace.selectedRepresentationId, QString("text"));
}

TEST(GuiUndoRedo, CommandPaletteWorkflowStaysAlignedWithCommands) {
    AppCore* core = AppCore::instance();
    ASSERT_NE(core, nullptr);
    auto* palette = core->palette();
    auto* notebook = core->notebook();
    ASSERT_NE(palette, nullptr);
    ASSERT_NE(notebook, nullptr);

    bool saw_reload_shortcut = false;
    bool saw_delete_session_shortcut = false;
    bool saw_delete_cell_shortcut = false;
    bool saw_run_and_advance = false;
    bool saw_insert_above = false;
    for (const auto& item : palette->commands()) {
        const auto command = item.toMap();
        const auto id = command.value("id").toString();
        const auto shortcut = command.value("shortcut").toString();
        if (id == "reload_workspace" && shortcut == "⌘⇧R") {
            saw_reload_shortcut = true;
        } else if (id == "delete_session" && shortcut == "⌘⌫") {
            saw_delete_session_shortcut = true;
        } else if (id == "delete_cell" && shortcut == "⌘⌦") {
            saw_delete_cell_shortcut = true;
        } else if (id == "run_and_advance" && shortcut == "⌘↵") {
            saw_run_and_advance = true;
        } else if (id == "insert_cell_above" && shortcut == "⌘⌥↑") {
            saw_insert_above = true;
        }
    }
    EXPECT_TRUE(saw_reload_shortcut);
    EXPECT_TRUE(saw_delete_session_shortcut);
    EXPECT_TRUE(saw_delete_cell_shortcut);
    EXPECT_TRUE(saw_run_and_advance);
    EXPECT_TRUE(saw_insert_above);

    while (notebook->cellCount() > 1) {
        notebook->setSelectedIndex(notebook->cellCount() - 1);
        core->deleteCurrentCell();
    }
    notebook->setSelectedIndex(0);
    EXPECT_FALSE(core->isCommandEnabled("merge_with_previous"));
    EXPECT_FALSE(core->isCommandEnabled("select_previous_cell"));
    EXPECT_FALSE(core->isCommandEnabled("merge_with_next"));
    EXPECT_FALSE(core->isCommandEnabled("select_next_cell"));

    core->insertCellBelow();
    EXPECT_TRUE(core->isCommandEnabled("merge_with_previous"));
    EXPECT_TRUE(core->isCommandEnabled("select_previous_cell"));
    EXPECT_FALSE(core->isCommandEnabled("merge_with_next"));
    EXPECT_FALSE(core->isCommandEnabled("select_next_cell"));
    core->selectPreviousCell();
    EXPECT_TRUE(core->isCommandEnabled("merge_with_next"));
    EXPECT_TRUE(core->isCommandEnabled("select_next_cell"));

    palette->setQuery("ricarica");
    EXPECT_EQ(palette->firstMatchingCommandId(), QString("reload_workspace"));

    palette->setVisible(true);
    core->invokeCommand("reload_workspace");
    EXPECT_FALSE(palette->visible());
}

int main(int argc, char** argv) {
    qputenv("QT_QPA_PLATFORM", "offscreen");
    QGuiApplication app(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

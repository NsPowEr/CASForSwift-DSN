#include <gtest/gtest.h>

#include <QCoreApplication>
#include <QElapsedTimer>
#include <QEventLoop>
#include <QGuiApplication>
#include <QTemporaryDir>
#include <QThread>

#include "viewmodels/AppCore.h"

namespace {

void pump_events_for_ms(int duration_ms) {
    QElapsedTimer timer;
    timer.start();
    while (timer.elapsed() < duration_ms) {
        QCoreApplication::processEvents(QEventLoop::AllEvents, 20);
        QThread::msleep(10);
    }
}

void reset_notebook(AppCore* core) {
    ASSERT_NE(core, nullptr);
    core->reloadWorkspaceFromDisk();
    pump_events_for_ms(400);

    auto* notebook = core->notebook();
    ASSERT_NE(notebook, nullptr);
    while (notebook->cellCount() > 1) {
        notebook->setSelectedIndex(notebook->cellCount() - 1);
        core->deleteCurrentCell();
    }
    notebook->setSelectedIndex(0);
    core->focusCell(0);
}

} // namespace

TEST(GuiNotebookInteraction, BackendSeparatesFocusEditingAndSelectionTransitions) {
    ASSERT_NE(QGuiApplication::instance(), nullptr);

    QTemporaryDir temp_dir;
    ASSERT_TRUE(temp_dir.isValid());
    qputenv("CAS_GUI_STORAGE_PATH", temp_dir.filePath("workspace.json").toUtf8());

    AppCore* core = AppCore::instance();
    ASSERT_NE(core, nullptr);
    reset_notebook(core);

    auto* notebook = core->notebook();
    ASSERT_NE(notebook, nullptr);
    auto* first = notebook->currentCell();
    ASSERT_NE(first, nullptr);

    EXPECT_TRUE(first->active());
    EXPECT_TRUE(first->focused());
    EXPECT_FALSE(first->editing());
    EXPECT_EQ(core->selectedCellInteractionState(), QString("focused"));

    core->beginEditingCell(0);
    EXPECT_TRUE(first->active());
    EXPECT_TRUE(first->focused());
    EXPECT_TRUE(first->editing());
    EXPECT_EQ(core->selectedCellInteractionState(), QString("editing"));

    core->insertCellBelow();
    ASSERT_EQ(notebook->selectedIndex(), 1);
    auto* second = notebook->currentCell();
    ASSERT_NE(second, nullptr);
    EXPECT_FALSE(first->editing());
    EXPECT_TRUE(second->active());
    EXPECT_TRUE(second->focused());
    EXPECT_TRUE(second->editing());

    core->focusCell(0);
    EXPECT_EQ(notebook->selectedIndex(), 0);
    EXPECT_TRUE(first->active());
    EXPECT_TRUE(first->focused());
    EXPECT_FALSE(first->editing());
    EXPECT_FALSE(second->focused());
    EXPECT_FALSE(second->editing());
    EXPECT_EQ(core->selectedCellInteractionState(), QString("focused"));
}

TEST(GuiNotebookInteraction, RunAndAdvanceStartsEditingNextCell) {
    ASSERT_NE(QGuiApplication::instance(), nullptr);

    QTemporaryDir temp_dir;
    ASSERT_TRUE(temp_dir.isValid());
    qputenv("CAS_GUI_STORAGE_PATH", temp_dir.filePath("workspace.json").toUtf8());

    AppCore* core = AppCore::instance();
    ASSERT_NE(core, nullptr);
    reset_notebook(core);

    auto* notebook = core->notebook();
    ASSERT_NE(notebook, nullptr);
    auto* first = notebook->currentCell();
    ASSERT_NE(first, nullptr);

    first->setInputLatex("a = 2");
    pump_events_for_ms(600);

    core->beginEditingCell(0);
    core->runAndAdvance();
    pump_events_for_ms(1400);

    ASSERT_EQ(notebook->selectedIndex(), 1);
    ASSERT_GE(notebook->cellCount(), 2);
    auto* second = notebook->currentCell();
    ASSERT_NE(second, nullptr);

    EXPECT_EQ(first->outputLatex(), QString("a = 2"));
    EXPECT_FALSE(first->editing());
    EXPECT_TRUE(second->active());
    EXPECT_TRUE(second->focused());
    EXPECT_TRUE(second->editing());
    EXPECT_TRUE(second->inputLatex().isEmpty());
    EXPECT_EQ(core->selectedCellInteractionState(), QString("editing"));
}

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

void reset_to_single_cell(AppCore* core) {
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
    auto* cell = notebook->currentCell();
    ASSERT_NE(cell, nullptr);
    cell->setInputLatex(QString());
    pump_events_for_ms(600);
}

} // namespace

TEST(GuiCommands, RunPlotAndSplitRequireNonEmptyInput) {
    ASSERT_NE(QGuiApplication::instance(), nullptr);

    QTemporaryDir temp_dir;
    ASSERT_TRUE(temp_dir.isValid());
    qputenv("CAS_GUI_STORAGE_PATH", temp_dir.filePath("workspace.json").toUtf8());

    AppCore* core = AppCore::instance();
    ASSERT_NE(core, nullptr);
    reset_to_single_cell(core);

    EXPECT_FALSE(core->canRunCell());
    EXPECT_FALSE(core->isCommandEnabled("run"));
    EXPECT_FALSE(core->isCommandEnabled("run_and_advance"));
    EXPECT_FALSE(core->isCommandEnabled("plot"));
    EXPECT_FALSE(core->isCommandEnabled("split_cell"));

    auto* cell = core->notebook()->currentCell();
    ASSERT_NE(cell, nullptr);
    cell->setInputLatex("x+0");
    pump_events_for_ms(600);

    EXPECT_TRUE(core->canRunCell());
    EXPECT_TRUE(core->isCommandEnabled("run"));
    EXPECT_TRUE(core->isCommandEnabled("run_and_advance"));
    EXPECT_TRUE(core->isCommandEnabled("plot"));
    EXPECT_TRUE(core->isCommandEnabled("split_cell"));
}

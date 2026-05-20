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

void reset_palette_state(AppCore* core) {
    ASSERT_NE(core, nullptr);
    core->reloadWorkspaceFromDisk();
    pump_events_for_ms(400);
    auto* palette = core->palette();
    ASSERT_NE(palette, nullptr);
    palette->setVisible(false);
    palette->setQuery(QString());
}

} // namespace

TEST(GuiCommandPalette, RecentCommandsAndQueryRankingStayUseful) {
    ASSERT_NE(QGuiApplication::instance(), nullptr);

    QTemporaryDir temp_dir;
    ASSERT_TRUE(temp_dir.isValid());
    qputenv("CAS_GUI_STORAGE_PATH", temp_dir.filePath("workspace.json").toUtf8());

    AppCore* core = AppCore::instance();
    ASSERT_NE(core, nullptr);
    reset_palette_state(core);

    auto* palette = core->palette();
    ASSERT_NE(palette, nullptr);

    palette->recordInvocation("reload_workspace");
    const QVariantList initial = palette->filteredCommands();
    ASSERT_FALSE(initial.isEmpty());
    EXPECT_EQ(initial.front().toMap().value(QStringLiteral("id")).toString(),
              QString("reload_workspace"));

    palette->setQuery("inspector");
    EXPECT_EQ(palette->firstMatchingCommandId(), QString("toggle_inspector"));

    palette->setQuery("next");
    EXPECT_EQ(palette->firstMatchingCommandId(), QString("run_and_advance"));
}

TEST(GuiCommandPalette, DynamicLabelsAndDisabledReasonsStayBackendDriven) {
    ASSERT_NE(QGuiApplication::instance(), nullptr);

    QTemporaryDir temp_dir;
    ASSERT_TRUE(temp_dir.isValid());
    qputenv("CAS_GUI_STORAGE_PATH", temp_dir.filePath("workspace.json").toUtf8());

    AppCore* core = AppCore::instance();
    ASSERT_NE(core, nullptr);
    reset_palette_state(core);

    auto* notebook = core->notebook();
    ASSERT_NE(notebook, nullptr);
    auto* cell = notebook->currentCell();
    ASSERT_NE(cell, nullptr);
    cell->setInputLatex(QString());
    pump_events_for_ms(600);

    EXPECT_EQ(core->commandDisplayLabel("toggle_sidebar"), QString("Nascondi sidebar"));
    core->toggleSidebar();
    EXPECT_EQ(core->commandDisplayLabel("toggle_sidebar"), QString("Mostra sidebar"));
    core->toggleSidebar();

    EXPECT_FALSE(core->isCommandEnabled("run"));
    EXPECT_EQ(core->commandDisabledReason("run"), QString("Serve una cella con input non vuoto"));
}

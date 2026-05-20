#include <gtest/gtest.h>

#include <QCoreApplication>
#include <QElapsedTimer>
#include <QEventLoop>
#include <QFile>
#include <QGuiApplication>
#include <QTemporaryDir>
#include <QThread>

#include "core/storage/SessionStore.h"
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

} // namespace

TEST(GuiWorkspaceLayout, PersistAndRestoreLayoutState) {
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

    core->setSidebarVisible(false);
    core->setInspectorVisible(true);
    core->setSidebarWidth(275.0);
    core->setInspectorWidth(410.0);
    core->setNotebookScrollY(84.0);
    pump_events_for_ms(800);

    const auto loaded = cas::gui::storage::SessionStore::load();
    ASSERT_TRUE(loaded.found);
    EXPECT_FALSE(loaded.workspace.ui.sidebarVisible);
    EXPECT_TRUE(loaded.workspace.ui.inspectorVisible);
    EXPECT_DOUBLE_EQ(loaded.workspace.ui.sidebarWidth, 275.0);
    EXPECT_DOUBLE_EQ(loaded.workspace.ui.inspectorWidth, 410.0);
    EXPECT_DOUBLE_EQ(loaded.workspace.ui.notebookScrollY, 84.0);

    core->setSidebarVisible(true);
    core->setInspectorVisible(false);
    core->setSidebarWidth(230.0);
    core->setInspectorWidth(320.0);
    core->setNotebookScrollY(0.0);

    core->reloadWorkspaceFromDisk();
    pump_events_for_ms(300);

    EXPECT_FALSE(core->sidebarVisible());
    EXPECT_TRUE(core->inspectorVisible());
    EXPECT_DOUBLE_EQ(core->sidebarWidth(), 275.0);
    EXPECT_DOUBLE_EQ(core->inspectorWidth(), 410.0);
    EXPECT_DOUBLE_EQ(core->notebookScrollY(), 84.0);

    core->resetLayoutState();
    EXPECT_TRUE(core->sidebarVisible());
    EXPECT_TRUE(core->inspectorVisible());
    EXPECT_DOUBLE_EQ(core->sidebarWidth(), 230.0);
    EXPECT_DOUBLE_EQ(core->inspectorWidth(), 320.0);
    EXPECT_DOUBLE_EQ(core->notebookScrollY(), 0.0);
}

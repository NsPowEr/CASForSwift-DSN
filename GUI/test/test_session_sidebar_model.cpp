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

QVariantMap session_item(const QVariantList& items, int index) {
    if (index < 0 || index >= items.size()) {
        return {};
    }
    return items.at(index).toMap();
}

} // namespace

TEST(GuiSessions, SidebarModelExposesSessionMetadata) {
    ASSERT_NE(QGuiApplication::instance(), nullptr);

    QTemporaryDir temp_dir;
    ASSERT_TRUE(temp_dir.isValid());
    qputenv("CAS_GUI_STORAGE_PATH", temp_dir.filePath("workspace.json").toUtf8());

    AppCore* core = AppCore::instance();
    ASSERT_NE(core, nullptr);
    core->reloadWorkspaceFromDisk();
    pump_events_for_ms(400);

    auto* notebook = core->notebook();
    auto* sessions = core->sessions();
    ASSERT_NE(notebook, nullptr);
    ASSERT_NE(sessions, nullptr);

    notebook->currentCell()->setInputLatex("a = 2");
    pump_events_for_ms(600);
    core->executeCurrentCell();
    pump_events_for_ms(1200);
    core->addEmptyCell();
    notebook->currentCell()->setInputLatex("b + 1");
    pump_events_for_ms(1200);

    const QVariantList items = sessions->items();
    ASSERT_FALSE(items.isEmpty());
    const QVariantMap first = session_item(items, sessions->activeIndex());
    ASSERT_FALSE(first.isEmpty());
    EXPECT_EQ(first.value(QStringLiteral("title")).toString(), notebook->sessionTitle());
    EXPECT_EQ(first.value(QStringLiteral("cellCount")).toInt(), notebook->cellCount());
    EXPECT_EQ(first.value(QStringLiteral("definitionCount")).toInt(), 1);
    EXPECT_TRUE(first.value(QStringLiteral("summary")).toString().contains("cell"));
    EXPECT_EQ(first.value(QStringLiteral("selectedPreview")).toString(), QString("b + 1"));
    EXPECT_EQ(first.value(QStringLiteral("selectedStatusLabel")).toString(), QString("draft"));
    EXPECT_EQ(first.value(QStringLiteral("errorCount")).toInt(), 0);
}

TEST(GuiSessions, ReorderSessionKeepsActiveSessionContent) {
    ASSERT_NE(QGuiApplication::instance(), nullptr);

    QTemporaryDir temp_dir;
    ASSERT_TRUE(temp_dir.isValid());
    qputenv("CAS_GUI_STORAGE_PATH", temp_dir.filePath("workspace.json").toUtf8());

    AppCore* core = AppCore::instance();
    ASSERT_NE(core, nullptr);
    core->reloadWorkspaceFromDisk();
    pump_events_for_ms(400);

    core->renameActiveSession("Base Session");
    core->newSession();
    core->renameActiveSession("Second Session");
    core->notebook()->currentCell()->setInputLatex("x+0");
    pump_events_for_ms(600);

    ASSERT_EQ(core->sessions()->activeIndex(), 1);
    core->moveSessionUp(1);

    EXPECT_EQ(core->sessions()->activeIndex(), 0);
    EXPECT_EQ(core->notebook()->sessionTitle(), QString("Second Session"));
    EXPECT_EQ(core->notebook()->currentCell()->inputLatex(), QString("x+0"));

    const QVariantList items = core->sessions()->items();
    ASSERT_EQ(items.size(), 2);
    EXPECT_EQ(session_item(items, 0).value(QStringLiteral("title")).toString(), QString("Second Session"));
    EXPECT_TRUE(session_item(items, 0).value(QStringLiteral("active")).toBool());
    EXPECT_EQ(session_item(items, 1).value(QStringLiteral("title")).toString(), QString("Base Session"));

    core->moveSessionDown(0);
    EXPECT_EQ(core->sessions()->activeIndex(), 1);
    EXPECT_EQ(core->notebook()->sessionTitle(), QString("Second Session"));
}

TEST(GuiSessions, SidebarModelSummarizesErrorState) {
    ASSERT_NE(QGuiApplication::instance(), nullptr);

    QTemporaryDir temp_dir;
    ASSERT_TRUE(temp_dir.isValid());
    qputenv("CAS_GUI_STORAGE_PATH", temp_dir.filePath("workspace.json").toUtf8());

    AppCore* core = AppCore::instance();
    ASSERT_NE(core, nullptr);
    core->reloadWorkspaceFromDisk();
    pump_events_for_ms(400);

    core->notebook()->currentCell()->setInputLatex("sin(");
    pump_events_for_ms(600);
    core->executeCurrentCell();
    pump_events_for_ms(1200);

    const QVariantMap current = session_item(core->sessions()->items(), core->sessions()->activeIndex());
    ASSERT_FALSE(current.isEmpty());
    EXPECT_EQ(current.value(QStringLiteral("selectedStatusLabel")).toString(), QString("error"));
    EXPECT_GE(current.value(QStringLiteral("errorCount")).toInt(), 1);
    EXPECT_TRUE(current.value(QStringLiteral("summary")).toString().contains("err"));
}

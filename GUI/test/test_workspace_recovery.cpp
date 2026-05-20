#include <gtest/gtest.h>

#include <QDateTime>
#include <QDir>
#include <QFile>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QTemporaryDir>

#include "core/storage/SessionStore.h"
#include "viewmodels/AppCore.h"

namespace {

void write_workspace_file(const QString& path, const QByteArray& payload) {
    QFile file(path);
    ASSERT_TRUE(file.open(QIODevice::WriteOnly | QIODevice::Truncate));
    ASSERT_EQ(file.write(payload), payload.size());
    file.close();
}

QJsonObject make_cell(const QString& input) {
    return QJsonObject{
        {"input", input},
        {"status", "ready"},
        {"output", QString()},
        {"meta", QString()},
        {"has_output", false},
        {"alternatives", QJsonArray()},
        {"steps", QJsonArray()},
    };
}

} // namespace

TEST(GuiWorkspaceRecovery, CorruptWorkspaceRecoversWithBackup) {
    QTemporaryDir temp_dir;
    ASSERT_TRUE(temp_dir.isValid());
    const QString storage_path = temp_dir.filePath("workspace.json");
    qputenv("CAS_GUI_STORAGE_PATH", storage_path.toUtf8());
    write_workspace_file(storage_path, QByteArrayLiteral("{ definitely not valid json"));

    AppCore* core = AppCore::instance();
    ASSERT_NE(core, nullptr);

    core->reloadWorkspaceFromDisk();

    EXPECT_TRUE(core->storageStatus().startsWith("Recovered workspace after load failure:"));
    EXPECT_EQ(core->storageSeverity(), 3);
    EXPECT_FALSE(core->storageNotification().isEmpty());
    EXPECT_TRUE(core->storageNotification().contains("load failed"));

    core->clearStorageNotification();
    EXPECT_TRUE(core->storageNotification().isEmpty());

    ASSERT_NE(core->sessions(), nullptr);
    ASSERT_NE(core->notebook(), nullptr);
    EXPECT_EQ(core->sessions()->titles().size(), 1);
    EXPECT_EQ(core->notebook()->cellCount(), 1);
    EXPECT_EQ(core->notebook()->sessionTitle(), QString("Sessione 1"));

    const QDir dir(temp_dir.path());
    const QStringList backups =
        dir.entryList(QStringList() << "workspace.corrupt-*.json", QDir::Files);
    ASSERT_EQ(backups.size(), 1);
    EXPECT_FALSE(QFile::exists(storage_path));
}

TEST(GuiWorkspaceRecovery, InvalidPersistedWorkspaceStateIsNormalized) {
    QTemporaryDir temp_dir;
    ASSERT_TRUE(temp_dir.isValid());
    const QString storage_path = temp_dir.filePath("workspace.json");
    qputenv("CAS_GUI_STORAGE_PATH", storage_path.toUtf8());

    const QJsonObject root{
        {"active_session_index", 99},
        {"plot",
         QJsonObject{
             {"variable", ""},
             {"x_min", 5.0},
             {"x_max", 5.0},
         }},
        {"ui",
         QJsonObject{
             {"sidebar_visible", false},
             {"inspector_visible", true},
             {"sidebar_width", 999.0},
             {"inspector_width", 42.0},
             {"notebook_scroll_y", -5.0},
         }},
        {"selected_inspector_tab", ""},
        {"selected_representation_id", "text"},
        {"sessions",
         QJsonArray{
             QJsonObject{
                 {"title", ""},
                 {"selected_index", 99},
                 {"definitions", QJsonArray{"not-an-object"}},
                 {"cells", QJsonArray{make_cell("x+0"), QStringLiteral("skip-me")}},
             },
         }},
    };
    write_workspace_file(storage_path, QJsonDocument(root).toJson(QJsonDocument::Compact));

    AppCore* core = AppCore::instance();
    ASSERT_NE(core, nullptr);

    core->reloadWorkspaceFromDisk();

    ASSERT_NE(core->sessions(), nullptr);
    ASSERT_NE(core->notebook(), nullptr);
    ASSERT_NE(core->plot(), nullptr);
    EXPECT_EQ(core->sessions()->activeIndex(), 0);
    EXPECT_EQ(core->notebook()->sessionTitle(), QString("Sessione"));
    EXPECT_EQ(core->notebook()->cellCount(), 1);
    EXPECT_EQ(core->notebook()->currentCell()->inputLatex(), QString("x+0"));
    EXPECT_EQ(core->notebook()->selectedIndex(), 0);
    EXPECT_EQ(core->plot()->variable(), QString("x"));
    EXPECT_DOUBLE_EQ(core->plot()->xMin(), -10.0);
    EXPECT_DOUBLE_EQ(core->plot()->xMax(), 10.0);
    EXPECT_FALSE(core->sidebarVisible());
    EXPECT_TRUE(core->inspectorVisible());
    EXPECT_DOUBLE_EQ(core->sidebarWidth(), 230.0);
    EXPECT_DOUBLE_EQ(core->inspectorWidth(), 320.0);
    EXPECT_DOUBLE_EQ(core->notebookScrollY(), 0.0);
    EXPECT_EQ(core->selectedInspectorTab(), QString("Result"));
    EXPECT_TRUE(core->selectedRepresentationId().isEmpty());
    EXPECT_TRUE(core->storageStatus().startsWith("Workspace restored with recovery:"));
    EXPECT_EQ(core->storageSeverity(), 2);
    EXPECT_FALSE(core->storageNotification().isEmpty());
    EXPECT_TRUE(core->storageNotification().contains("restored with some fixes"));

    core->clearStorageNotification();
    EXPECT_TRUE(core->storageNotification().isEmpty());

    const auto loaded = cas::gui::storage::SessionStore::load();
    ASSERT_TRUE(loaded.found);
    EXPECT_TRUE(loaded.warning.contains("Recovered empty session title"));
    EXPECT_TRUE(loaded.warning.contains("Skipped non-object cell entry"));
    EXPECT_TRUE(loaded.warning.contains("Skipped non-object definition entry"));
    EXPECT_TRUE(loaded.warning.contains("Recovered empty plot variable"));
    EXPECT_TRUE(loaded.warning.contains("Recovered invalid plot range"));
    EXPECT_TRUE(loaded.warning.contains("Recovered invalid sidebar width"));
    EXPECT_TRUE(loaded.warning.contains("Recovered invalid inspector width"));
    EXPECT_TRUE(loaded.warning.contains("Recovered invalid notebook scroll position"));
    EXPECT_TRUE(loaded.warning.contains("Clamped invalid selected cell index"));
    EXPECT_TRUE(loaded.warning.contains("Clamped invalid active session index"));
}

#include <gtest/gtest.h>

#include <QDir>
#include <QFileInfo>
#include <QApplication>
#include <QCoreApplication>
#include <QEventLoop>
#include <QQmlApplicationEngine>
#include <QQmlComponent>
#include <QQmlError>
#include <QFont>
#include <QMetaObject>
#include <QStringList>

#include "viewmodels/AppCore.h"

namespace {

QStringList warning_strings(const QList<QQmlError>& warnings) {
    QStringList lines;
    lines.reserve(warnings.size());
    for (const auto& warning : warnings) {
        lines.push_back(warning.toString());
    }
    return lines;
}

void add_gui_import_paths(QQmlApplicationEngine& engine) {
    const QDir app_dir(QCoreApplication::applicationDirPath());
    const QStringList candidates = {
        app_dir.filePath("GUI"),
        app_dir.absolutePath(),
        QDir::current().filePath("build-gui-arm64-ninja/GUI"),
        QDir::current().filePath("build-gui/GUI"),
    };

    for (const auto& candidate : candidates) {
        engine.addImportPath(candidate);
    }
}

} // namespace

TEST(GuiQmlSmoke, MainModuleLoadsOffscreen) {
    ASSERT_NE(QApplication::instance(), nullptr);

    QQmlApplicationEngine engine;
    QList<QQmlError> warnings;
    QObject::connect(
        &engine,
        &QQmlApplicationEngine::warnings,
        [&warnings](const QList<QQmlError>& current) { warnings = current; });

    const QDir app_dir(QCoreApplication::applicationDirPath());
    ASSERT_TRUE(QFileInfo::exists(app_dir.filePath("GUI/CAS/qmldir")) ||
                QFileInfo::exists(QDir::current().filePath("build-gui-arm64-ninja/GUI/CAS/qmldir")));
    add_gui_import_paths(engine);
    engine.loadFromModule("CAS", "Main");

    ASSERT_TRUE(warnings.isEmpty()) << warning_strings(warnings).join('\n').toStdString();
    ASSERT_FALSE(engine.rootObjects().isEmpty());

    QObject* root = engine.rootObjects().front();
    ASSERT_NE(root, nullptr);
    EXPECT_EQ(root->objectName(), QString("mainWindow"));

    EXPECT_NE(root->findChild<QObject*>("titleBar"), nullptr);
    EXPECT_NE(root->findChild<QObject*>("sidebar"), nullptr);
    EXPECT_NE(root->findChild<QObject*>("notebook"), nullptr);
    EXPECT_NE(root->findChild<QObject*>("inspector"), nullptr);
    EXPECT_NE(root->findChild<QObject*>("commandPalette"), nullptr);
    EXPECT_NE(root->findChild<QObject*>("inspectorTabBar"), nullptr);
}

TEST(GuiQmlSmoke, MainWindowRegistersCriticalShortcuts) {
    ASSERT_NE(QApplication::instance(), nullptr);

    QQmlApplicationEngine engine;
    add_gui_import_paths(engine);
    engine.loadFromModule("CAS", "Main");
    ASSERT_FALSE(engine.rootObjects().isEmpty());

    QObject* root = engine.rootObjects().front();
    ASSERT_NE(root, nullptr);

    const auto* reload = root->findChild<QObject*>("shortcutReloadWorkspace");
    const auto* plot = root->findChild<QObject*>("shortcutPlot");
    const auto* delete_session = root->findChild<QObject*>("shortcutDeleteSession");
    const auto* delete_cell = root->findChild<QObject*>("shortcutDeleteCell");
    const auto* select_prev = root->findChild<QObject*>("shortcutSelectPrev");
    const auto* select_next = root->findChild<QObject*>("shortcutSelectNext");
    const auto* insert_above = root->findChild<QObject*>("shortcutInsertCellAbove");
    const auto* toggle_sidebar = root->findChild<QObject*>("shortcutToggleSidebar");
    const auto* toggle_inspector = root->findChild<QObject*>("shortcutToggleInspector");
    const auto* reset_layout = root->findChild<QObject*>("shortcutResetLayout");
    const auto* undo = root->findChild<QObject*>("shortcutUndo");
    const auto* redo = root->findChild<QObject*>("shortcutRedo");

    ASSERT_NE(reload, nullptr);
    ASSERT_NE(plot, nullptr);
    ASSERT_NE(delete_session, nullptr);
    ASSERT_NE(delete_cell, nullptr);
    ASSERT_NE(select_prev, nullptr);
    ASSERT_NE(select_next, nullptr);
    ASSERT_NE(insert_above, nullptr);
    ASSERT_NE(toggle_sidebar, nullptr);
    ASSERT_NE(toggle_inspector, nullptr);
    ASSERT_NE(reset_layout, nullptr);
    ASSERT_NE(undo, nullptr);
    ASSERT_NE(redo, nullptr);

    EXPECT_EQ(reload->property("sequence").toString(), QString("Meta+Shift+R"));
    EXPECT_EQ(plot->property("sequence").toString(), QString("Meta+Shift+P"));
    EXPECT_EQ(delete_session->property("sequence").toString(), QString("Meta+Backspace"));
    EXPECT_EQ(delete_cell->property("sequence").toString(), QString("Meta+Delete"));
    EXPECT_EQ(select_prev->property("sequence").toString(), QString("Ctrl+Up"));
    EXPECT_EQ(select_next->property("sequence").toString(), QString("Ctrl+Down"));
    EXPECT_EQ(insert_above->property("sequence").toString(), QString("Meta+Alt+Up"));
    EXPECT_EQ(toggle_sidebar->property("sequence").toString(), QString("Meta+Alt+S"));
    EXPECT_EQ(toggle_inspector->property("sequence").toString(), QString("Meta+Alt+I"));
    EXPECT_EQ(reset_layout->property("sequence").toString(), QString("Meta+Alt+0"));
    EXPECT_EQ(undo->property("sequence").toString(), QString("Meta+Z"));
    EXPECT_EQ(redo->property("sequence").toString(), QString("Meta+Shift+Z"));
}

TEST(GuiQmlSmoke, InspectorTabsFollowAppCoreState) {
    ASSERT_NE(QApplication::instance(), nullptr);

    QQmlApplicationEngine engine;
    add_gui_import_paths(engine);
    engine.loadFromModule("CAS", "Main");
    ASSERT_FALSE(engine.rootObjects().isEmpty());

    QObject* root = engine.rootObjects().front();
    ASSERT_NE(root, nullptr);

    auto* tab_bar = root->findChild<QObject*>("inspectorTabBar");
    auto* stack = root->findChild<QObject*>("inspectorStack");
    ASSERT_NE(tab_bar, nullptr);
    ASSERT_NE(stack, nullptr);

    auto* core = AppCore::instance();
    ASSERT_NE(core, nullptr);

    core->setSelectedInspectorTab("Plot");
    QCoreApplication::processEvents(QEventLoop::AllEvents, 50);
    EXPECT_EQ(tab_bar->property("currentIndex").toInt(), 4);
    EXPECT_EQ(stack->property("currentIndex").toInt(), 4);

    core->setSelectedInspectorTab("Steps");
    QCoreApplication::processEvents(QEventLoop::AllEvents, 50);
    EXPECT_EQ(tab_bar->property("currentIndex").toInt(), 2);
    EXPECT_EQ(stack->property("currentIndex").toInt(), 2);
}

TEST(GuiQmlSmoke, LayoutVisibilityFollowsAppCoreState) {
    ASSERT_NE(QApplication::instance(), nullptr);

    QQmlApplicationEngine engine;
    add_gui_import_paths(engine);
    engine.loadFromModule("CAS", "Main");
    ASSERT_FALSE(engine.rootObjects().isEmpty());

    QObject* root = engine.rootObjects().front();
    ASSERT_NE(root, nullptr);

    auto* sidebar = root->findChild<QObject*>("sidebar");
    auto* inspector = root->findChild<QObject*>("inspector");
    auto* split = root->findChild<QObject*>("workspaceSplit");
    ASSERT_NE(sidebar, nullptr);
    ASSERT_NE(inspector, nullptr);
    ASSERT_NE(split, nullptr);

    auto* core = AppCore::instance();
    ASSERT_NE(core, nullptr);

    core->setSidebarVisible(false);
    core->setInspectorVisible(false);
    QCoreApplication::processEvents(QEventLoop::AllEvents, 50);
    EXPECT_FALSE(sidebar->property("visible").toBool());
    EXPECT_FALSE(inspector->property("visible").toBool());

    core->setSidebarVisible(true);
    core->setInspectorVisible(true);
    core->setSidebarWidth(260.0);
    core->setInspectorWidth(390.0);
    QCoreApplication::processEvents(QEventLoop::AllEvents, 50);
    EXPECT_TRUE(sidebar->property("visible").toBool());
    EXPECT_TRUE(inspector->property("visible").toBool());
    EXPECT_EQ(split->objectName(), QString("workspaceSplit"));
}

TEST(GuiQmlSmoke, VariablesInspectorExposesFilterField) {
    ASSERT_NE(QApplication::instance(), nullptr);

    QQmlApplicationEngine engine;
    add_gui_import_paths(engine);
    engine.loadFromModule("CAS", "Main");
    ASSERT_FALSE(engine.rootObjects().isEmpty());

    QObject* root = engine.rootObjects().front();
    ASSERT_NE(root, nullptr);

    auto* search = root->findChild<QObject*>("variableSearchField");
    ASSERT_NE(search, nullptr);

    auto* core = AppCore::instance();
    ASSERT_NE(core, nullptr);

    core->setVariableSearchQuery("alpha");
    QCoreApplication::processEvents(QEventLoop::AllEvents, 50);
    EXPECT_EQ(search->property("text").toString(), QString("alpha"));

    core->setVariableSearchQuery(QString());
    QCoreApplication::processEvents(QEventLoop::AllEvents, 50);
    EXPECT_TRUE(search->property("text").toString().isEmpty());
}

TEST(GuiQmlSmoke, CellDelegateTracksBackendInteractionState) {
    ASSERT_NE(QApplication::instance(), nullptr);

    QQmlApplicationEngine engine;
    add_gui_import_paths(engine);

    auto* core = AppCore::instance();
    ASSERT_NE(core, nullptr);
    ASSERT_NE(core->notebook()->currentCell(), nullptr);

    const QStringList candidates = {
        QDir::current().filePath("GUI/macos-cpp/qml/workspace/CellView.qml"),
        QDir::current().filePath("../GUI/macos-cpp/qml/workspace/CellView.qml"),
        QDir::current().filePath("GUI/CAS/workspace/CellView.qml"),
        QDir(QCoreApplication::applicationDirPath()).filePath("GUI/macos-cpp/qml/workspace/CellView.qml"),
        QDir(QCoreApplication::applicationDirPath()).filePath("../GUI/macos-cpp/qml/workspace/CellView.qml"),
        QDir(QCoreApplication::applicationDirPath()).filePath("GUI/CAS/workspace/CellView.qml"),
    };

    QString cell_view_path;
    for (const auto& candidate : candidates) {
        if (QFileInfo::exists(candidate)) {
            cell_view_path = candidate;
            break;
        }
    }
    ASSERT_FALSE(cell_view_path.isEmpty());

    QQmlComponent component(&engine, QUrl::fromLocalFile(cell_view_path));
    ASSERT_TRUE(component.isReady()) << warning_strings(component.errors()).join('\n').toStdString();

    QObject* cell_view = component.createWithInitialProperties(
        {{"cell", QVariant::fromValue(core->notebook()->currentCell())}});
    ASSERT_NE(cell_view, nullptr);

    auto* input = cell_view->findChild<QObject*>("cellInput");
    ASSERT_NE(input, nullptr);

    core->focusCell(0);
    QCoreApplication::processEvents(QEventLoop::AllEvents, 50);
    EXPECT_EQ(core->selectedCellInteractionState(), QString("focused"));
    EXPECT_FALSE(input->property("focus").toBool());

    core->beginEditingCell(0);
    QCoreApplication::processEvents(QEventLoop::AllEvents, 100);
    EXPECT_EQ(core->selectedCellInteractionState(), QString("editing"));
    EXPECT_TRUE(input->property("focus").toBool());

    core->endEditingCell(0);
    QCoreApplication::processEvents(QEventLoop::AllEvents, 100);
    EXPECT_EQ(core->selectedCellInteractionState(), QString("focused"));
    EXPECT_FALSE(input->property("focus").toBool());

    delete cell_view;
}

int main(int argc, char** argv) {
    qputenv("QT_QPA_PLATFORM", "offscreen");
    QApplication app(argc, argv);
    app.setFont(QFont("Helvetica Neue", 13));
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

#include <gtest/gtest.h>

#include <QCoreApplication>
#include <QClipboard>
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

QVariantMap find_variable(const QVariantList& variables, const QString& name) {
    for (const auto& item : variables) {
        const QVariantMap map = item.toMap();
        if (map.value(QStringLiteral("name")).toString() == name) {
            return map;
        }
    }
    return {};
}

void prepare_variable_session(AppCore* core, const QString& definition) {
    ASSERT_NE(core, nullptr);
    core->reloadWorkspaceFromDisk();
    pump_events_for_ms(400);

    auto* notebook = core->notebook();
    ASSERT_NE(notebook, nullptr);

    while (notebook->cellCount() > 1) {
        notebook->setSelectedIndex(notebook->cellCount() - 1);
        core->deleteCurrentCell();
    }

    auto* cell = notebook->currentCell();
    ASSERT_NE(cell, nullptr);
    cell->setInputLatex(definition);
    pump_events_for_ms(600);
    core->runAndAdvance();
    pump_events_for_ms(1200);
}

} // namespace

TEST(GuiVariables, FilterAndInsertActionsStayOperational) {
    ASSERT_NE(QGuiApplication::instance(), nullptr);

    QTemporaryDir temp_dir;
    ASSERT_TRUE(temp_dir.isValid());
    qputenv("CAS_GUI_STORAGE_PATH", temp_dir.filePath("workspace.json").toUtf8());

    AppCore* core = AppCore::instance();
    ASSERT_NE(core, nullptr);
    prepare_variable_session(core, "alphaVar = 42");

    const QVariantMap variable = find_variable(core->variables(), QStringLiteral("alphaVar"));
    ASSERT_FALSE(variable.isEmpty());
    EXPECT_EQ(variable.value(QStringLiteral("value")).toString(), QString("42"));

    core->setVariableSearchQuery("alpha");
    const QVariantList filtered_by_name = core->filteredVariables();
    ASSERT_EQ(filtered_by_name.size(), 1);
    EXPECT_EQ(filtered_by_name.front().toMap().value(QStringLiteral("name")).toString(),
              QString("alphaVar"));

    core->setVariableSearchQuery("42");
    const QVariantList filtered_by_value = core->filteredVariables();
    ASSERT_EQ(filtered_by_value.size(), 1);
    EXPECT_EQ(filtered_by_value.front().toMap().value(QStringLiteral("name")).toString(),
              QString("alphaVar"));

    auto* notebook = core->notebook();
    ASSERT_NE(notebook, nullptr);
    auto* current = notebook->currentCell();
    ASSERT_NE(current, nullptr);
    EXPECT_TRUE(current->inputLatex().isEmpty());

    core->insertVariableName("alphaVar");
    EXPECT_EQ(current->inputLatex(), QString("alphaVar"));

    core->insertVariableValue("alphaVar");
    EXPECT_EQ(current->inputLatex(), QString("alphaVar 42"));
}

TEST(GuiVariables, ClipboardActionsCopyLiveVariableData) {
    ASSERT_NE(QGuiApplication::instance(), nullptr);

    QTemporaryDir temp_dir;
    ASSERT_TRUE(temp_dir.isValid());
    qputenv("CAS_GUI_STORAGE_PATH", temp_dir.filePath("workspace.json").toUtf8());

    AppCore* core = AppCore::instance();
    ASSERT_NE(core, nullptr);
    prepare_variable_session(core, "betaVar = x+1");
    const QVariantMap variable = find_variable(core->variables(), QStringLiteral("betaVar"));
    ASSERT_FALSE(variable.isEmpty());

    auto* clipboard = QGuiApplication::clipboard();
    ASSERT_NE(clipboard, nullptr);

    core->copyVariableName("betaVar");
    EXPECT_EQ(clipboard->text(), QString("betaVar"));

    core->copyVariableValue("betaVar");
    EXPECT_EQ(clipboard->text(), variable.value(QStringLiteral("value")).toString());
}

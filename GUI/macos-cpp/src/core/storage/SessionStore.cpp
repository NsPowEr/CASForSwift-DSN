#include "SessionStore.h"

#include <QDateTime>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QSaveFile>
#include <QStandardPaths>
#include <QtGlobal>

namespace cas::gui::storage {
namespace {

QJsonObject to_json(const StoredCell& cell) {
    return QJsonObject{
        {"input", cell.input},
        {"status", cell.status},
        {"output", cell.output},
        {"meta", cell.meta},
        {"has_output", cell.hasOutput},
        {"alternatives", QJsonArray::fromVariantList(cell.alternatives)},
        {"steps", QJsonArray::fromVariantList(cell.steps)},
    };
}

StoredCell cell_from_json(const QJsonObject& object) {
    StoredCell cell;
    cell.input = object.value("input").toString();
    cell.status = object.value("status").toString("ready");
    cell.output = object.value("output").toString();
    cell.meta = object.value("meta").toString();
    cell.hasOutput = object.value("has_output").toBool(false);
    cell.alternatives = object.value("alternatives").toArray().toVariantList();
    cell.steps = object.value("steps").toArray().toVariantList();
    return cell;
}

QJsonObject to_json(const StoredDefinition& definition) {
    return QJsonObject{
        {"name", definition.name},
        {"value_text", definition.valueText},
    };
}

StoredDefinition definition_from_json(const QJsonObject& object) {
    StoredDefinition definition;
    definition.name = object.value("name").toString();
    definition.valueText = object.value("value_text").toString();
    return definition;
}

QJsonObject to_json(const StoredSession& session) {
    QJsonArray cells;
    for (const auto& cell : session.cells) {
        cells.push_back(to_json(cell));
    }
    QJsonArray definitions;
    for (const auto& definition : session.definitions) {
        definitions.push_back(to_json(definition));
    }

    return QJsonObject{
        {"title", session.title},
        {"selected_index", session.selectedIndex},
        {"definitions", definitions},
        {"cells", cells},
    };
}

StoredSession session_from_json(const QJsonObject& object, QStringList* warnings) {
    StoredSession session;
    session.title = object.value("title").toString("Sessione");
    session.selectedIndex = object.value("selected_index").toInt(0);

    if (session.title.trimmed().isEmpty()) {
        session.title = QStringLiteral("Sessione");
        if (warnings) {
            warnings->push_back(QStringLiteral("Recovered empty session title"));
        }
    }

    const QJsonArray cells = object.value("cells").toArray();
    session.cells.reserve(static_cast<std::size_t>(cells.size()));
    for (const auto& value : cells) {
        if (value.isObject()) {
            session.cells.push_back(cell_from_json(value.toObject()));
        } else if (warnings) {
            warnings->push_back(QStringLiteral("Skipped non-object cell entry"));
        }
    }
    const QJsonArray definitions = object.value("definitions").toArray();
    session.definitions.reserve(static_cast<std::size_t>(definitions.size()));
    for (const auto& value : definitions) {
        if (value.isObject()) {
            session.definitions.push_back(definition_from_json(value.toObject()));
        } else if (warnings) {
            warnings->push_back(QStringLiteral("Skipped non-object definition entry"));
        }
    }

    if (session.cells.empty()) {
        session.cells.push_back(StoredCell{});
        session.selectedIndex = 0;
        if (warnings) {
            warnings->push_back(QStringLiteral("Recovered empty session cells"));
        }
    } else if (session.selectedIndex < 0 ||
               session.selectedIndex >= static_cast<int>(session.cells.size())) {
        session.selectedIndex =
            qBound(0, session.selectedIndex, static_cast<int>(session.cells.size()) - 1);
        if (warnings) {
            warnings->push_back(QStringLiteral("Clamped invalid selected cell index"));
        }
    }

    return session;
}

QJsonObject to_json(const StoredPlotSettings& plot) {
    return QJsonObject{
        {"variable", plot.variable},
        {"x_min", plot.xMin},
        {"x_max", plot.xMax},
    };
}

QJsonObject to_json(const StoredUiState& ui) {
    return QJsonObject{
        {"sidebar_visible", ui.sidebarVisible},
        {"inspector_visible", ui.inspectorVisible},
        {"sidebar_width", ui.sidebarWidth},
        {"inspector_width", ui.inspectorWidth},
        {"notebook_scroll_y", ui.notebookScrollY},
    };
}

StoredPlotSettings plot_from_json(const QJsonObject& object, QStringList* warnings) {
    StoredPlotSettings plot;
    plot.variable = object.value("variable").toString("x");
    plot.xMin = object.value("x_min").toDouble(-10.0);
    plot.xMax = object.value("x_max").toDouble(10.0);

    if (plot.variable.trimmed().isEmpty()) {
        plot.variable = QStringLiteral("x");
        if (warnings) {
            warnings->push_back(QStringLiteral("Recovered empty plot variable"));
        }
    }
    if (!qIsFinite(plot.xMin) || !qIsFinite(plot.xMax) || !(plot.xMin < plot.xMax)) {
        plot.xMin = -10.0;
        plot.xMax = 10.0;
        if (warnings) {
            warnings->push_back(QStringLiteral("Recovered invalid plot range"));
        }
    }

    return plot;
}

StoredUiState ui_from_json(const QJsonObject& object, QStringList* warnings) {
    StoredUiState ui;
    ui.sidebarVisible = object.value("sidebar_visible").toBool(true);
    ui.inspectorVisible = object.value("inspector_visible").toBool(true);
    ui.sidebarWidth = object.value("sidebar_width").toDouble(230.0);
    ui.inspectorWidth = object.value("inspector_width").toDouble(320.0);
    ui.notebookScrollY = object.value("notebook_scroll_y").toDouble(0.0);

    if (!qIsFinite(ui.sidebarWidth) || ui.sidebarWidth < 180.0 || ui.sidebarWidth > 420.0) {
        ui.sidebarWidth = 230.0;
        if (warnings) {
            warnings->push_back(QStringLiteral("Recovered invalid sidebar width"));
        }
    }
    if (!qIsFinite(ui.inspectorWidth) || ui.inspectorWidth < 260.0 || ui.inspectorWidth > 520.0) {
        ui.inspectorWidth = 320.0;
        if (warnings) {
            warnings->push_back(QStringLiteral("Recovered invalid inspector width"));
        }
    }
    if (!qIsFinite(ui.notebookScrollY) || ui.notebookScrollY < 0.0) {
        ui.notebookScrollY = 0.0;
        if (warnings) {
            warnings->push_back(QStringLiteral("Recovered invalid notebook scroll position"));
        }
    }

    return ui;
}

QString backup_corrupt_workspace(const QString& path) {
    const QFileInfo info(path);
    if (!info.exists()) {
        return {};
    }

    const QString suffix = info.suffix().isEmpty() ? QStringLiteral("json") : info.suffix();
    const QString backup_path =
        info.dir().filePath(QString("%1.corrupt-%2.%3")
                                .arg(info.completeBaseName())
                                .arg(QDateTime::currentDateTimeUtc().toString("yyyyMMdd-hhmmsszzz"))
                                .arg(suffix));
    QFile::remove(backup_path);
    if (QFile::rename(path, backup_path)) {
        return backup_path;
    }
    return {};
}

} // namespace

QString SessionStore::storagePath() {
    const QString override_path = qEnvironmentVariable("CAS_GUI_STORAGE_PATH");
    if (!override_path.trimmed().isEmpty()) {
        return override_path;
    }

    const QString app_data_dir =
        QStandardPaths::writableLocation(QStandardPaths::AppDataLocation);
    return app_data_dir + "/workspace.json";
}

LoadResult SessionStore::load() {
    LoadResult result;
    QStringList warnings;
    const QString path = storagePath();
    QFile file(path);
    if (!file.exists()) {
        return result;
    }

    result.found = true;
    if (!file.open(QIODevice::ReadOnly)) {
        result.error = QString("Cannot open workspace: %1").arg(file.errorString());
        return result;
    }

    const QByteArray payload = file.readAll();
    file.close();
    QJsonParseError parse_error;
    const QJsonDocument doc = QJsonDocument::fromJson(payload, &parse_error);
    if (parse_error.error != QJsonParseError::NoError || !doc.isObject()) {
        const QString backup_path = backup_corrupt_workspace(path);
        result.error = QString("Invalid workspace JSON: %1").arg(parse_error.errorString());
        if (!backup_path.isEmpty()) {
            result.error += QString(" (backup: %1)").arg(QFileInfo(backup_path).fileName());
        }
        return result;
    }

    const QJsonObject root = doc.object();
    result.workspace.activeSessionIndex = root.value("active_session_index").toInt(0);
    if (root.value("plot").isObject()) {
        result.workspace.plot = plot_from_json(root.value("plot").toObject(), &warnings);
    }
    if (root.value("ui").isObject()) {
        result.workspace.ui = ui_from_json(root.value("ui").toObject(), &warnings);
    }
    result.workspace.selectedInspectorTab =
        root.value("selected_inspector_tab").toString("Result");
    result.workspace.selectedRepresentationId =
        root.value("selected_representation_id").toString();

    const QJsonArray sessions = root.value("sessions").toArray();
    result.workspace.sessions.reserve(static_cast<std::size_t>(sessions.size()));
    for (const auto& value : sessions) {
        if (value.isObject()) {
            result.workspace.sessions.push_back(session_from_json(value.toObject(), &warnings));
        } else {
            warnings.push_back(QStringLiteral("Skipped non-object session entry"));
        }
    }

    if (result.workspace.sessions.empty()) {
        result.found = false;
        result.warning = warnings.join(QStringLiteral("; "));
        return result;
    }

    if (result.workspace.activeSessionIndex < 0 ||
        result.workspace.activeSessionIndex >= static_cast<int>(result.workspace.sessions.size())) {
        result.workspace.activeSessionIndex =
            qBound(0,
                   result.workspace.activeSessionIndex,
                   static_cast<int>(result.workspace.sessions.size()) - 1);
        warnings.push_back(QStringLiteral("Clamped invalid active session index"));
    }

    result.warning = warnings.join(QStringLiteral("; "));
    return result;
}

bool SessionStore::save(const StoredWorkspace& workspace, QString* error) {
    const QString path = storagePath();
    const QFileInfo info(path);
    QDir dir = info.dir();
    if (!dir.exists() && !dir.mkpath(".")) {
        if (error) {
            *error = QString("Cannot create storage directory: %1").arg(dir.path());
        }
        return false;
    }

    QJsonArray sessions;
    for (const auto& session : workspace.sessions) {
        sessions.push_back(to_json(session));
    }

    const QJsonObject root{
        {"schema_version", 1},
        {"active_session_index", workspace.activeSessionIndex},
        {"plot", to_json(workspace.plot)},
        {"ui", to_json(workspace.ui)},
        {"selected_inspector_tab", workspace.selectedInspectorTab},
        {"selected_representation_id", workspace.selectedRepresentationId},
        {"sessions", sessions},
    };

    QSaveFile file(path);
    if (!file.open(QIODevice::WriteOnly)) {
        if (error) {
            *error = QString("Cannot write workspace: %1").arg(file.errorString());
        }
        return false;
    }

    const QByteArray payload = QJsonDocument(root).toJson(QJsonDocument::Indented);
    if (file.write(payload) != payload.size()) {
        if (error) {
            *error = QString("Short write while saving workspace");
        }
        file.cancelWriting();
        return false;
    }
    if (!file.commit()) {
        if (error) {
            *error = QString("Cannot commit workspace: %1").arg(file.errorString());
        }
        return false;
    }
    return true;
}

} // namespace cas::gui::storage

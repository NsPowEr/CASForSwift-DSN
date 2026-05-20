#pragma once

#include <QString>
#include <QVariantList>
#include <vector>

namespace cas::gui::storage {

struct StoredCell {
    QString input;
    QString status;
    QString output;
    QString meta;
    bool hasOutput = false;
    QVariantList alternatives;
    QVariantList steps;
};

struct StoredDefinition {
    QString name;
    QString valueText;
};

struct StoredSession {
    QString title;
    std::vector<StoredCell> cells;
    std::vector<StoredDefinition> definitions;
    int selectedIndex = 0;
};

struct StoredPlotSettings {
    QString variable = "x";
    double xMin = -10.0;
    double xMax = 10.0;
};

struct StoredUiState {
    bool sidebarVisible = true;
    bool inspectorVisible = true;
    double sidebarWidth = 230.0;
    double inspectorWidth = 320.0;
    double notebookScrollY = 0.0;
};

struct StoredWorkspace {
    std::vector<StoredSession> sessions;
    int activeSessionIndex = 0;
    StoredPlotSettings plot;
    StoredUiState ui;
    QString selectedInspectorTab = "Result";
    QString selectedRepresentationId;
};

struct LoadResult {
    bool found = false;
    StoredWorkspace workspace;
    QString error;
    QString warning;
};

class SessionStore final {
public:
    static QString storagePath();
    static LoadResult load();
    static bool save(const StoredWorkspace& workspace, QString* error = nullptr);
};

} // namespace cas::gui::storage

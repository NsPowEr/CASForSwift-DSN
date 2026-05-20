// AppCore.h — root C++ singleton, esposto a QML
// Espone tutti i ViewModel come Q_PROPERTY per binding dichiarativo.

#pragma once
#include <QObject>
#include <QQmlEngine>
#include <QTimer>
#include <memory>
#include <vector>

#include "ThemeVM.h"
#include "SessionListVM.h"
#include "NotebookVM.h"
#include "KeyboardVM.h"
#include "CommandPaletteVM.h"
#include "PlotVM.h"
#include "../core/kernel/KernelDispatcher.h"
#include "../core/storage/SessionStore.h"

using namespace cas::gui;

class AppCore : public QObject {
    Q_OBJECT
    QML_ELEMENT
    QML_SINGLETON

    Q_PROPERTY(ThemeVM*           theme           READ theme           CONSTANT)
    Q_PROPERTY(SessionListVM*     sessions        READ sessions        CONSTANT)
    Q_PROPERTY(NotebookVM*        notebook        READ notebook        CONSTANT)
    Q_PROPERTY(KeyboardVM*        keyboard        READ keyboard        CONSTANT)
    Q_PROPERTY(CommandPaletteVM*  palette         READ palette         CONSTANT)
    Q_PROPERTY(PlotVM*            plot            READ plot            CONSTANT)
    Q_PROPERTY(QString            kernelVersion   READ kernelVersion   NOTIFY kernelChanged)
    Q_PROPERTY(QString            kernelMode      READ kernelMode      NOTIFY kernelChanged)
    Q_PROPERTY(QString            selectedCellInput READ selectedCellInput NOTIFY selectedCellChanged)
    Q_PROPERTY(QString            selectedCellStatus READ selectedCellStatus NOTIFY selectedCellChanged)
    Q_PROPERTY(QString            selectedCellOutput READ selectedCellOutput NOTIFY selectedCellChanged)
    Q_PROPERTY(QString            selectedCellMeta READ selectedCellMeta NOTIFY selectedCellChanged)
    Q_PROPERTY(QString            selectedCellState READ selectedCellState NOTIFY selectedCellChanged)
    Q_PROPERTY(QString            selectedCellInteractionState READ selectedCellInteractionState NOTIFY selectedCellChanged)
    Q_PROPERTY(bool               selectedCellHasError READ selectedCellHasError NOTIFY selectedCellChanged)
    Q_PROPERTY(bool               selectedCellRestored READ selectedCellRestored NOTIFY selectedCellChanged)
    Q_PROPERTY(bool               selectedCellHasOutput READ selectedCellHasOutput NOTIFY selectedCellChanged)
    Q_PROPERTY(QVariantList       selectedCellAlternatives READ selectedCellAlternatives NOTIFY selectedCellChanged)
    Q_PROPERTY(QVariantList       selectedCellSteps READ selectedCellSteps NOTIFY selectedCellChanged)
    Q_PROPERTY(QStringList        inspectorTabs READ inspectorTabs CONSTANT)
    Q_PROPERTY(QString            selectedInspectorTab READ selectedInspectorTab WRITE setSelectedInspectorTab NOTIFY selectedInspectorTabChanged)
    Q_PROPERTY(QString            selectedRepresentationId READ selectedRepresentationId WRITE setSelectedRepresentationId NOTIFY selectedRepresentationChanged)
    Q_PROPERTY(bool               canUndo READ canUndo NOTIFY undoAvailabilityChanged)
    Q_PROPERTY(bool               canRedo READ canRedo NOTIFY undoAvailabilityChanged)
    Q_PROPERTY(bool               canDeleteSession READ canDeleteSession NOTIFY commandStatesChanged)
    Q_PROPERTY(bool               canMoveSessionUp READ canMoveSessionUp NOTIFY commandStatesChanged)
    Q_PROPERTY(bool               canMoveSessionDown READ canMoveSessionDown NOTIFY commandStatesChanged)
    Q_PROPERTY(bool               canDeleteCell READ canDeleteCell NOTIFY commandStatesChanged)
    Q_PROPERTY(bool               canMoveCellUp READ canMoveCellUp NOTIFY commandStatesChanged)
    Q_PROPERTY(bool               canMoveCellDown READ canMoveCellDown NOTIFY commandStatesChanged)
    Q_PROPERTY(bool               canRunCell READ canRunCell NOTIFY commandStatesChanged)
    Q_PROPERTY(QString            undoActionLabel READ undoActionLabel NOTIFY undoAvailabilityChanged)
    Q_PROPERTY(QString            redoActionLabel READ redoActionLabel NOTIFY undoAvailabilityChanged)
    Q_PROPERTY(QString            storageStatus READ storageStatus NOTIFY storageStatusChanged)
    Q_PROPERTY(int                storageSeverity READ storageSeverity NOTIFY storageStatusChanged)
    Q_PROPERTY(QString            storageNotification READ storageNotification NOTIFY storageNotificationChanged)
    Q_PROPERTY(bool               isKernelBusy READ isKernelBusy NOTIFY kernelBusyChanged)
    Q_PROPERTY(bool               sidebarVisible READ sidebarVisible WRITE setSidebarVisible NOTIFY layoutStateChanged)
    Q_PROPERTY(bool               inspectorVisible READ inspectorVisible WRITE setInspectorVisible NOTIFY layoutStateChanged)
    Q_PROPERTY(qreal              sidebarWidth READ sidebarWidth WRITE setSidebarWidth NOTIFY layoutStateChanged)
    Q_PROPERTY(qreal              inspectorWidth READ inspectorWidth WRITE setInspectorWidth NOTIFY layoutStateChanged)
    Q_PROPERTY(qreal              notebookScrollY READ notebookScrollY WRITE setNotebookScrollY NOTIFY layoutStateChanged)
    Q_PROPERTY(QString            helpContent READ helpContent NOTIFY helpChanged)
    Q_PROPERTY(QVariantList       variables READ variables NOTIFY variablesChanged)
    Q_PROPERTY(QString            variableSearchQuery READ variableSearchQuery WRITE setVariableSearchQuery NOTIFY variableSearchQueryChanged)
    Q_PROPERTY(QVariantList       filteredVariables READ filteredVariables NOTIFY variablesChanged)
    Q_PROPERTY(QStringList        functions READ functions CONSTANT)

public:
    static AppCore* instance();
    static AppCore* create(QQmlEngine* engine, QJSEngine*) {
        AppCore* core = instance();
        if (engine) {
            QQmlEngine::setObjectOwnership(core, QQmlEngine::CppOwnership);
        }
        return core;
    }

    ThemeVM*          theme()    const { return m_theme.get(); }
    SessionListVM*    sessions() const { return m_sessions.get(); }
    NotebookVM*       notebook() const { return m_notebook.get(); }
    KeyboardVM*       keyboard() const { return m_keyboard.get(); }
    CommandPaletteVM* palette()  const { return m_palette.get(); }
    PlotVM*           plot()     const { return m_plot.get(); }

    QString kernelVersion() const;
    QString kernelMode()    const;  // "embedded" | "remote" | "hybrid"
    bool isKernelBusy() const { return m_kernelBusy; }
    QString selectedCellInput() const;
    QString selectedCellStatus() const;
    QString selectedCellOutput() const;
    QString selectedCellMeta() const;
    QString selectedCellState() const;
    QString selectedCellInteractionState() const;
    bool selectedCellHasError() const;
    bool selectedCellRestored() const;
    bool selectedCellHasOutput() const;
    QVariantList selectedCellAlternatives() const;
    QVariantList selectedCellSteps() const;
    QStringList inspectorTabs() const {
        return {QStringLiteral("Result"),
                QStringLiteral("Representations"),
                QStringLiteral("Steps"),
                QStringLiteral("Variables"),
                QStringLiteral("Plot")};
    }
    QString selectedInspectorTab() const { return m_selectedInspectorTab; }
    void setSelectedInspectorTab(const QString& value);
    QString selectedRepresentationId() const { return m_selectedRepresentationId; }
    void setSelectedRepresentationId(const QString& value);
    void setSidebarVisible(bool visible);
    void setInspectorVisible(bool visible);
    void setSidebarWidth(qreal width);
    void setInspectorWidth(qreal width);
    void setNotebookScrollY(qreal scrollY);
    bool canUndo() const { return !m_undoStack.empty(); }
    bool canRedo() const { return !m_redoStack.empty(); }
    bool canDeleteSession() const { return m_sessionStorage.size() > 1; }
    bool canMoveSessionUp() const;
    bool canMoveSessionDown() const;
    bool canDeleteCell() const;
    bool canMoveCellUp() const;
    bool canMoveCellDown() const;
    bool canRunCell() const;
    QString undoActionLabel() const;
    QString redoActionLabel() const;
    QString storageStatus() const { return m_storageStatus; }
    int storageSeverity() const { return m_storageSeverity; }
    QString storageNotification() const { return m_storageNotification; }
    bool sidebarVisible() const { return m_sidebarVisible; }
    bool inspectorVisible() const { return m_inspectorVisible; }
    qreal sidebarWidth() const { return m_sidebarWidth; }
    qreal inspectorWidth() const { return m_inspectorWidth; }
    qreal notebookScrollY() const { return m_notebookScrollY; }
    QString helpContent() const { return m_helpContent; }
    QVariantList variables() const { return m_kernel ? m_kernel->listVariables() : QVariantList(); }
    QString variableSearchQuery() const { return m_variableSearchQuery; }
    void setVariableSearchQuery(const QString& value);
    QVariantList filteredVariables() const;
    QStringList functions() const { return m_kernel ? m_kernel->listFunctions() : QStringList(); }

    Q_INVOKABLE void clearStorageNotification() { setStorageNotification(QString()); }
    Q_INVOKABLE void toggleSidebar();
    Q_INVOKABLE void toggleInspector();
    Q_INVOKABLE void resetLayoutState();
    Q_INVOKABLE void newSession();
    Q_INVOKABLE void openSession(int index);
    Q_INVOKABLE void renameActiveSession(const QString& title);
    Q_INVOKABLE void duplicateSession(int index = -1);
    Q_INVOKABLE void deleteSession(int index = -1);
    Q_INVOKABLE void moveSessionUp(int index = -1);
    Q_INVOKABLE void moveSessionDown(int index = -1);
    Q_INVOKABLE void addEmptyCell();
    Q_INVOKABLE void insertCellAbove();
    Q_INVOKABLE void insertCellBelow();
    Q_INVOKABLE void duplicateCurrentCell();
    Q_INVOKABLE void deleteCurrentCell();
    Q_INVOKABLE void moveCurrentCellUp();
    Q_INVOKABLE void moveCurrentCellDown();
    Q_INVOKABLE void selectPreviousCell();
    Q_INVOKABLE void selectNextCell();
    Q_INVOKABLE void focusCell(int index = -1);
    Q_INVOKABLE void beginEditingCell(int index = -1);
    Q_INVOKABLE void endEditingCell(int index = -1);
    Q_INVOKABLE void undo();
    Q_INVOKABLE void redo();
    Q_INVOKABLE void executeCurrentCell();
    Q_INVOKABLE void runAndAdvance();
    Q_INVOKABLE void splitCurrentCell(int cursorPosition = -1);
    Q_INVOKABLE void mergeWithPreviousCell();
    Q_INVOKABLE void mergeWithNextCell();
    Q_INVOKABLE void plotSelectedCell();
    Q_INVOKABLE void zoomPlotIn();
    Q_INVOKABLE void zoomPlotOut();
    Q_INVOKABLE void panPlotLeft();
    Q_INVOKABLE void panPlotRight();
    Q_INVOKABLE void resetPlotRange();
    Q_INVOKABLE void applyPlotPreset(const QString& presetId);
    Q_INVOKABLE bool isCommandEnabled(const QString& commandId) const;
    Q_INVOKABLE QString commandDisplayLabel(const QString& commandId) const;
    Q_INVOKABLE QString commandDisabledReason(const QString& commandId) const;
    Q_INVOKABLE void invokeCommand(const QString& commandId);
    Q_INVOKABLE void exportToHtml(const QString& filePath);
    Q_INVOKABLE void interruptKernel();
    Q_INVOKABLE void lookupHelp(const QString& query);
    Q_INVOKABLE void reloadWorkspaceFromDisk();
    Q_INVOKABLE void clearWorkspace();
    Q_INVOKABLE void quit();
    Q_INVOKABLE void copyVariableName(const QString& name);
    Q_INVOKABLE void copyVariableValue(const QString& name);
    Q_INVOKABLE void insertVariableName(const QString& name);
    Q_INVOKABLE void insertVariableValue(const QString& name);

    // Exposed for testing
    void beginUndoTransaction(const QString& label);
    void endUndoTransaction(bool commit = true);
    bool isUndoTransactionActive() const { return m_undoTransactionDepth > 0; }

signals:
    void kernelChanged();
    void selectedCellChanged();
    void storageStatusChanged();
    void helpChanged();
    void undoAvailabilityChanged();
    void commandStatesChanged();
    void variablesChanged();
    void selectedInspectorTabChanged();
    void selectedRepresentationChanged();
    void storageNotificationChanged();
    void kernelBusyChanged();
    void layoutStateChanged();
    void variableSearchQueryChanged();

private:
    explicit AppCore(QObject* parent = nullptr);
    CellVM* currentCell() const;
    CellVM* addTrackedCell();
    void attachCell(CellVM* cell);
    void initializeDefaultSession();
    void captureCurrentSession();
    void loadSession(int index);
    void restoreWorkspace();
    void restoreWorkspaceSnapshot(const storage::StoredWorkspace& workspace);
    void schedulePersist(const QString& pendingStatus = QString("Autosave pending"), int severity = 1);
    void persistWorkspace();
    void syncSessionTitles();
    void setStorageStatus(const QString& status, int severity = 0);
    void setStorageNotification(const QString& message);
    QString ensureUniqueSessionTitle(const QString& baseTitle, int skipIndex = -1) const;
    storage::StoredWorkspace snapshotWorkspace() const;
    bool workspaceEquals(const storage::StoredWorkspace& lhs, const storage::StoredWorkspace& rhs) const;
    void finalizeUndoableChange(const QString& label, const storage::StoredWorkspace& before);
    void pushUndoStep(const QString& label,
                      const storage::StoredWorkspace& before,
                      const storage::StoredWorkspace& after);
    storage::StoredWorkspace snapshotBeforeUndoableAction();
    void markPendingEditUndo(const QString& label);
    void commitPendingEditUndo();
    void ensureSelectedRepresentationStillValid();
    QVariantMap variableByName(const QString& name) const;
    bool insertTextIntoCurrentCell(const QString& text);
    void syncCellInteractionState();
    void setFocusedCellIndex(int index);
    void setEditingCellIndex(int index);

    struct UndoStep {
        QString label;
        storage::StoredWorkspace before;
        storage::StoredWorkspace after;
    };

    std::unique_ptr<KernelDispatcher>  m_kernel;
    std::unique_ptr<ThemeVM>           m_theme;
    std::unique_ptr<SessionListVM>     m_sessions;
    std::unique_ptr<NotebookVM>        m_notebook;
    std::unique_ptr<KeyboardVM>        m_keyboard;
    std::unique_ptr<CommandPaletteVM>  m_palette;
    std::unique_ptr<PlotVM>            m_plot;
    std::vector<storage::StoredSession> m_sessionStorage;
    QTimer                             m_autosaveTimer;
    QTimer                             m_editUndoTimer;
    QString                            m_storageStatus = "Storage idle";
    int                                m_storageSeverity = 0;
    QString                            m_storageNotification;
    bool                               m_kernelBusy = false;
    bool                               m_sidebarVisible = true;
    bool                               m_inspectorVisible = true;
    qreal                              m_sidebarWidth = 230.0;
    qreal                              m_inspectorWidth = 320.0;
    qreal                              m_notebookScrollY = 0.0;
    QString                            m_helpContent;
    QString                            m_variableSearchQuery;
    QString                            m_selectedInspectorTab = "Result";
    QString                            m_selectedRepresentationId;
    int                                m_focusedCellIndex = -1;
    int                                m_editingCellIndex = -1;
    bool                               m_restoringWorkspace = false;
    std::vector<UndoStep>              m_undoStack;
    std::vector<UndoStep>              m_redoStack;
    bool                               m_hasPendingEditUndo = false;
    QString                            m_pendingEditUndoLabel = "Edit cell";
    storage::StoredWorkspace           m_pendingEditBefore;
    int                                m_undoTransactionDepth = 0;
    bool                               m_undoTransactionCommit = true;
    QString                            m_undoTransactionLabel;
    storage::StoredWorkspace           m_undoTransactionBefore;
};

#include "AppCore.h"

#include <QCoreApplication>

AppCore* AppCore::instance() {
    static AppCore inst;
    return &inst;
}

AppCore::AppCore(QObject* parent) : QObject(parent) {
    m_kernel = std::make_unique<KernelDispatcher>(this);
    m_theme = std::make_unique<ThemeVM>(this);
    m_sessions = std::make_unique<SessionListVM>(this);
    m_notebook = std::make_unique<NotebookVM>(this);
    m_keyboard = std::make_unique<KeyboardVM>(this);
    m_palette = std::make_unique<CommandPaletteVM>(this);
    m_plot = std::make_unique<PlotVM>(this);
    m_autosaveTimer.setSingleShot(true);
    m_autosaveTimer.setInterval(350);
    connect(&m_autosaveTimer, &QTimer::timeout, this, &AppCore::persistWorkspace);
    m_editUndoTimer.setSingleShot(true);
    m_editUndoTimer.setInterval(800);
    connect(&m_editUndoTimer, &QTimer::timeout, this, &AppCore::commitPendingEditUndo);
    connect(m_notebook.get(), &NotebookVM::selectedIndexChanged, this, [this]() {
        const int selected_index = m_notebook ? m_notebook->selectedIndex() : -1;
        if (selected_index < 0) {
            m_focusedCellIndex = -1;
            m_editingCellIndex = -1;
        } else {
            if (m_focusedCellIndex < 0 || m_focusedCellIndex >= m_notebook->cellCount()) {
                m_focusedCellIndex = selected_index;
            }
            if (m_editingCellIndex >= 0 && m_editingCellIndex != selected_index) {
                m_editingCellIndex = -1;
            }
        }
        syncCellInteractionState();
        ensureSelectedRepresentationStillValid();
        emit selectedCellChanged();
        emit commandStatesChanged();
        if (!m_restoringWorkspace && !isUndoTransactionActive()) {
            captureCurrentSession();
            schedulePersist();
        }
    });
    connect(m_notebook.get(), &NotebookVM::cellsChanged, this, &AppCore::commandStatesChanged);
    connect(m_notebook.get(), &NotebookVM::sessionTitleChanged, this, [this]() {
        if (!m_restoringWorkspace && !isUndoTransactionActive()) {
            captureCurrentSession();
            syncSessionTitles();
            schedulePersist();
        }
    });
    connect(m_plot.get(), &PlotVM::settingsAboutToChange, this, [this]() {
        if (!m_restoringWorkspace && !isUndoTransactionActive()) {
            markPendingEditUndo("Edit plot");
        }
    });
    connect(m_plot.get(), &PlotVM::settingsChanged, this, [this]() {
        if (!m_restoringWorkspace && !isUndoTransactionActive()) {
            schedulePersist();
        }
    });
    connect(qApp, &QCoreApplication::aboutToQuit, this, &AppCore::persistWorkspace);
    restoreWorkspace();
}

QString AppCore::kernelVersion() const {
    return m_kernel ? m_kernel->version() : "unavailable";
}

QString AppCore::kernelMode() const {
    return m_kernel ? m_kernel->mode() : "unavailable";
}

CellVM* AppCore::currentCell() const {
    return m_notebook ? m_notebook->currentCell() : nullptr;
}

bool AppCore::canDeleteCell() const {
    return m_notebook && m_notebook->cellCount() > 1;
}

bool AppCore::canMoveSessionUp() const {
    return m_sessions && m_sessions->activeIndex() > 0;
}

bool AppCore::canMoveSessionDown() const {
    return m_sessions &&
           m_sessions->activeIndex() >= 0 &&
           m_sessions->activeIndex() < (static_cast<int>(m_sessionStorage.size()) - 1);
}

bool AppCore::canMoveCellUp() const {
    return m_notebook && m_notebook->selectedIndex() > 0;
}

bool AppCore::canMoveCellDown() const {
    return m_notebook &&
           m_notebook->selectedIndex() >= 0 &&
           m_notebook->selectedIndex() < (m_notebook->cellCount() - 1);
}

bool AppCore::canRunCell() const {
    auto* cell = currentCell();
    return cell && !cell->inputLatex().trimmed().isEmpty();
}

QString AppCore::selectedCellInteractionState() const {
    auto* cell = currentCell();
    return cell ? cell->interactionState() : QStringLiteral("idle");
}

QString AppCore::undoActionLabel() const {
    if (m_undoStack.empty()) return QStringLiteral("Undo");
    return QStringLiteral("Undo ") + m_undoStack.back().label;
}

QString AppCore::redoActionLabel() const {
    if (m_redoStack.empty()) return QStringLiteral("Redo");
    return QStringLiteral("Redo ") + m_redoStack.back().label;
}

bool AppCore::isCommandEnabled(const QString& commandId) const {
    if (commandId == QStringLiteral("undo") ||
        commandId == QStringLiteral("undo_delete_cell") ||
        commandId == QStringLiteral("undo_delete_session")) {
        return canUndo();
    }
    if (commandId == QStringLiteral("redo")) {
        return canRedo();
    }
    if (commandId == QStringLiteral("delete_session")) {
        return canDeleteSession();
    }
    if (commandId == QStringLiteral("delete_cell")) {
        return canDeleteCell();
    }
    if (commandId == QStringLiteral("move_cell_up")) {
        return canMoveCellUp();
    }
    if (commandId == QStringLiteral("move_cell_down")) {
        return canMoveCellDown();
    }
    if (commandId == QStringLiteral("select_previous_cell") ||
        commandId == QStringLiteral("merge_with_previous")) {
        return m_notebook && m_notebook->selectedIndex() > 0;
    }
    if (commandId == QStringLiteral("select_next_cell") ||
        commandId == QStringLiteral("merge_with_next")) {
        return m_notebook &&
               m_notebook->selectedIndex() >= 0 &&
               m_notebook->selectedIndex() < (m_notebook->cellCount() - 1);
    }
    if (commandId == QStringLiteral("toggle_sidebar") ||
        commandId == QStringLiteral("toggle_inspector") ||
        commandId == QStringLiteral("reset_layout")) {
        return true;
    }
    if (commandId == QStringLiteral("duplicate_cell")) {
        return currentCell() != nullptr;
    }
    if (commandId == QStringLiteral("run") ||
        commandId == QStringLiteral("run_and_advance") ||
        commandId == QStringLiteral("plot")) {
        return canRunCell();
    }
    if (commandId == QStringLiteral("insert_cell_above") ||
        commandId == QStringLiteral("insert_cell_below")) {
        return currentCell() != nullptr;
    }
    if (commandId == QStringLiteral("split_cell")) {
        return canRunCell();
    }
    return true;
}

QString AppCore::commandDisplayLabel(const QString& commandId) const {
    if (commandId == QStringLiteral("undo")) {
        return undoActionLabel();
    }
    if (commandId == QStringLiteral("redo")) {
        return redoActionLabel();
    }
    if (commandId == QStringLiteral("toggle_sidebar")) {
        return sidebarVisible() ? QStringLiteral("Nascondi sidebar")
                                : QStringLiteral("Mostra sidebar");
    }
    if (commandId == QStringLiteral("toggle_inspector")) {
        return inspectorVisible() ? QStringLiteral("Nascondi inspector")
                                  : QStringLiteral("Mostra inspector");
    }
    if (commandId == QStringLiteral("run")) {
        return QStringLiteral("Esegui cella");
    }
    if (commandId == QStringLiteral("run_and_advance")) {
        return QStringLiteral("Esegui e vai alla prossima");
    }
    if (commandId == QStringLiteral("plot")) {
        return QStringLiteral("Plot cella");
    }
    if (commandId == QStringLiteral("reload_workspace")) {
        return QStringLiteral("Ricarica workspace da disco");
    }
    if (commandId == QStringLiteral("new_session")) {
        return QStringLiteral("Nuova sessione");
    }
    if (commandId == QStringLiteral("duplicate_session")) {
        return QStringLiteral("Duplica sessione");
    }
    if (commandId == QStringLiteral("delete_session")) {
        return QStringLiteral("Elimina sessione");
    }
    if (commandId == QStringLiteral("insert_cell_above")) {
        return QStringLiteral("Inserisci cella sopra");
    }
    if (commandId == QStringLiteral("insert_cell_below")) {
        return QStringLiteral("Inserisci cella sotto");
    }
    if (commandId == QStringLiteral("duplicate_cell")) {
        return QStringLiteral("Duplica cella");
    }
    if (commandId == QStringLiteral("delete_cell")) {
        return QStringLiteral("Elimina cella");
    }
    if (commandId == QStringLiteral("split_cell")) {
        return QStringLiteral("Dividi cella");
    }
    if (commandId == QStringLiteral("merge_with_previous")) {
        return QStringLiteral("Unisci con la precedente");
    }
    if (commandId == QStringLiteral("merge_with_next")) {
        return QStringLiteral("Unisci con la successiva");
    }
    if (commandId == QStringLiteral("select_previous_cell")) {
        return QStringLiteral("Seleziona cella precedente");
    }
    if (commandId == QStringLiteral("select_next_cell")) {
        return QStringLiteral("Seleziona cella successiva");
    }
    if (commandId == QStringLiteral("move_cell_up")) {
        return QStringLiteral("Sposta cella su");
    }
    if (commandId == QStringLiteral("move_cell_down")) {
        return QStringLiteral("Sposta cella giu");
    }
    if (commandId == QStringLiteral("reset_layout")) {
        return QStringLiteral("Ripristina layout workspace");
    }
    return commandId;
}

QString AppCore::commandDisabledReason(const QString& commandId) const {
    if (isCommandEnabled(commandId)) {
        return {};
    }

    if (commandId == QStringLiteral("run") ||
        commandId == QStringLiteral("run_and_advance") ||
        commandId == QStringLiteral("plot") ||
        commandId == QStringLiteral("split_cell")) {
        return QStringLiteral("Serve una cella con input non vuoto");
    }
    if (commandId == QStringLiteral("merge_with_previous") ||
        commandId == QStringLiteral("select_previous_cell") ||
        commandId == QStringLiteral("move_cell_up")) {
        return QStringLiteral("Seleziona una cella che non sia la prima");
    }
    if (commandId == QStringLiteral("merge_with_next") ||
        commandId == QStringLiteral("select_next_cell") ||
        commandId == QStringLiteral("move_cell_down")) {
        return QStringLiteral("Seleziona una cella che non sia l'ultima");
    }
    if (commandId == QStringLiteral("delete_cell")) {
        return QStringLiteral("L'ultima cella non puo essere eliminata");
    }
    if (commandId == QStringLiteral("delete_session")) {
        return QStringLiteral("L'ultima sessione non puo essere eliminata");
    }
    if (commandId == QStringLiteral("undo") ||
        commandId == QStringLiteral("undo_delete_cell") ||
        commandId == QStringLiteral("undo_delete_session")) {
        return QStringLiteral("Nessuna azione da annullare");
    }
    if (commandId == QStringLiteral("redo")) {
        return QStringLiteral("Nessuna azione da ripetere");
    }
    return QStringLiteral("Comando non disponibile nel contesto corrente");
}

#include <QGuiApplication>
#include <QDateTime>
#include <QLoggingCategory>

Q_LOGGING_CATEGORY(casGuiAuditLog, "cas.gui.audit", QtWarningMsg)

void AppCore::undo() {
    commitPendingEditUndo();
    if (m_undoStack.empty()) {
        setStorageStatus("Nothing to undo", 1);
        return;
    }


    const UndoStep step = m_undoStack.back();
    m_undoStack.pop_back();
    
    qCDebug(casGuiAuditLog).noquote()
        << "[Audit] [" << QDateTime::currentDateTime().toString(Qt::ISODateWithMs)
        << "] undo:" << step.label;

    restoreWorkspaceSnapshot(step.before);
    m_redoStack.push_back(step);
    emit undoAvailabilityChanged();
    setStorageStatus(QString("Undo: %1").arg(step.label));
    schedulePersist();
}

void AppCore::redo() {
    commitPendingEditUndo();
    if (m_redoStack.empty()) {
        setStorageStatus("Nothing to redo", 1);
        return;
    }


    const UndoStep step = m_redoStack.back();
    m_redoStack.pop_back();
    
    qCDebug(casGuiAuditLog).noquote()
        << "[Audit] [" << QDateTime::currentDateTime().toString(Qt::ISODateWithMs)
        << "] redo:" << step.label;

    restoreWorkspaceSnapshot(step.after);
    m_undoStack.push_back(step);
    emit undoAvailabilityChanged();
    setStorageStatus(QString("Redo: %1").arg(step.label));
    schedulePersist();
}

#include <QFile>
#include <QTextStream>

void AppCore::exportToHtml(const QString& filePath) {
    QFile file(filePath);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        setStorageStatus("Export failed: cannot open file");
        return;
    }

    captureCurrentSession();

    QTextStream out(&file);
    out << "<!DOCTYPE html>\n<html>\n<head>\n";
    out << "<meta charset=\"UTF-8\">\n";
    out << "<title>CAS Workspace Export</title>\n";
    out << "<script src=\"https://polyfill.io/v3/polyfill.min.js?features=es6\"></script>\n";
    out << "<script id=\"MathJax-script\" async src=\"https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js\"></script>\n";
    out << "<style>\n";
    out << "body { font-family: -apple-system, sans-serif; max-width: 800px; margin: 40px auto; line-height: 1.6; color: #333; }\n";
    out << ".session { border-bottom: 2px solid #eee; padding-bottom: 20px; margin-bottom: 40px; }\n";
    out << ".cell { margin-bottom: 20px; padding: 15px; background: #f9f9f9; border-radius: 8px; }\n";
    out << ".input { font-family: monospace; color: #555; margin-bottom: 10px; border-left: 3px solid #ccc; padding-left: 10px; }\n";
    out << ".output { font-size: 1.2em; }\n";
    out << "h1 { color: #007aff; }\n";
    out << "</style>\n</head>\n<body>\n";
    out << "<h1>CAS Workspace</h1>\n";

    for (const auto& session : m_sessionStorage) {
        out << "<div class=\"session\">\n";
        out << "<h2>Session: " << session.title << "</h2>\n";
        for (const auto& cell : session.cells) {
            out << "<div class=\"cell\">\n";
            out << "<div class=\"input\">In: " << cell.input << "</div>\n";
            if (!cell.output.isEmpty()) {
                out << "<div class=\"output\">$$\n" << cell.output << "\n$$</div>\n";
            }
            out << "</div>\n";
        }
        out << "</div>\n";
    }

    out << "</body>\n</html>\n";
    file.close();
    setStorageStatus("Export complete: " + filePath);
}

void AppCore::lookupHelp(const QString& query) {
    static const QMap<QString, QString> helpDb = {
        {"sin", "<b>sin(x)</b>: Calcola il seno di x (in radianti)."},
        {"cos", "<b>cos(x)</b>: Calcola il coseno di x (in radianti)."},
        {"tan", "<b>tan(x)</b>: Calcola la tangente di x."},
        {"integrate", "<b>integrate(f, x)</b>: Calcola l'integrale indefinito di f rispetto a x.<br><b>integrate(f, x, a, b)</b>: Calcola l'integrale definito di f da a a b."},
        {"diff", "<b>diff(f, x)</b>: Calcola la derivata di f rispetto a x.<br><b>diff(f, x, n)</b>: Calcola la derivata n-esima."},
        {"solve", "<b>solve(eq, x)</b>: Risolve l'equazione eq rispetto a x."},
        {"limit", "<b>limit(f, x, a)</b>: Calcola il limite di f per x che tende ad a."},
        {"expand", "<b>expand(expr)</b>: Espande i prodotti e le potenze nell'espressione."},
        {"simplify", "<b>simplify(expr)</b>: Tenta di semplificare l'espressione."}
    };

    QString key = query.trimmed().toLower();
    if (key.startsWith('?')) {
        key = key.mid(1).trimmed();
    }

    if (helpDb.contains(key)) {
        m_helpContent = helpDb[key];
    } else {
        m_helpContent = "Nessun aiuto trovato per: " + key;
    }
    emit helpChanged();
}

void AppCore::interruptKernel() {
    if (m_kernel) {
        m_kernel->interrupt();
        setStorageStatus("Interrupting calculation...", 1);
    }
}

void AppCore::reloadWorkspaceFromDisk() {
    commitPendingEditUndo();
    restoreWorkspace();
}

void AppCore::clearWorkspace() {
    commitPendingEditUndo();
    const auto before = snapshotBeforeUndoableAction();
    
    m_restoringWorkspace = true;
    m_sessionStorage.clear();
    initializeDefaultSession();
    if (m_plot) {
        m_plot->clear();
    }
    m_focusedCellIndex = -1;
    m_editingCellIndex = -1;
    syncCellInteractionState();
    resetLayoutState();
    m_selectedInspectorTab = QStringLiteral("Result");
    m_selectedRepresentationId.clear();
    emit selectedInspectorTabChanged();
    emit selectedRepresentationChanged();
    
    m_restoringWorkspace = false;
    finalizeUndoableChange("Clear workspace", before);
    setStorageStatus("Workspace cleared");
}

void AppCore::quit() {
    QGuiApplication::quit();
}

void AppCore::setFocusedCellIndex(int index) {
    if (!m_notebook) {
        m_focusedCellIndex = -1;
        return;
    }

    if (index < 0 || index >= m_notebook->cellCount()) {
        m_focusedCellIndex = -1;
        return;
    }
    m_focusedCellIndex = index;
}

void AppCore::setEditingCellIndex(int index) {
    if (!m_notebook) {
        m_editingCellIndex = -1;
        return;
    }

    if (index < 0 || index >= m_notebook->cellCount()) {
        m_editingCellIndex = -1;
        return;
    }
    m_editingCellIndex = index;
    if (m_editingCellIndex >= 0) {
        m_focusedCellIndex = m_editingCellIndex;
    }
}

void AppCore::syncCellInteractionState() {
    if (!m_notebook) {
        return;
    }

    const int cell_count = m_notebook->cellCount();
    if (m_focusedCellIndex >= cell_count) {
        m_focusedCellIndex = -1;
    }
    if (m_editingCellIndex >= cell_count) {
        m_editingCellIndex = -1;
    }

    const auto& cells = m_notebook->rawCells();
    for (int i = 0; i < cells.size(); ++i) {
        auto* cell = cells.at(i);
        cell->setFocused(i == m_focusedCellIndex);
        cell->setEditing(i == m_editingCellIndex);
    }
}

void AppCore::focusCell(int index) {
    if (!m_notebook || m_notebook->cellCount() == 0) {
        return;
    }

    const int target = index >= 0 ? index : m_notebook->selectedIndex();
    if (target < 0 || target >= m_notebook->cellCount()) {
        return;
    }

    if (m_notebook->selectedIndex() != target) {
        m_notebook->setSelectedIndex(target);
    }
    setFocusedCellIndex(target);
    setEditingCellIndex(-1);
    syncCellInteractionState();
    emit selectedCellChanged();
}

void AppCore::beginEditingCell(int index) {
    if (!m_notebook || m_notebook->cellCount() == 0) {
        return;
    }

    const int target = index >= 0 ? index : m_notebook->selectedIndex();
    if (target < 0 || target >= m_notebook->cellCount()) {
        return;
    }

    if (m_notebook->selectedIndex() != target) {
        m_notebook->setSelectedIndex(target);
    }
    setEditingCellIndex(target);
    syncCellInteractionState();
    emit selectedCellChanged();
}

void AppCore::endEditingCell(int index) {
    if (!m_notebook || m_notebook->cellCount() == 0) {
        return;
    }

    const int target = index >= 0 ? index : m_notebook->selectedIndex();
    if (target < 0 || target >= m_notebook->cellCount()) {
        return;
    }

    if (m_editingCellIndex != target) {
        return;
    }

    setEditingCellIndex(-1);
    setFocusedCellIndex(target);
    syncCellInteractionState();
    emit selectedCellChanged();
}

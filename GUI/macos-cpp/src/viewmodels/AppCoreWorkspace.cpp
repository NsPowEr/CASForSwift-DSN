#include "AppCore.h"

#include <algorithm>
#include <QDateTime>

namespace {

QString summarize_session_preview(const QString& input) {
    const QString normalized = input.simplified();
    if (normalized.isEmpty()) {
        return QStringLiteral("Empty cell");
    }
    constexpr int kMaxPreviewLength = 56;
    if (normalized.size() <= kMaxPreviewLength) {
        return normalized;
    }
    return normalized.left(kMaxPreviewLength - 1) + QStringLiteral("…");
}

QString selected_cell_state_label(const storage::StoredCell& cell) {
    if (!cell.hasOutput) {
        return cell.input.trimmed().isEmpty() ? QStringLiteral("empty")
                                              : QStringLiteral("draft");
    }
    return cell.status == QStringLiteral("ok") ? QStringLiteral("ok")
                                               : QStringLiteral("error");
}

QVariantMap build_session_sidebar_item(const storage::StoredSession& session, int index, int activeIndex) {
    QVariantMap item;
    item[QStringLiteral("title")] = session.title;
    item[QStringLiteral("index")] = index;
    item[QStringLiteral("active")] = index == activeIndex;
    item[QStringLiteral("cellCount")] = static_cast<int>(session.cells.size());
    item[QStringLiteral("definitionCount")] = static_cast<int>(session.definitions.size());
    const int selected_index =
        std::clamp(session.selectedIndex, 0, static_cast<int>(session.cells.size()) - 1);
    item[QStringLiteral("selectedCellNumber")] = selected_index + 1;

    int outputCount = 0;
    int errorCount = 0;
    for (const auto& cell : session.cells) {
        if (cell.hasOutput) {
            ++outputCount;
            if (cell.status != QStringLiteral("ok")) {
                ++errorCount;
            }
        }
    }
    item[QStringLiteral("outputCount")] = outputCount;
    item[QStringLiteral("errorCount")] = errorCount;

    if (!session.cells.empty()) {
        const auto& selected = session.cells[static_cast<std::size_t>(selected_index)];
        item[QStringLiteral("selectedPreview")] = summarize_session_preview(selected.input);
        item[QStringLiteral("selectedStatusLabel")] = selected_cell_state_label(selected);
    } else {
        item[QStringLiteral("selectedPreview")] = QStringLiteral("Empty cell");
        item[QStringLiteral("selectedStatusLabel")] = QStringLiteral("empty");
    }

    QStringList summaryParts;
    summaryParts.push_back(QString("%1 cell%2")
                               .arg(static_cast<int>(session.cells.size()))
                               .arg(session.cells.size() == 1 ? "" : "e"));
    if (!session.definitions.empty()) {
        summaryParts.push_back(QString("%1 var").arg(static_cast<int>(session.definitions.size())));
    }
    if (outputCount > 0) {
        summaryParts.push_back(QString("%1 out").arg(outputCount));
    }
    if (errorCount > 0) {
        summaryParts.push_back(QString("%1 err").arg(errorCount));
    }
    item[QStringLiteral("summary")] = summaryParts.join(QStringLiteral(" · "));
    return item;
}

}

void AppCore::initializeDefaultSession() {
    storage::StoredSession session;
    session.title = QStringLiteral("Sessione 1");
    session.cells.push_back(storage::StoredCell{});
    m_sessionStorage.push_back(std::move(session));
    m_selectedInspectorTab = QStringLiteral("Result");
    m_selectedRepresentationId.clear();
    m_focusedCellIndex = -1;
    m_editingCellIndex = -1;
    syncSessionTitles();
    m_sessions->setActiveIndex(0);
    loadSession(0);
    syncSessionTitles();
}

void AppCore::syncSessionTitles() {
    QStringList titles;
    QVariantList items;
    titles.reserve(static_cast<qsizetype>(m_sessionStorage.size()));
    items.reserve(static_cast<qsizetype>(m_sessionStorage.size()));
    const int active_index = m_sessions ? m_sessions->activeIndex() : -1;
    int index = 0;
    for (const auto& session : m_sessionStorage) {
        titles.push_back(session.title);
        items.push_back(build_session_sidebar_item(session, index, active_index));
        ++index;
    }
    m_sessions->setTitles(titles);
    m_sessions->setItems(items);
}

void AppCore::captureCurrentSession() {
    const int active_index = m_sessions->activeIndex();
    if (active_index < 0 || active_index >= static_cast<int>(m_sessionStorage.size())) {
        return;
    }

    auto& session = m_sessionStorage[static_cast<std::size_t>(active_index)];
    session.title = m_notebook->sessionTitle();
    session.selectedIndex = m_notebook->selectedIndex();
    session.cells.clear();
    session.definitions.clear();
    for (const auto* cell : m_notebook->rawCells()) {
        session.cells.push_back(storage::StoredCell{
            cell->inputLatex(),
            cell->statusText(),
            cell->outputLatex(),
            cell->outputMeta(),
            cell->hasOutput(),
            cell->alternatives(),
            cell->steps(),
        });
    }
    if (m_kernel) {
        for (const auto& definition : m_kernel->snapshotDefinitions()) {
            session.definitions.push_back(storage::StoredDefinition{
                QString::fromStdString(definition.name),
                QString::fromStdString(definition.value_text),
            });
        }
    }

    if (session.cells.empty()) {
        session.cells.push_back(storage::StoredCell{});
        session.selectedIndex = 0;
    }

    syncSessionTitles();
}

void AppCore::loadSession(int index) {
    if (index < 0 || index >= static_cast<int>(m_sessionStorage.size())) {
        return;
    }

    m_notebook->clearAllCells();
    m_focusedCellIndex = -1;
    m_editingCellIndex = -1;
    const auto& session = m_sessionStorage[static_cast<std::size_t>(index)];
    m_notebook->setSessionTitle(session.title);
    if (m_kernel) {
        std::vector<cas::gui::CasGuiSession::StoredDefinition> definitions;
        definitions.reserve(session.definitions.size());
        for (const auto& definition : session.definitions) {
            definitions.push_back(cas::gui::CasGuiSession::StoredDefinition{
                definition.name.toStdString(),
                definition.valueText.toStdString(),
            });
        }
        auto restored = m_kernel->restoreDefinitions(definitions);
        if (restored.is_error()) {
            setStorageStatus(QString("Definition restore failed: %1")
                                 .arg(QString::fromStdString(restored.error().message)));
        }
    }

    for (const auto& stored : session.cells) {
        auto* cell = addTrackedCell();
        cell->setInputLatex(stored.input);
        if (stored.hasOutput) {
            cell->setOutputState(
                stored.status,
                true,
                stored.output,
                stored.meta,
                stored.alternatives,
                stored.steps,
                true);
        } else {
            cell->clearOutput();
            cell->setRestoredFromWorkspace(true);
        }
    }

    if (m_notebook->cellCount() == 0) {
        addTrackedCell();
    }
    const int target_index = std::clamp(session.selectedIndex, 0, m_notebook->cellCount() - 1);
    m_notebook->setSelectedIndex(target_index);
    syncCellInteractionState();
    emit selectedCellChanged();
    emit variablesChanged();
}

void AppCore::restoreWorkspace() {
    m_restoringWorkspace = true;
    m_hasPendingEditUndo = false;
    m_editUndoTimer.stop();
    m_undoStack.clear();
    m_redoStack.clear();
    emit undoAvailabilityChanged();
    setStorageNotification(QString());

    const auto loaded = storage::SessionStore::load();
    if (!loaded.error.isEmpty()) {
        m_sessionStorage.clear();
        initializeDefaultSession();
        m_restoringWorkspace = false;
        schedulePersist(
            QString("Recovered workspace after load failure: %1; autosave pending")
                .arg(loaded.error), 3);
        setStorageNotification(QString("Workspace load failed: %1. Recovered to default.").arg(loaded.error));
        return;
    }

    if (!loaded.found) {
        m_sessionStorage.clear();
        initializeDefaultSession();
        m_restoringWorkspace = false;
        if (loaded.warning.isEmpty()) {
            schedulePersist(QStringLiteral("New workspace; autosave pending"));
        } else {
            schedulePersist(
                QString("Recovered empty workspace: %1; autosave pending").arg(loaded.warning), 2);
            setStorageNotification(QString("Workspace was empty or corrupted: %1").arg(loaded.warning));
        }
        return;
    }

    restoreWorkspaceSnapshot(loaded.workspace);
    if (loaded.warning.isEmpty()) {
        setStorageStatus(QStringLiteral("Workspace restored"), 0);
    } else {
        setStorageStatus(QString("Workspace restored with recovery: %1").arg(loaded.warning), 2);
        setStorageNotification(QString("Workspace restored with some fixes: %1").arg(loaded.warning));
    }
}

void AppCore::restoreWorkspaceSnapshot(const storage::StoredWorkspace& workspace) {
    m_restoringWorkspace = true;
    m_sessionStorage = workspace.sessions;
    if (m_sessionStorage.empty()) {
        initializeDefaultSession();
        m_restoringWorkspace = false;
        return;
    }

    syncSessionTitles();
    const int active_index =
        std::clamp(workspace.activeSessionIndex, 0, static_cast<int>(m_sessionStorage.size()) - 1);
    m_sessions->setActiveIndex(active_index);
    loadSession(active_index);
    syncSessionTitles();
    if (m_plot) {
        m_plot->setVariable(workspace.plot.variable);
        m_plot->setXMin(workspace.plot.xMin);
        m_plot->setXMax(workspace.plot.xMax);
        m_plot->setIdle(QStringLiteral("Plot restored from workspace"));
    }
    m_sidebarVisible = workspace.ui.sidebarVisible;
    m_inspectorVisible = workspace.ui.inspectorVisible;
    m_sidebarWidth = workspace.ui.sidebarWidth;
    m_inspectorWidth = workspace.ui.inspectorWidth;
    m_notebookScrollY = workspace.ui.notebookScrollY;
    m_selectedInspectorTab = workspace.selectedInspectorTab.isEmpty()
                                 ? QStringLiteral("Result")
                                 : workspace.selectedInspectorTab;
    m_selectedRepresentationId = workspace.selectedRepresentationId;
    emit layoutStateChanged();
    emit selectedInspectorTabChanged();
    emit selectedRepresentationChanged();
    ensureSelectedRepresentationStillValid();
    m_restoringWorkspace = false;
}

void AppCore::schedulePersist(const QString& pendingStatus, int severity) {
    setStorageStatus(pendingStatus, severity);
    m_autosaveTimer.start();
}

void AppCore::persistWorkspace() {
    if (m_restoringWorkspace || m_sessionStorage.empty()) {
        return;
    }

    commitPendingEditUndo();
    captureCurrentSession();
    storage::StoredWorkspace workspace = snapshotWorkspace();

    QString error;
    if (!storage::SessionStore::save(workspace, &error)) {
        setStorageStatus(QString("Save failed: %1").arg(error), 3);
        return;
    }
    setStorageStatus("Saved", 0);
}

void AppCore::setStorageStatus(const QString& status, int severity) {
    if (m_storageStatus == status && m_storageSeverity == severity) {
        return;
    }
    m_storageStatus = status;
    m_storageSeverity = severity;
    emit storageStatusChanged();
}

void AppCore::setStorageNotification(const QString& message) {
    if (m_storageNotification == message) {
        return;
    }
    m_storageNotification = message;
    emit storageNotificationChanged();
}

QString AppCore::ensureUniqueSessionTitle(const QString& baseTitle, int skipIndex) const {
    QString candidate = baseTitle.trimmed();
    if (candidate.isEmpty()) {
        candidate = "Sessione";
    }

    auto is_used = [&](const QString& title) {
        for (std::size_t i = 0; i < m_sessionStorage.size(); ++i) {
            if (static_cast<int>(i) == skipIndex) {
                continue;
            }
            if (m_sessionStorage[i].title == title) {
                return true;
            }
        }
        return false;
    };

    if (!is_used(candidate)) {
        return candidate;
    }

    for (int suffix = 2; suffix < 10'000; ++suffix) {
        const QString next = QString("%1 (%2)").arg(candidate).arg(suffix);
        if (!is_used(next)) {
            return next;
        }
    }
    return QString("%1 (%2)").arg(candidate).arg(QDateTime::currentMSecsSinceEpoch());
}

storage::StoredWorkspace AppCore::snapshotWorkspace() const {
    storage::StoredWorkspace workspace;
    workspace.sessions = m_sessionStorage;
    workspace.activeSessionIndex = m_sessions ? m_sessions->activeIndex() : 0;
    if (m_plot) {
        workspace.plot.variable = m_plot->variable();
        workspace.plot.xMin = m_plot->xMin();
        workspace.plot.xMax = m_plot->xMax();
    }
    workspace.ui.sidebarVisible = m_sidebarVisible;
    workspace.ui.inspectorVisible = m_inspectorVisible;
    workspace.ui.sidebarWidth = m_sidebarWidth;
    workspace.ui.inspectorWidth = m_inspectorWidth;
    workspace.ui.notebookScrollY = m_notebookScrollY;
    workspace.selectedInspectorTab = m_selectedInspectorTab;
    workspace.selectedRepresentationId = m_selectedRepresentationId;
    return workspace;
}

bool AppCore::workspaceEquals(const storage::StoredWorkspace& lhs,
                              const storage::StoredWorkspace& rhs) const {
    if (lhs.plot.variable != rhs.plot.variable ||
        lhs.plot.xMin != rhs.plot.xMin ||
        lhs.plot.xMax != rhs.plot.xMax ||
        lhs.ui.sidebarVisible != rhs.ui.sidebarVisible ||
        lhs.ui.inspectorVisible != rhs.ui.inspectorVisible ||
        lhs.ui.sidebarWidth != rhs.ui.sidebarWidth ||
        lhs.ui.inspectorWidth != rhs.ui.inspectorWidth ||
        lhs.ui.notebookScrollY != rhs.ui.notebookScrollY ||
        lhs.selectedInspectorTab != rhs.selectedInspectorTab ||
        lhs.selectedRepresentationId != rhs.selectedRepresentationId) {
        return false;
    }

    if (lhs.activeSessionIndex != rhs.activeSessionIndex ||
        lhs.sessions.size() != rhs.sessions.size()) {
        return false;
    }

    for (std::size_t s = 0; s < lhs.sessions.size(); ++s) {
        const auto& a = lhs.sessions[s];
        const auto& b = rhs.sessions[s];
        if (a.title != b.title ||
            a.selectedIndex != b.selectedIndex ||
            a.definitions.size() != b.definitions.size() ||
            a.cells.size() != b.cells.size()) {
            return false;
        }
        for (std::size_t d = 0; d < a.definitions.size(); ++d) {
            const auto& ad = a.definitions[d];
            const auto& bd = b.definitions[d];
            if (ad.name != bd.name || ad.valueText != bd.valueText) {
                return false;
            }
        }
        for (std::size_t c = 0; c < a.cells.size(); ++c) {
            const auto& ac = a.cells[c];
            const auto& bc = b.cells[c];
            if (ac.input != bc.input ||
                ac.status != bc.status ||
                ac.output != bc.output ||
                ac.meta != bc.meta ||
                ac.alternatives != bc.alternatives ||
                ac.steps != bc.steps ||
                ac.hasOutput != bc.hasOutput) {
                return false;
            }
        }
    }
    return true;
}

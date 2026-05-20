#include "AppCore.h"

#include <QDateTime>
#include <QLoggingCategory>

namespace {
constexpr std::size_t kUndoLimit = 200;
}

Q_DECLARE_LOGGING_CATEGORY(casGuiAuditLog)

void AppCore::finalizeUndoableChange(const QString& label, const storage::StoredWorkspace& before) {
    if (isUndoTransactionActive()) {
        return;
    }
    captureCurrentSession();
    const auto after = snapshotWorkspace();
    if (workspaceEquals(before, after)) {
        return;
    }

    pushUndoStep(label, before, after);
    schedulePersist();
}

void AppCore::pushUndoStep(const QString& label,
                           const storage::StoredWorkspace& before,
                           const storage::StoredWorkspace& after) {
    qCDebug(casGuiAuditLog).noquote()
        << "[Audit] [" << QDateTime::currentDateTime().toString(Qt::ISODateWithMs)
        << "] pushUndoStep:" << label;
    m_undoStack.push_back(UndoStep{label, before, after});
    if (m_undoStack.size() > kUndoLimit) {
        m_undoStack.erase(m_undoStack.begin());
    }
    m_redoStack.clear();
    emit undoAvailabilityChanged();
}

storage::StoredWorkspace AppCore::snapshotBeforeUndoableAction() {
    if (isUndoTransactionActive()) {
        return m_undoTransactionBefore;
    }
    commitPendingEditUndo();
    return snapshotWorkspace();
}

void AppCore::markPendingEditUndo(const QString& label) {
    if (m_restoringWorkspace || isUndoTransactionActive()) {
        return;
    }
    if (!m_hasPendingEditUndo) {
        m_pendingEditBefore = snapshotWorkspace();
        m_hasPendingEditUndo = true;
        m_pendingEditUndoLabel = label;
    } else if (m_pendingEditUndoLabel != label) {
        commitPendingEditUndo();
        m_pendingEditBefore = snapshotWorkspace();
        m_hasPendingEditUndo = true;
        m_pendingEditUndoLabel = label;
    }
    m_editUndoTimer.start();
}

void AppCore::commitPendingEditUndo() {
    if (isUndoTransactionActive()) {
        return;
    }
    if (!m_hasPendingEditUndo) {
        return;
    }
    m_editUndoTimer.stop();
    m_hasPendingEditUndo = false;
    captureCurrentSession();
    const auto after = snapshotWorkspace();
    if (!workspaceEquals(m_pendingEditBefore, after)) {
        pushUndoStep(m_pendingEditUndoLabel, m_pendingEditBefore, after);
    }
    m_pendingEditUndoLabel = "Edit cell";
}

void AppCore::beginUndoTransaction(const QString& label) {
    if (m_undoTransactionDepth == 0) {
        commitPendingEditUndo();
        m_undoTransactionBefore = snapshotWorkspace();
        m_undoTransactionLabel = label;
        m_undoTransactionCommit = true;
    }
    ++m_undoTransactionDepth;
}

void AppCore::endUndoTransaction(bool commit) {
    if (m_undoTransactionDepth <= 0) {
        return;
    }
    if (!commit) {
        m_undoTransactionCommit = false;
    }
    --m_undoTransactionDepth;
    if (m_undoTransactionDepth > 0) {
        return;
    }

    if (m_undoTransactionCommit) {
        captureCurrentSession();
        const auto after = snapshotWorkspace();
        if (!workspaceEquals(m_undoTransactionBefore, after)) {
            pushUndoStep(m_undoTransactionLabel, m_undoTransactionBefore, after);
            schedulePersist();
        }
    }

    m_undoTransactionCommit = true;
    m_undoTransactionLabel.clear();
}

#include "AppCore.h"

#include <QtGlobal>

namespace {

constexpr qreal kSidebarWidthDefault = 230.0;
constexpr qreal kSidebarWidthMin = 180.0;
constexpr qreal kSidebarWidthMax = 420.0;
constexpr qreal kInspectorWidthDefault = 320.0;
constexpr qreal kInspectorWidthMin = 260.0;
constexpr qreal kInspectorWidthMax = 520.0;

} // namespace

void AppCore::setSidebarVisible(bool visible) {
    if (m_sidebarVisible == visible) {
        return;
    }
    m_sidebarVisible = visible;
    emit layoutStateChanged();
    if (!m_restoringWorkspace && !isUndoTransactionActive()) {
        schedulePersist();
    }
}

void AppCore::setInspectorVisible(bool visible) {
    if (m_inspectorVisible == visible) {
        return;
    }
    m_inspectorVisible = visible;
    emit layoutStateChanged();
    if (!m_restoringWorkspace && !isUndoTransactionActive()) {
        schedulePersist();
    }
}

void AppCore::setSidebarWidth(qreal width) {
    const qreal clamped = qBound(kSidebarWidthMin, width, kSidebarWidthMax);
    if (qFuzzyCompare(m_sidebarWidth + 1.0, clamped + 1.0)) {
        return;
    }
    m_sidebarWidth = clamped;
    emit layoutStateChanged();
    if (!m_restoringWorkspace && !isUndoTransactionActive()) {
        schedulePersist();
    }
}

void AppCore::setInspectorWidth(qreal width) {
    const qreal clamped = qBound(kInspectorWidthMin, width, kInspectorWidthMax);
    if (qFuzzyCompare(m_inspectorWidth + 1.0, clamped + 1.0)) {
        return;
    }
    m_inspectorWidth = clamped;
    emit layoutStateChanged();
    if (!m_restoringWorkspace && !isUndoTransactionActive()) {
        schedulePersist();
    }
}

void AppCore::setNotebookScrollY(qreal scrollY) {
    const qreal clamped = scrollY < 0.0 ? 0.0 : scrollY;
    if (qFuzzyCompare(m_notebookScrollY + 1.0, clamped + 1.0)) {
        return;
    }
    m_notebookScrollY = clamped;
    emit layoutStateChanged();
    if (!m_restoringWorkspace && !isUndoTransactionActive()) {
        schedulePersist();
    }
}

void AppCore::toggleSidebar() {
    setSidebarVisible(!m_sidebarVisible);
}

void AppCore::toggleInspector() {
    setInspectorVisible(!m_inspectorVisible);
}

void AppCore::resetLayoutState() {
    const bool changed = m_sidebarVisible != true ||
                         m_inspectorVisible != true ||
                         !qFuzzyCompare(m_sidebarWidth + 1.0, kSidebarWidthDefault + 1.0) ||
                         !qFuzzyCompare(m_inspectorWidth + 1.0, kInspectorWidthDefault + 1.0) ||
                         !qFuzzyCompare(m_notebookScrollY + 1.0, 1.0);
    m_sidebarVisible = true;
    m_inspectorVisible = true;
    m_sidebarWidth = kSidebarWidthDefault;
    m_inspectorWidth = kInspectorWidthDefault;
    m_notebookScrollY = 0.0;
    if (!changed) {
        return;
    }
    emit layoutStateChanged();
    if (!m_restoringWorkspace && !isUndoTransactionActive()) {
        schedulePersist();
    }
}

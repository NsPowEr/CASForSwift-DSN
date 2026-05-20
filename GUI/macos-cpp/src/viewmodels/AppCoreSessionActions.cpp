#include "AppCore.h"
#include <algorithm>

void AppCore::newSession() {
    beginUndoTransaction("New session");
    captureCurrentSession();
    storage::StoredSession session;
    session.title = ensureUniqueSessionTitle(QString("Sessione %1").arg(m_sessionStorage.size() + 1));
    session.cells.push_back(storage::StoredCell{});
    m_sessionStorage.push_back(std::move(session));
    syncSessionTitles();
    const int new_index = static_cast<int>(m_sessionStorage.size()) - 1;
    m_sessions->setActiveIndex(new_index);
    loadSession(new_index);
    syncSessionTitles();
    endUndoTransaction(true);
}

void AppCore::openSession(int index) {
    if (index == m_sessions->activeIndex()) {
        return;
    }
    captureCurrentSession();
    m_sessions->setActiveIndex(index);
    loadSession(index);
    syncSessionTitles();
}

void AppCore::renameActiveSession(const QString& title) {
    const int active_index = m_sessions->activeIndex();
    if (active_index < 0 || active_index >= static_cast<int>(m_sessionStorage.size())) {
        return;
    }

    QString normalized = title.trimmed();
    if (normalized.isEmpty()) {
        normalized = QString("Sessione %1").arg(active_index + 1);
    }
    normalized = ensureUniqueSessionTitle(normalized, active_index);

    beginUndoTransaction("Rename session");
    m_notebook->setSessionTitle(normalized);
    endUndoTransaction(true);
}

void AppCore::duplicateSession(int index) {
    if (m_sessionStorage.empty()) {
        return;
    }

    if (index < 0) {
        index = m_sessions->activeIndex();
    }
    if (index < 0 || index >= static_cast<int>(m_sessionStorage.size())) {
        return;
    }

    beginUndoTransaction("Duplicate session");
    captureCurrentSession();
    storage::StoredSession duplicated = m_sessionStorage[static_cast<std::size_t>(index)];
    duplicated.title = ensureUniqueSessionTitle(QString("%1 copy").arg(duplicated.title));
    if (duplicated.cells.empty()) {
        duplicated.cells.push_back(storage::StoredCell{});
        duplicated.selectedIndex = 0;
    } else {
        duplicated.selectedIndex = std::clamp(duplicated.selectedIndex, 0, static_cast<int>(duplicated.cells.size()) - 1);
    }

    const int insert_index = index + 1;
    m_sessionStorage.insert(m_sessionStorage.begin() + insert_index, std::move(duplicated));
    syncSessionTitles();
    m_sessions->setActiveIndex(insert_index);
    loadSession(insert_index);
    syncSessionTitles();
    endUndoTransaction(true);
}

void AppCore::deleteSession(int index) {
    if (m_sessionStorage.size() <= 1) {
        setStorageStatus("Cannot delete last session");
        return;
    }

    if (index < 0) {
        index = m_sessions->activeIndex();
    }
    if (index < 0 || index >= static_cast<int>(m_sessionStorage.size())) {
        return;
    }

    beginUndoTransaction("Delete session");
    captureCurrentSession();
    m_sessionStorage.erase(m_sessionStorage.begin() + index);
    const int target_index = std::clamp(index, 0, static_cast<int>(m_sessionStorage.size()) - 1);
    syncSessionTitles();
    m_sessions->setActiveIndex(target_index);
    loadSession(target_index);
    syncSessionTitles();
    endUndoTransaction(true);
}

void AppCore::moveSessionUp(int index) {
    if (m_sessionStorage.size() <= 1) {
        return;
    }

    if (index < 0) {
        index = m_sessions->activeIndex();
    }
    if (index <= 0 || index >= static_cast<int>(m_sessionStorage.size())) {
        return;
    }

    beginUndoTransaction("Move session up");
    captureCurrentSession();
    std::swap(m_sessionStorage[static_cast<std::size_t>(index)],
              m_sessionStorage[static_cast<std::size_t>(index - 1)]);
    syncSessionTitles();
    m_sessions->setActiveIndex(index - 1);
    loadSession(index - 1);
    syncSessionTitles();
    endUndoTransaction(true);
}

void AppCore::moveSessionDown(int index) {
    if (m_sessionStorage.size() <= 1) {
        return;
    }

    if (index < 0) {
        index = m_sessions->activeIndex();
    }
    if (index < 0 || index >= static_cast<int>(m_sessionStorage.size()) - 1) {
        return;
    }

    beginUndoTransaction("Move session down");
    captureCurrentSession();
    std::swap(m_sessionStorage[static_cast<std::size_t>(index)],
              m_sessionStorage[static_cast<std::size_t>(index + 1)]);
    syncSessionTitles();
    m_sessions->setActiveIndex(index + 1);
    loadSession(index + 1);
    syncSessionTitles();
    endUndoTransaction(true);
}

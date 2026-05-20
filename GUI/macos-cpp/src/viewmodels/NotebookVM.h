#pragma once
#include <algorithm>
#include <QObject>
#include <QQmlListProperty>
#include <vector>
#include "CellVM.h"

namespace cas::gui {

class NotebookVM : public QObject {
    Q_OBJECT
    Q_PROPERTY(QQmlListProperty<cas::gui::CellVM> cells READ cells NOTIFY cellsChanged)
    Q_PROPERTY(int selectedIndex READ selectedIndex WRITE setSelectedIndex NOTIFY selectedIndexChanged)
    Q_PROPERTY(QString sessionTitle READ sessionTitle NOTIFY sessionTitleChanged)

public:
    explicit NotebookVM(QObject* parent = nullptr) : QObject(parent) {}

    QQmlListProperty<CellVM> cells() {
        return QQmlListProperty<CellVM>(this, &m_cellsList,
            &NotebookVM::appendCell, &NotebookVM::cellsCount,
            &NotebookVM::cellAt, &NotebookVM::clearCells);
    }

    int selectedIndex() const { return m_selectedIndex; }
    void setSelectedIndex(int index) {
        if (index == m_selectedIndex) return;
        if (m_selectedIndex >= 0 && m_selectedIndex < m_cellsList.size())
            m_cellsList[m_selectedIndex]->setActive(false);
        m_selectedIndex = index;
        if (m_selectedIndex >= 0 && m_selectedIndex < m_cellsList.size())
            m_cellsList[m_selectedIndex]->setActive(true);
        emit selectedIndexChanged();
    }

    QString sessionTitle() const { return m_title; }
    void setSessionTitle(const QString& title) {
        if (m_title == title) return;
        m_title = title;
        emit sessionTitleChanged();
    }

    CellVM* addEmptyCell() {
        auto* cell = new CellVM(m_cellsList.size() + 1, this);
        m_cellsList.append(cell);
        reindexCells();
        emit cellsChanged();
        setSelectedIndex(m_cellsList.size() - 1);
        return cell;
    }

    CellVM* insertCell(int index,
                       const QString& input,
                       const QString& status,
                       bool hasOutput,
                       const QString& output,
                       const QString& meta,
                       const QVariantList& alternatives = {},
                       const QVariantList& steps = {}) {
        const int clamped_index = std::clamp(index, 0, static_cast<int>(m_cellsList.size()));
        auto* cell = new CellVM(0, this);
        cell->setInputLatex(input);
        if (hasOutput) {
            cell->setOutputState(status, true, output, meta, alternatives, steps);
        } else {
            cell->clearOutput();
        }
        m_cellsList.insert(clamped_index, cell);
        reindexCells();
        emit cellsChanged();
        setSelectedIndex(clamped_index);
        return cell;
    }

    CellVM* insertEmptyCellAt(int index) {
        return insertCell(index, QString(), QStringLiteral("ready"), false, QString(), QString());
    }

    CellVM* currentCell() const {
        if (m_selectedIndex >= 0 && m_selectedIndex < m_cellsList.size())
            return m_cellsList[m_selectedIndex];
        return nullptr;
    }
    const QList<CellVM*>& rawCells() const { return m_cellsList; }
    int cellCount() const { return m_cellsList.size(); }

    CellVM* duplicateCell(int index) {
        if (index < 0 || index >= m_cellsList.size()) {
            return nullptr;
        }
        auto* source = m_cellsList[index];
        auto* cell = new CellVM(0, this);
        cell->setInputLatex(source->inputLatex());
        cell->setOutputState(
            source->statusText(),
            source->hasOutput(),
            source->outputLatex(),
            source->outputMeta(),
            source->alternatives(),
            source->steps());
        m_cellsList.insert(index + 1, cell);
        reindexCells();
        emit cellsChanged();
        setSelectedIndex(index + 1);
        return cell;
    }

    bool deleteCell(int index) {
        if (index < 0 || index >= m_cellsList.size() || m_cellsList.size() <= 1) {
            return false;
        }
        auto* cell = m_cellsList.takeAt(index);
        delete cell;
        reindexCells();
        emit cellsChanged();
        const int max_index = static_cast<int>(m_cellsList.size()) - 1;
        const int target = std::clamp(index, 0, max_index);
        setSelectedIndex(target);
        return true;
    }

    bool moveCellUp(int index) {
        if (index <= 0 || index >= m_cellsList.size()) {
            return false;
        }
        m_cellsList.swapItemsAt(index, index - 1);
        reindexCells();
        emit cellsChanged();
        setSelectedIndex(index - 1);
        return true;
    }

    bool moveCellDown(int index) {
        if (index < 0 || index >= m_cellsList.size() - 1) {
            return false;
        }
        m_cellsList.swapItemsAt(index, index + 1);
        reindexCells();
        emit cellsChanged();
        setSelectedIndex(index + 1);
        return true;
    }

    bool selectPreviousCell() {
        if (m_cellsList.isEmpty()) {
            return false;
        }
        const int max_index = static_cast<int>(m_cellsList.size()) - 1;
        const int target = std::clamp(m_selectedIndex - 1, 0, max_index);
        if (target == m_selectedIndex) {
            return false;
        }
        setSelectedIndex(target);
        return true;
    }

    bool selectNextCell() {
        if (m_cellsList.isEmpty()) {
            return false;
        }
        const int max_index = static_cast<int>(m_cellsList.size()) - 1;
        const int target = std::clamp(m_selectedIndex + 1, 0, max_index);
        if (target == m_selectedIndex) {
            return false;
        }
        setSelectedIndex(target);
        return true;
    }

    void clearAllCells() {
        qDeleteAll(m_cellsList);
        m_cellsList.clear();
        m_selectedIndex = -1;
        emit cellsChanged();
        emit selectedIndexChanged();
    }

signals:
    void cellsChanged();
    void selectedIndexChanged();
    void sessionTitleChanged();

private:
    void reindexCells() {
        for (int i = 0; i < m_cellsList.size(); ++i) {
            m_cellsList[i]->setIndex(i + 1);
        }
    }

    static void appendCell(QQmlListProperty<CellVM>* p, CellVM* c) {
        static_cast<QList<CellVM*>*>(p->data)->append(c);
    }
    static qsizetype cellsCount(QQmlListProperty<CellVM>* p) {
        return static_cast<QList<CellVM*>*>(p->data)->size();
    }
    static CellVM* cellAt(QQmlListProperty<CellVM>* p, qsizetype i) {
        return static_cast<QList<CellVM*>*>(p->data)->at(i);
    }
    static void clearCells(QQmlListProperty<CellVM>* p) {
        static_cast<QList<CellVM*>*>(p->data)->clear();
    }

    QList<CellVM*> m_cellsList;
    int m_selectedIndex = -1;
    QString m_title = "Nuova Sessione";
};

} // namespace cas::gui

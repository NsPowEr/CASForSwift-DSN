#include "AppCore.h"

namespace {

QString join_cell_inputs(const QString& left, const QString& right) {
    if (left.isEmpty()) {
        return right;
    }
    if (right.isEmpty()) {
        return left;
    }
    return left + QLatin1Char('\n') + right;
}

}

void AppCore::addEmptyCell() {
    beginUndoTransaction("Add cell");
    addTrackedCell();
    beginEditingCell();
    endUndoTransaction(true);
    emit selectedCellChanged();
}

void AppCore::insertCellAbove() {
    if (!m_notebook) {
        return;
    }
    beginUndoTransaction("Insert cell above");
    const int index = std::max(0, m_notebook->selectedIndex());
    auto* cell = m_notebook->insertEmptyCellAt(index);
    attachCell(cell);
    beginEditingCell();
    endUndoTransaction(true);
    emit selectedCellChanged();
}

void AppCore::insertCellBelow() {
    if (!m_notebook) {
        return;
    }
    beginUndoTransaction("Insert cell below");
    const int index = std::max(0, m_notebook->selectedIndex());
    auto* cell = m_notebook->insertEmptyCellAt(index + 1);
    attachCell(cell);
    beginEditingCell();
    endUndoTransaction(true);
    emit selectedCellChanged();
}

void AppCore::duplicateCurrentCell() {
    if (!m_notebook) {
        return;
    }
    beginUndoTransaction("Duplicate cell");
    const int index = m_notebook->selectedIndex();
    auto* duplicated = m_notebook->duplicateCell(index);
    if (!duplicated) {
        endUndoTransaction(false);
        return;
    }
    attachCell(duplicated);
    focusCell();
    endUndoTransaction(true);
    emit selectedCellChanged();
}

void AppCore::deleteCurrentCell() {
    if (!m_notebook) {
        return;
    }
    beginUndoTransaction("Delete cell");
    const int index = m_notebook->selectedIndex();
    if (!m_notebook->deleteCell(index)) {
        endUndoTransaction(false);
        setStorageStatus("Cannot delete last cell");
        return;
    }
    focusCell();
    endUndoTransaction(true);
    emit selectedCellChanged();
}

void AppCore::moveCurrentCellUp() {
    if (!m_notebook) {
        return;
    }
    beginUndoTransaction("Move cell up");
    if (!m_notebook->moveCellUp(m_notebook->selectedIndex())) {
        endUndoTransaction(false);
        return;
    }
    focusCell();
    endUndoTransaction(true);
    emit selectedCellChanged();
}

void AppCore::moveCurrentCellDown() {
    if (!m_notebook) {
        return;
    }
    beginUndoTransaction("Move cell down");
    if (!m_notebook->moveCellDown(m_notebook->selectedIndex())) {
        endUndoTransaction(false);
        return;
    }
    focusCell();
    endUndoTransaction(true);
    emit selectedCellChanged();
}

void AppCore::selectPreviousCell() {
    if (!m_notebook) {
        return;
    }
    if (m_notebook->selectPreviousCell()) {
        focusCell();
        emit selectedCellChanged();
    }
}

void AppCore::selectNextCell() {
    if (!m_notebook) {
        return;
    }
    if (m_notebook->selectNextCell()) {
        focusCell();
        emit selectedCellChanged();
    }
}

void AppCore::executeCurrentCell() {
    auto* current = m_notebook ? m_notebook->currentCell() : nullptr;
    if (!current && m_notebook && m_notebook->cellCount() == 0) {
        current = addTrackedCell();
    }

    if (!current || !m_kernel) {
        return;
    }

    const QString input = current->inputLatex().trimmed();
    if (input.isEmpty()) {
        return;
    }

    if (input.startsWith('?')) {
        lookupHelp(input);
        return;
    }

    beginUndoTransaction("Execute cell");
    m_kernelBusy = true;
    emit kernelBusyChanged();

    if (m_plot) {
        m_plot->setIdle(QStringLiteral("Plot idle"));
    }

    m_kernel->evaluateAsync(input, [this, current](ComputeResult result) {
        QMetaObject::invokeMethod(this, [this, current, result]() {
            m_kernelBusy = false;
            emit kernelBusyChanged();

            if (result.interrupted) {
                current->setOutputState("interrupted", false, QString(), "Calculation stopped by user");
                setStorageStatus("Interrupted", 2);
            } else {
                current->setResult(result);
            }

            ensureSelectedRepresentationStillValid();
            endUndoTransaction(true);
            emit selectedCellChanged();
            emit variablesChanged();
        });
    });
}

void AppCore::runAndAdvance() {
    auto* current = m_notebook->currentCell();
    if (!current && m_notebook->cellCount() == 0) {
        current = addTrackedCell();
    }

    if (current && m_kernel) {
        QString input = current->inputLatex().trimmed();
        if (input.isEmpty()) {
            return;
        }

        if (input.startsWith('?')) {
            lookupHelp(input);
            return;
        }

        beginUndoTransaction("Run and advance");

        if (m_plot) {
            m_plot->setIdle(QStringLiteral("Plot idle"));
        }

        m_kernel->evaluateAsync(input, [this, current](ComputeResult result) {
            QMetaObject::invokeMethod(this, [this, current, result]() {
                current->setResult(result);
                ensureSelectedRepresentationStillValid();
                if (m_notebook->currentCell() == current) {
                    const int current_index = m_notebook->selectedIndex();
                    if (current_index >= 0 && current_index < (m_notebook->cellCount() - 1)) {
                        m_notebook->setSelectedIndex(current_index + 1);
                    } else {
                        addTrackedCell();
                    }
                    beginEditingCell();
                }
                endUndoTransaction(true);
                emit selectedCellChanged();
                emit variablesChanged();
            });
        });
    }
}

void AppCore::splitCurrentCell(int cursorPosition) {
    if (!m_notebook) {
        return;
    }

    auto* current = currentCell();
    if (!current) {
        return;
    }

    const QString input = current->inputLatex();
    const int input_size = static_cast<int>(input.size());
    const int split_index =
        std::clamp(cursorPosition >= 0 ? cursorPosition : input_size, 0, input_size);
    const QString head = input.left(split_index);
    const QString tail = input.mid(split_index);

    beginUndoTransaction("Split cell");
    current->setInputLatex(head);
    auto* inserted = m_notebook->insertEmptyCellAt(m_notebook->selectedIndex() + 1);
    attachCell(inserted);
    inserted->setInputLatex(tail);
    beginEditingCell();
    endUndoTransaction(true);
    emit selectedCellChanged();
}

void AppCore::mergeWithPreviousCell() {
    if (!m_notebook) {
        return;
    }
    const int index = m_notebook->selectedIndex();
    if (index <= 0) {
        return;
    }

    auto* previous = m_notebook->rawCells().at(index - 1);
    auto* current = m_notebook->rawCells().at(index);
    beginUndoTransaction("Merge with previous");
    previous->setInputLatex(join_cell_inputs(previous->inputLatex(), current->inputLatex()));
    if (current->hasOutput() && !previous->hasOutput()) {
        previous->setOutputState(current->statusText(),
                                 current->hasOutput(),
                                 current->outputLatex(),
                                 current->outputMeta(),
                                 current->alternatives(),
                                 current->steps());
    }
    m_notebook->deleteCell(index);
    beginEditingCell(index - 1);
    endUndoTransaction(true);
    emit selectedCellChanged();
}

void AppCore::mergeWithNextCell() {
    if (!m_notebook) {
        return;
    }
    const int index = m_notebook->selectedIndex();
    if (index < 0 || index >= (m_notebook->cellCount() - 1)) {
        return;
    }

    auto* current = m_notebook->rawCells().at(index);
    auto* next = m_notebook->rawCells().at(index + 1);
    beginUndoTransaction("Merge with next");
    current->setInputLatex(join_cell_inputs(current->inputLatex(), next->inputLatex()));
    if (next->hasOutput() && !current->hasOutput()) {
        current->setOutputState(next->statusText(),
                                next->hasOutput(),
                                next->outputLatex(),
                                next->outputMeta(),
                                next->alternatives(),
                                next->steps());
    }
    m_notebook->deleteCell(index + 1);
    m_notebook->setSelectedIndex(index);
    beginEditingCell(index);
    endUndoTransaction(true);
    emit selectedCellChanged();
}

CellVM* AppCore::addTrackedCell() {
    auto* cell = m_notebook->addEmptyCell();
    attachCell(cell);
    return cell;
}

void AppCore::attachCell(CellVM* cell) {
    if (!cell) {
        return;
    }

    connect(cell, &CellVM::inputAboutToChange, this, [this]() {
        if (!m_restoringWorkspace && !isUndoTransactionActive()) {
            markPendingEditUndo("Edit cell");
        }
    });
    connect(cell, &CellVM::inputChanged, this, [this, cell]() {
        if (!m_restoringWorkspace && !isUndoTransactionActive()) {
            captureCurrentSession();
            schedulePersist();
        }
        if (currentCell() == cell) {
            emit selectedCellChanged();
        }
    });
    connect(cell, &CellVM::statusChanged, this, [this, cell]() {
        if (!m_restoringWorkspace && !isUndoTransactionActive()) {
            captureCurrentSession();
            schedulePersist();
        }
        if (currentCell() == cell) {
            ensureSelectedRepresentationStillValid();
            emit selectedCellChanged();
        }
    });
    connect(cell, &CellVM::outputChanged, this, [this, cell]() {
        if (!m_restoringWorkspace && !isUndoTransactionActive()) {
            captureCurrentSession();
            schedulePersist();
        }
        if (currentCell() == cell) {
            ensureSelectedRepresentationStillValid();
            emit selectedCellChanged();
        }
    });
}

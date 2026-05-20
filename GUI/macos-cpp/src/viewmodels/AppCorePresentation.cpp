#include "AppCore.h"

#include <QGuiApplication>
#include <QClipboard>

QString AppCore::selectedCellInput() const {
    auto* cell = currentCell();
    return cell ? cell->inputLatex() : QString();
}

QString AppCore::selectedCellStatus() const {
    auto* cell = currentCell();
    return cell ? cell->statusText() : QStringLiteral("No cell selected");
}

QString AppCore::selectedCellOutput() const {
    auto* cell = currentCell();
    return cell ? cell->outputLatex() : QString();
}

QString AppCore::selectedCellMeta() const {
    auto* cell = currentCell();
    return cell ? cell->outputMeta() : QString();
}

QString AppCore::selectedCellState() const {
    auto* cell = currentCell();
    if (!cell) {
        return QStringLiteral("no-selection");
    }
    if (cell->hasError()) {
        return QStringLiteral("error");
    }
    if (cell->restoredFromWorkspace()) {
        return QStringLiteral("restored");
    }
    if (cell->hasOutput()) {
        return QStringLiteral("result");
    }
    return QStringLiteral("ready");
}

bool AppCore::selectedCellHasError() const {
    auto* cell = currentCell();
    return cell && cell->hasError();
}

bool AppCore::selectedCellRestored() const {
    auto* cell = currentCell();
    return cell && cell->restoredFromWorkspace();
}

bool AppCore::selectedCellHasOutput() const {
    auto* cell = currentCell();
    return cell && !cell->outputLatex().isEmpty();
}

QVariantList AppCore::selectedCellAlternatives() const {
    auto* cell = currentCell();
    return cell ? cell->alternatives() : QVariantList{};
}

QVariantList AppCore::selectedCellSteps() const {
    auto* cell = currentCell();
    return cell ? cell->steps() : QVariantList{};
}

void AppCore::setVariableSearchQuery(const QString& value) {
    if (m_variableSearchQuery == value) {
        return;
    }
    m_variableSearchQuery = value;
    emit variableSearchQueryChanged();
    emit variablesChanged();
}

QVariantList AppCore::filteredVariables() const {
    const QString query = m_variableSearchQuery.trimmed().toLower();
    const QVariantList all_variables = variables();
    if (query.isEmpty()) {
        return all_variables;
    }

    QVariantList filtered;
    for (const auto& item : all_variables) {
        const QVariantMap map = item.toMap();
        const QString haystack =
            (map.value(QStringLiteral("name")).toString() + QLatin1Char(' ') +
             map.value(QStringLiteral("value")).toString())
                .toLower();
        if (haystack.contains(query)) {
            filtered.push_back(map);
        }
    }
    return filtered;
}

void AppCore::setSelectedInspectorTab(const QString& value) {
    if (value.isEmpty() || m_selectedInspectorTab == value) {
        return;
    }
    m_selectedInspectorTab = value;
    emit selectedInspectorTabChanged();
    if (!m_restoringWorkspace && !isUndoTransactionActive()) {
        schedulePersist();
    }
}

void AppCore::setSelectedRepresentationId(const QString& value) {
    if (m_selectedRepresentationId == value) {
        return;
    }
    m_selectedRepresentationId = value;
    emit selectedRepresentationChanged();
    if (!m_restoringWorkspace && !isUndoTransactionActive()) {
        schedulePersist();
    }
}

void AppCore::plotSelectedCell() {
    auto* cell = currentCell();
    if (!cell || !m_kernel || !m_plot) {
        return;
    }

    const QString input = cell->inputLatex().trimmed();
    if (input.isEmpty()) {
        m_plot->setError(QStringLiteral("Plot skipped: empty input"));
        return;
    }
    const QString variable = m_plot->variable().trimmed();
    if (variable.isEmpty()) {
        m_plot->setError(QStringLiteral("Plot skipped: empty variable"));
        return;
    }
    if (!(m_plot->xMin() < m_plot->xMax())) {
        m_plot->setError(QStringLiteral("Plot skipped: xMin must be smaller than xMax"));
        return;
    }

    m_plot->beginLoading();
    m_kernel->sample2DAsync(input, variable, m_plot->xMin(), m_plot->xMax(), [this, input](auto sampled) {
        QMetaObject::invokeMethod(this, [this, sampled, input]() {
            if (sampled.is_error()) {
                m_plot->setError(QString::fromStdString(sampled.error().message));
                return;
            }

            QVariantList points;
            for (const auto& point : sampled.value()) {
                points.append(QPointF(point.x, point.y));
            }
            m_plot->addSeries(input, points);
        });
    });
}

void AppCore::zoomPlotIn() {
    if (m_plot) {
        m_plot->zoomIn();
    }
}

void AppCore::zoomPlotOut() {
    if (m_plot) {
        m_plot->zoomOut();
    }
}

void AppCore::panPlotLeft() {
    if (m_plot) {
        m_plot->panLeft();
    }
}

void AppCore::panPlotRight() {
    if (m_plot) {
        m_plot->panRight();
    }
}

void AppCore::resetPlotRange() {
    if (m_plot) {
        m_plot->resetRange();
    }
}

void AppCore::applyPlotPreset(const QString& presetId) {
    if (m_plot) {
        m_plot->applyPreset(presetId);
    }
}

QVariantMap AppCore::variableByName(const QString& name) const {
    if (name.trimmed().isEmpty()) {
        return {};
    }
    for (const auto& item : variables()) {
        const QVariantMap map = item.toMap();
        if (map.value(QStringLiteral("name")).toString() == name) {
            return map;
        }
    }
    return {};
}

bool AppCore::insertTextIntoCurrentCell(const QString& text) {
    auto* cell = currentCell();
    if (!cell) {
        setStorageStatus(QStringLiteral("No cell selected"), 1);
        return false;
    }

    const QString trimmed = text.trimmed();
    if (trimmed.isEmpty()) {
        setStorageStatus(QStringLiteral("Nothing to insert"), 1);
        return false;
    }

    QString next = cell->inputLatex();
    if (!next.isEmpty()) {
        const QChar last = next.back();
        if (!last.isSpace() &&
            (last.isLetterOrNumber() || last == QLatin1Char(')') || last == QLatin1Char('}'))) {
            next += QLatin1Char(' ');
        }
    }
    next += trimmed;
    cell->setInputLatex(next);
    setStorageStatus(QStringLiteral("Inserted into current cell"), 0);
    return true;
}

void AppCore::copyVariableName(const QString& name) {
    const QVariantMap map = variableByName(name);
    if (map.isEmpty()) {
        setStorageStatus(QStringLiteral("Variable not found"), 1);
        return;
    }

    if (auto* clipboard = QGuiApplication::clipboard()) {
        clipboard->setText(map.value(QStringLiteral("name")).toString());
        setStorageStatus(QStringLiteral("Variable name copied"), 0);
    }
}

void AppCore::copyVariableValue(const QString& name) {
    const QVariantMap map = variableByName(name);
    if (map.isEmpty()) {
        setStorageStatus(QStringLiteral("Variable not found"), 1);
        return;
    }

    if (auto* clipboard = QGuiApplication::clipboard()) {
        clipboard->setText(map.value(QStringLiteral("value")).toString());
        setStorageStatus(QStringLiteral("Variable value copied"), 0);
    }
}

void AppCore::insertVariableName(const QString& name) {
    const QVariantMap map = variableByName(name);
    if (map.isEmpty()) {
        setStorageStatus(QStringLiteral("Variable not found"), 1);
        return;
    }
    insertTextIntoCurrentCell(map.value(QStringLiteral("name")).toString());
}

void AppCore::insertVariableValue(const QString& name) {
    const QVariantMap map = variableByName(name);
    if (map.isEmpty()) {
        setStorageStatus(QStringLiteral("Variable not found"), 1);
        return;
    }
    insertTextIntoCurrentCell(map.value(QStringLiteral("value")).toString());
}

void AppCore::invokeCommand(const QString& commandId) {
    if (!isCommandEnabled(commandId)) {
        setStorageStatus(commandDisabledReason(commandId).isEmpty()
                             ? QStringLiteral("Command unavailable")
                             : commandDisabledReason(commandId),
                         1);
        return;
    }

    if (commandId == QStringLiteral("run")) {
        executeCurrentCell();
    } else if (commandId == QStringLiteral("run_and_advance")) {
        runAndAdvance();
    } else if (commandId == QStringLiteral("plot")) {
        plotSelectedCell();
    } else if (commandId == QStringLiteral("reload_workspace")) {
        reloadWorkspaceFromDisk();
    } else if (commandId == QStringLiteral("new_session")) {
        newSession();
    } else if (commandId == QStringLiteral("duplicate_session")) {
        duplicateSession(m_sessions ? m_sessions->activeIndex() : -1);
    } else if (commandId == QStringLiteral("delete_session")) {
        deleteSession(m_sessions ? m_sessions->activeIndex() : -1);
    } else if (commandId == QStringLiteral("duplicate_cell")) {
        duplicateCurrentCell();
    } else if (commandId == QStringLiteral("delete_cell")) {
        deleteCurrentCell();
    } else if (commandId == QStringLiteral("move_cell_up")) {
        moveCurrentCellUp();
    } else if (commandId == QStringLiteral("move_cell_down")) {
        moveCurrentCellDown();
    } else if (commandId == QStringLiteral("insert_cell_above")) {
        insertCellAbove();
    } else if (commandId == QStringLiteral("insert_cell_below")) {
        insertCellBelow();
    } else if (commandId == QStringLiteral("split_cell")) {
        splitCurrentCell();
    } else if (commandId == QStringLiteral("merge_with_previous")) {
        mergeWithPreviousCell();
    } else if (commandId == QStringLiteral("merge_with_next")) {
        mergeWithNextCell();
    } else if (commandId == QStringLiteral("select_previous_cell")) {
        selectPreviousCell();
    } else if (commandId == QStringLiteral("select_next_cell")) {
        selectNextCell();
    } else if (commandId == QStringLiteral("toggle_sidebar")) {
        toggleSidebar();
    } else if (commandId == QStringLiteral("toggle_inspector")) {
        toggleInspector();
    } else if (commandId == QStringLiteral("reset_layout")) {
        resetLayoutState();
    } else if (commandId == QStringLiteral("undo") ||
               commandId == QStringLiteral("undo_delete_cell") ||
               commandId == QStringLiteral("undo_delete_session")) {
        undo();
    } else if (commandId == QStringLiteral("redo")) {
        redo();
    }

    if (m_palette) {
        m_palette->recordInvocation(commandId);
        m_palette->setVisible(false);
    }
}

void AppCore::ensureSelectedRepresentationStillValid() {
    const QVariantList alternatives = selectedCellAlternatives();
    if (alternatives.isEmpty()) {
        if (!m_selectedRepresentationId.isEmpty()) {
            m_selectedRepresentationId.clear();
            emit selectedRepresentationChanged();
        }
        return;
    }

    for (const auto& item : alternatives) {
        const auto map = item.toMap();
        if (map.value(QStringLiteral("id")).toString() == m_selectedRepresentationId) {
            return;
        }
    }

    const QString fallback = alternatives.front().toMap().value(QStringLiteral("id")).toString();
    if (m_selectedRepresentationId != fallback) {
        m_selectedRepresentationId = fallback;
        emit selectedRepresentationChanged();
    }
}

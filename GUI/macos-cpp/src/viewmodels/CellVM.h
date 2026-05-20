#pragma once
#include <QObject>
#include <QString>
#include <QVariantList>
#include "../../../src/CasGuiSession.hpp"

namespace cas::gui {

class CellVM : public QObject {
    Q_OBJECT
    Q_PROPERTY(int index READ index NOTIFY indexChanged)
    Q_PROPERTY(QString inputLatex READ inputLatex WRITE setInputLatex NOTIFY inputChanged)
    Q_PROPERTY(QString statusText READ statusText NOTIFY statusChanged)
    Q_PROPERTY(bool hasError READ hasError NOTIFY statusChanged)
    Q_PROPERTY(bool hasOutput READ hasOutput NOTIFY outputChanged)
    Q_PROPERTY(QString outputLatex READ outputLatex NOTIFY outputChanged)
    Q_PROPERTY(QString outputMeta READ outputMeta NOTIFY outputChanged)
    Q_PROPERTY(bool active READ active WRITE setActive NOTIFY activeChanged)
    Q_PROPERTY(bool focused READ focused WRITE setFocused NOTIFY focusChanged)
    Q_PROPERTY(bool editing READ editing WRITE setEditing NOTIFY editingChanged)
    Q_PROPERTY(QString interactionState READ interactionState NOTIFY interactionStateChanged)
    Q_PROPERTY(bool restoredFromWorkspace READ restoredFromWorkspace NOTIFY restoredStateChanged)
    Q_PROPERTY(QVariantList alternatives READ alternatives NOTIFY outputChanged)
    Q_PROPERTY(QVariantList steps READ steps NOTIFY outputChanged)

public:
    explicit CellVM(int idx, QObject* parent = nullptr) : QObject(parent), m_index(idx) {}

    int index() const { return m_index; }
    void setIndex(int value) {
        if (m_index == value) return;
        m_index = value;
        emit indexChanged();
    }
    QString inputLatex() const { return m_input; }
    void setInputLatex(const QString& s) {
        if(m_input == s) return;
        emit inputAboutToChange();
        m_input = s;
        setRestoredFromWorkspace(false);
        if (m_hasError) {
            m_hasError = false;
            m_status = "ready";
            emit statusChanged();
        }
        emit inputChanged();
    }

    QString statusText() const { return m_status; }
    bool hasError() const { return m_hasError; }
    bool hasOutput() const { return m_hasOutput; }
    QString outputLatex() const { return m_output; }
    QString outputMeta() const { return m_meta; }
    bool active() const { return m_active; }
    void setActive(bool a) {
        if (m_active == a) return;
        m_active = a;
        emit activeChanged();
        emit interactionStateChanged();
    }
    bool focused() const { return m_focused; }
    void setFocused(bool value) {
        if (m_focused == value) return;
        m_focused = value;
        emit focusChanged();
        emit interactionStateChanged();
    }
    bool editing() const { return m_editing; }
    void setEditing(bool value) {
        if (m_editing == value) return;
        m_editing = value;
        emit editingChanged();
        emit interactionStateChanged();
    }
    QString interactionState() const {
        if (m_editing) {
            return QStringLiteral("editing");
        }
        if (m_focused) {
            return QStringLiteral("focused");
        }
        if (m_active) {
            return QStringLiteral("selected");
        }
        return QStringLiteral("idle");
    }
    bool restoredFromWorkspace() const { return m_restoredFromWorkspace; }
    void setRestoredFromWorkspace(bool value) {
        if (m_restoredFromWorkspace == value) return;
        m_restoredFromWorkspace = value;
        emit restoredStateChanged();
    }
    QVariantList alternatives() const { return m_alternatives; }
    QVariantList steps() const { return m_steps; }
    void clearOutput() {
        m_status = "ready";
        m_output.clear();
        m_meta.clear();
        m_hasOutput = false;
        m_alternatives.clear();
        m_steps.clear();
        emit statusChanged();
        emit outputChanged();
    }
    void setOutputState(const QString& status,
                        bool has_output,
                        const QString& output,
                        const QString& meta,
                        const QVariantList& alternatives = {},
                        const QVariantList& steps = {},
                        bool restored = false) {
        m_status = status;
        m_hasOutput = has_output;
        m_output = output;
        m_meta = meta;
        m_alternatives = alternatives;
        m_steps = steps;
        setRestoredFromWorkspace(restored);
        emit statusChanged();
        emit outputChanged();
    }

    void setResult(const ComputeResult& res) {
        setRestoredFromWorkspace(false);
        m_hasOutput = true;
        m_hasError = !res.ok;
        m_status = res.ok ? "ok" : QString::fromStdString(res.error);
        m_output = QString::fromStdString(res.latex);
        m_meta = res.ok
            ? QString("Computation complete · %1 step(s)").arg(res.steps.size())
            : "Error";
        m_alternatives.clear();
        for (const auto& representation : res.representations) {
            QVariantMap item;
            item["id"] = QString::fromStdString(representation.id);
            item["label"] = QString::fromStdString(representation.label);
            item["value"] = QString::fromStdString(representation.value);
            m_alternatives.push_back(item);
        }
        m_steps.clear();
        for (const auto& step : res.steps) {
            QVariantMap item;
            item["ruleId"] = static_cast<int>(step.rule_id);
            item["depth"] = static_cast<int>(step.depth);
            item["beforeLatex"] = QString::fromStdString(step.before_latex);
            item["afterLatex"] = QString::fromStdString(step.after_latex);
            item["rootLatex"] = QString::fromStdString(step.root_latex);
            m_steps.push_back(item);
        }
        emit statusChanged();
        emit outputChanged();
    }

signals:
    void indexChanged();
    void inputAboutToChange();
    void inputChanged();
    void statusChanged();
    void outputChanged();
    void activeChanged();
    void focusChanged();
    void editingChanged();
    void interactionStateChanged();
    void restoredStateChanged();

private:
    int m_index;
    QString m_input;
    QString m_status = "ready";
    bool m_hasError = false;
    QString m_output;
    QString m_meta;
    bool m_hasOutput = false;
    bool m_active = false;
    bool m_focused = false;
    bool m_editing = false;
    bool m_restoredFromWorkspace = false;
    QVariantList m_alternatives;
    QVariantList m_steps;
};

} // namespace cas::gui

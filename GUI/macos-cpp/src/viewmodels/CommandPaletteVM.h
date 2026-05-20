#pragma once
#include <QObject>
#include <QString>
#include <QStringList>
#include <QVariantList>
#include <QVector>
#include <algorithm>

namespace cas::gui {

class CommandPaletteVM : public QObject {
    Q_OBJECT
    Q_PROPERTY(bool visible READ visible WRITE setVisible NOTIFY visibleChanged)
    Q_PROPERTY(QString query READ query WRITE setQuery NOTIFY queryChanged)
    Q_PROPERTY(QVariantList commands READ commands CONSTANT)
    Q_PROPERTY(QVariantList filteredCommands READ filteredCommands NOTIFY filteredCommandsChanged)

public:
    explicit CommandPaletteVM(QObject* parent = nullptr) : QObject(parent) {}

    bool visible() const { return m_visible; }
    void setVisible(bool v) {
        if (v == m_visible) return;
        m_visible = v;
        emit visibleChanged();
    }

    QString query() const { return m_query; }
    void setQuery(const QString& value) {
        if (m_query == value) return;
        m_query = value;
        emit queryChanged();
        emit filteredCommandsChanged();
    }

    Q_INVOKABLE void toggle() { setVisible(!m_visible); }

    Q_INVOKABLE void recordInvocation(const QString& commandId) {
        if (commandId.isEmpty()) {
            return;
        }

        m_recentCommandIds.removeAll(commandId);
        m_recentCommandIds.prepend(commandId);
        while (m_recentCommandIds.size() > 6) {
            m_recentCommandIds.removeLast();
        }
        emit filteredCommandsChanged();
    }

    QVariantList commands() const {
        QVariantList list;
        for (const auto& command : catalog()) {
            list.push_back(command.toVariant());
        }
        return list;
    }

    QVariantList filteredCommands() const {
        const QString normalized_query = normalize(m_query);

        struct RankedCommand {
            int score = 0;
            CommandEntry command;
        };

        QVector<RankedCommand> ranked;
        ranked.reserve(catalog().size());

        for (const auto& command : catalog()) {
            const int score = scoreCommand(command, normalized_query);
            if (!normalized_query.isEmpty() && score <= 0) {
                continue;
            }
            ranked.push_back(RankedCommand{score, command});
        }

        std::stable_sort(ranked.begin(), ranked.end(), [](const RankedCommand& lhs, const RankedCommand& rhs) {
            if (lhs.score != rhs.score) {
                return lhs.score > rhs.score;
            }
            if (lhs.command.category != rhs.command.category) {
                return lhs.command.category < rhs.command.category;
            }
            if (lhs.command.priority != rhs.command.priority) {
                return lhs.command.priority < rhs.command.priority;
            }
            return lhs.command.id < rhs.command.id;
        });

        QVariantList result;
        result.reserve(ranked.size());
        for (const auto& item : ranked) {
            QVariantMap map = item.command.toVariant();
            map.insert(QStringLiteral("score"), item.score);
            result.push_back(map);
        }
        return result;
    }

    Q_INVOKABLE QString firstMatchingCommandId() const {
        const QVariantList filtered = filteredCommands();
        if (filtered.isEmpty()) {
            return {};
        }
        return filtered.front().toMap().value(QStringLiteral("id")).toString();
    }

signals:
    void visibleChanged();
    void queryChanged();
    void filteredCommandsChanged();

private:
    struct CommandEntry {
        QString id;
        QString label;
        QString shortcut;
        QString category;
        QStringList keywords;
        int priority = 0;

        QVariantMap toVariant() const {
            return QVariantMap{
                {QStringLiteral("id"), id},
                {QStringLiteral("label"), label},
                {QStringLiteral("shortcut"), shortcut},
                {QStringLiteral("category"), category},
                {QStringLiteral("keywords"), keywords},
            };
        }
    };

    static const QVector<CommandEntry>& catalog() {
        static const QVector<CommandEntry> entries = {
            {QStringLiteral("run"), QStringLiteral("Esegui cella"), QStringLiteral("⇧↵"), QStringLiteral("Notebook"),
             {QStringLiteral("execute"), QStringLiteral("compute"), QStringLiteral("eval")}, 10},
            {QStringLiteral("run_and_advance"), QStringLiteral("Esegui e vai alla prossima"), QStringLiteral("⌘↵"), QStringLiteral("Notebook"),
             {QStringLiteral("next"), QStringLiteral("advance"), QStringLiteral("compute")}, 11},
            {QStringLiteral("plot"), QStringLiteral("Plot cella"), QStringLiteral("⌘⇧P"), QStringLiteral("Plot"),
             {QStringLiteral("graph"), QStringLiteral("sample"), QStringLiteral("chart")}, 20},
            {QStringLiteral("reload_workspace"), QStringLiteral("Ricarica workspace da disco"), QStringLiteral("⌘⇧R"), QStringLiteral("Session"),
             {QStringLiteral("reload"), QStringLiteral("restore"), QStringLiteral("disk")}, 30},
            {QStringLiteral("new_session"), QStringLiteral("Nuova sessione"), QStringLiteral("⌘N"), QStringLiteral("Session"),
             {QStringLiteral("workspace"), QStringLiteral("new")}, 31},
            {QStringLiteral("duplicate_session"), QStringLiteral("Duplica sessione"), QStringLiteral("⌘D"), QStringLiteral("Session"),
             {QStringLiteral("copy"), QStringLiteral("workspace")}, 32},
            {QStringLiteral("delete_session"), QStringLiteral("Elimina sessione"), QStringLiteral("⌘⌫"), QStringLiteral("Session"),
             {QStringLiteral("remove"), QStringLiteral("workspace")}, 33},
            {QStringLiteral("insert_cell_above"), QStringLiteral("Inserisci cella sopra"), QStringLiteral("⌘⌥↑"), QStringLiteral("Notebook"),
             {QStringLiteral("insert"), QStringLiteral("above"), QStringLiteral("new cell")}, 12},
            {QStringLiteral("insert_cell_below"), QStringLiteral("Inserisci cella sotto"), QStringLiteral("⌘⌥↓"), QStringLiteral("Notebook"),
             {QStringLiteral("insert"), QStringLiteral("below"), QStringLiteral("new cell")}, 13},
            {QStringLiteral("duplicate_cell"), QStringLiteral("Duplica cella"), QStringLiteral("⌘⇧D"), QStringLiteral("Notebook"),
             {QStringLiteral("copy"), QStringLiteral("clone")}, 14},
            {QStringLiteral("delete_cell"), QStringLiteral("Elimina cella"), QStringLiteral("⌘⌦"), QStringLiteral("Notebook"),
             {QStringLiteral("remove"), QStringLiteral("trash")}, 15},
            {QStringLiteral("split_cell"), QStringLiteral("Dividi cella"), QStringLiteral("⌘⇧\\"), QStringLiteral("Notebook"),
             {QStringLiteral("split"), QStringLiteral("cursor")}, 16},
            {QStringLiteral("merge_with_previous"), QStringLiteral("Unisci con la precedente"), QStringLiteral("⌘⌥←"), QStringLiteral("Notebook"),
             {QStringLiteral("merge"), QStringLiteral("previous"), QStringLiteral("join")}, 17},
            {QStringLiteral("merge_with_next"), QStringLiteral("Unisci con la successiva"), QStringLiteral("⌘⌥→"), QStringLiteral("Notebook"),
             {QStringLiteral("merge"), QStringLiteral("next"), QStringLiteral("join")}, 18},
            {QStringLiteral("select_previous_cell"), QStringLiteral("Seleziona cella precedente"), QStringLiteral("⌃↑"), QStringLiteral("Notebook"),
             {QStringLiteral("navigate"), QStringLiteral("previous"), QStringLiteral("focus")}, 19},
            {QStringLiteral("select_next_cell"), QStringLiteral("Seleziona cella successiva"), QStringLiteral("⌃↓"), QStringLiteral("Notebook"),
             {QStringLiteral("navigate"), QStringLiteral("next"), QStringLiteral("focus")}, 20},
            {QStringLiteral("move_cell_up"), QStringLiteral("Sposta cella su"), QStringLiteral("⌥↑"), QStringLiteral("Notebook"),
             {QStringLiteral("move"), QStringLiteral("reorder"), QStringLiteral("up")}, 21},
            {QStringLiteral("move_cell_down"), QStringLiteral("Sposta cella giu"), QStringLiteral("⌥↓"), QStringLiteral("Notebook"),
             {QStringLiteral("move"), QStringLiteral("reorder"), QStringLiteral("down")}, 22},
            {QStringLiteral("toggle_sidebar"), QStringLiteral("Mostra o nascondi sidebar"), QStringLiteral("⌘⌥S"), QStringLiteral("Layout"),
             {QStringLiteral("layout"), QStringLiteral("sidebar"), QStringLiteral("panel")}, 40},
            {QStringLiteral("toggle_inspector"), QStringLiteral("Mostra o nascondi inspector"), QStringLiteral("⌘⌥I"), QStringLiteral("Layout"),
             {QStringLiteral("layout"), QStringLiteral("inspector"), QStringLiteral("panel")}, 41},
            {QStringLiteral("reset_layout"), QStringLiteral("Ripristina layout workspace"), QStringLiteral("⌘⌥0"), QStringLiteral("Layout"),
             {QStringLiteral("layout"), QStringLiteral("reset"), QStringLiteral("workspace")}, 42},
            {QStringLiteral("undo"), QStringLiteral("Annulla ultima azione"), QStringLiteral("⌘Z"), QStringLiteral("History"),
             {QStringLiteral("history"), QStringLiteral("rollback")}, 50},
            {QStringLiteral("redo"), QStringLiteral("Ripeti azione"), QStringLiteral("⌘⇧Z"), QStringLiteral("History"),
             {QStringLiteral("history"), QStringLiteral("repeat")}, 51},
        };
        return entries;
    }

    QString normalizedRecentAt(int index) const {
        if (index < 0 || index >= m_recentCommandIds.size()) {
            return {};
        }
        return m_recentCommandIds.at(index);
    }

    static QString normalize(const QString& value) {
        return value.trimmed().toLower();
    }

    int recentBonus(const QString& commandId) const {
        for (int i = 0; i < m_recentCommandIds.size(); ++i) {
            if (m_recentCommandIds.at(i) == commandId) {
                return std::max(0, 40 - (i * 6));
            }
        }
        return 0;
    }

    static int scoreField(const QString& normalized_field, const QString& query) {
        if (query.isEmpty()) {
            return 0;
        }
        if (normalized_field == query) {
            return 220;
        }
        if (normalized_field.startsWith(query)) {
            const int delta = static_cast<int>(normalized_field.size() - query.size());
            return 180 - std::min(40, delta);
        }
        const int contains_index = normalized_field.indexOf(query);
        if (contains_index >= 0) {
            return 120 - std::min(60, contains_index);
        }

        int query_index = 0;
        int matched = 0;
        for (const QChar ch : normalized_field) {
            if (query_index < query.size() && ch == query.at(query_index)) {
                ++query_index;
                ++matched;
            }
        }
        if (matched == query.size()) {
            const int delta = static_cast<int>(normalized_field.size() - query.size());
            return 70 - std::min(20, delta);
        }
        return 0;
    }

    int scoreCommand(const CommandEntry& command, const QString& query) const {
        int score = 0;
        if (query.isEmpty()) {
            score = 100 - command.priority;
        } else {
            score = std::max(score, scoreField(normalize(command.label), query));
            score = std::max(score, scoreField(normalize(command.id), query) - 5);
            score = std::max(score, scoreField(normalize(command.category), query) - 10);
            for (const auto& keyword : command.keywords) {
                score = std::max(score, scoreField(normalize(keyword), query) - 4);
            }
        }

        score += recentBonus(command.id);
        return score;
    }

    bool m_visible = false;
    QString m_query;
    QStringList m_recentCommandIds;
};

} // namespace cas::gui

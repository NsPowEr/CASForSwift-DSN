#pragma once

#include <QObject>
#include <QStringList>
#include <QVariantList>
#include <QtQml>

namespace cas::gui {

class AutocompleteVM : public QObject {
    Q_OBJECT
    QML_ELEMENT
    Q_PROPERTY(QString prefix READ prefix WRITE setPrefix NOTIFY prefixChanged)
    Q_PROPERTY(QVariantList suggestions READ suggestions NOTIFY suggestionsChanged)
    Q_PROPERTY(bool visible READ visible NOTIFY visibleChanged)

public:
    explicit AutocompleteVM(QObject* parent = nullptr);

    QString prefix() const { return m_prefix; }
    void setPrefix(const QString& p);

    QVariantList suggestions() const { return m_suggestions; }
    bool visible() const { return !m_suggestions.isEmpty() && !m_prefix.isEmpty(); }

    Q_INVOKABLE void updateLibrary(const QStringList& functions, const QVariantList& variables);

signals:
    void prefixChanged();
    void suggestionsChanged();
    void visibleChanged();

private:
    void refreshSuggestions();
    int fuzzyScore(const QString& source, const QString& target) const;
    QString normalizedLookupKey(const QString& value) const;

    QString m_prefix;
    QVariantList m_suggestions;
    QStringList m_functions;
    QVariantList m_variables;
    
    struct LibraryItem {
        QString text;
        QString type;
        QString insertText;
        QString detail;
    };
    std::vector<LibraryItem> m_library;
};

} // namespace cas::gui

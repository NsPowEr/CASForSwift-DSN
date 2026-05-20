#pragma once
#include <QObject>
#include <QStringList>
#include <QVariantList>

namespace cas::gui {

class SessionListVM : public QObject {
    Q_OBJECT
    Q_PROPERTY(QStringList titles READ titles NOTIFY titlesChanged)
    Q_PROPERTY(QVariantList items READ items NOTIFY itemsChanged)
    Q_PROPERTY(int activeIndex READ activeIndex WRITE setActiveIndex NOTIFY activeIndexChanged)

public:
    explicit SessionListVM(QObject* parent = nullptr) : QObject(parent) {}

    QStringList titles() const { return m_titles; }
    QVariantList items() const { return m_items; }
    int activeIndex() const { return m_activeIndex; }

    void setActiveIndex(int index) {
        if (index == m_activeIndex || index < 0 || index >= m_titles.size()) return;
        m_activeIndex = index;
        emit activeIndexChanged();
    }
    void setTitles(const QStringList& titles) {
        if (m_titles == titles) return;
        m_titles = titles;
        if (m_activeIndex >= m_titles.size()) {
            m_activeIndex = m_titles.isEmpty() ? -1 : m_titles.size() - 1;
            emit activeIndexChanged();
        }
        emit titlesChanged();
    }

    void setItems(const QVariantList& items) {
        if (m_items == items) {
            return;
        }
        m_items = items;
        emit itemsChanged();
    }

    Q_INVOKABLE void addSession(const QString& title) {
        m_titles.append(title);
        emit titlesChanged();
        setActiveIndex(m_titles.size() - 1);
    }

signals:
    void itemsChanged();
    void titlesChanged();
    void activeIndexChanged();

private:
    QVariantList m_items;
    QStringList m_titles{ "Sessione 1" };
    int m_activeIndex = 0;
};

} // namespace cas::gui

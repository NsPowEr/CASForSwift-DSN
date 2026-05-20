#pragma once
#include <QObject>
#include <QString>
#include <QNetworkAccessManager>
#include <QWebSocket>
#include "cas/result.hpp"

namespace cas::gui {

/**
 * @brief Client per il kernel remoto (Backend).
 * Gestisce REST per calcoli brevi e WebSocket per streaming di step/plot.
 */
class RemoteKernel : public QObject {
    Q_OBJECT
public:
    explicit RemoteKernel(const QString& url, QObject* parent = nullptr);

    void evaluate(const QString& input);

signals:
    void resultReady(const QString& result);
    void errorOccurred(const QString& error);

private:
    QString m_baseUrl;
    QNetworkAccessManager m_network;
    QWebSocket m_socket;
};

} // namespace cas::gui

#pragma once
#include <QObject>
#include <QString>
#include <QtConcurrent>
#include <functional>
#include <mutex>
#include "CasGuiSession.hpp"

namespace cas::gui {

class KernelDispatcher : public QObject {
    Q_OBJECT
public:
    explicit KernelDispatcher(QObject* parent = nullptr) 
        : QObject(parent) {}
    
    ComputeResult evaluate(const QString& input) {
        std::lock_guard<std::mutex> lock(mutex_);
        return m_session.simplify(input.toStdString());
    }

    void evaluateAsync(const QString& input, std::function<void(ComputeResult)> callback) {
        (void)QtConcurrent::run([this, input, callback]() {
            ComputeResult res;
            {
                std::lock_guard<std::mutex> lock(mutex_);
                res = m_session.simplify(input.toStdString());
            }
            callback(res);
        });
    }

    QString mode() const { return "embedded"; }
    QString version() const { return "session-bridge"; }

    QStringList listFunctions() const {
        QStringList result;
        std::lock_guard<std::mutex> lock(mutex_);
        for (const auto& fn : m_session.list_functions()) {
            result.append(QString::fromStdString(fn));
        }
        return result;
    }

    QVariantList listVariables() const {
        QVariantList result;
        std::lock_guard<std::mutex> lock(mutex_);
        for (const auto& [name, value] : m_session.list_variables()) {
            QVariantMap map;
            map["name"] = QString::fromStdString(name);
            map["value"] = QString::fromStdString(value);
            result.append(map);
        }
        return result;
    }

    std::vector<CasGuiSession::StoredDefinition> snapshotDefinitions() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return m_session.snapshot_definitions();
    }

    Result<void> restoreDefinitions(const std::vector<CasGuiSession::StoredDefinition>& definitions) {
        std::lock_guard<std::mutex> lock(mutex_);
        return m_session.restore_definitions(definitions);
    }

    void interrupt() {
        m_session.interrupt();
    }

    Result<std::vector<PlotSample>> sample2D(const QString& input,
                                             const QString& variable,
                                             double xMin,
                                             double xMax) {
        std::lock_guard<std::mutex> lock(mutex_);
        return m_session.sample_2d(
            input.toStdString(),
            variable.toStdString(),
            xMin,
            xMax);
    }

    void sample2DAsync(const QString& input, const QString& variable, double xMin, double xMax, std::function<void(Result<std::vector<PlotSample>>)> callback) {
        (void)QtConcurrent::run([this, input, variable, xMin, xMax, callback]() {
            auto res = [&]() -> Result<std::vector<PlotSample>> {
                std::lock_guard<std::mutex> lock(mutex_);
                return m_session.sample_2d(
                    input.toStdString(),
                    variable.toStdString(),
                    xMin,
                    xMax);
            }();
            callback(res);
        });
    }

private:
    mutable std::mutex mutex_;
    CasGuiSession m_session;
};

} // namespace cas::gui

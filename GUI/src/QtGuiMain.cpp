#include "CasGuiSession.hpp"

#include <QGuiApplication>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QQmlApplicationEngine>
#include <QQmlContext>
#include <QObject>
#include <QString>

namespace {

class GuiBackend final : public QObject {
    Q_OBJECT
    Q_PROPERTY(QString input READ input WRITE setInput NOTIFY inputChanged)
    Q_PROPERTY(QString status READ status NOTIFY outputChanged)
    Q_PROPERTY(QString textResult READ textResult NOTIFY outputChanged)
    Q_PROPERTY(QString latexResult READ latexResult NOTIFY outputChanged)
    Q_PROPERTY(QString asciiResult READ asciiResult NOTIFY outputChanged)
    Q_PROPERTY(QString numericResult READ numericResult NOTIFY outputChanged)
    Q_PROPERTY(QString plotJson READ plotJson NOTIFY plotChanged)

public:
    QString input() const { return input_; }
    QString status() const { return status_; }
    QString textResult() const { return text_result_; }
    QString latexResult() const { return latex_result_; }
    QString asciiResult() const { return ascii_result_; }
    QString numericResult() const { return numeric_result_; }
    QString plotJson() const { return plot_json_; }

    void setInput(const QString& value) {
        if (input_ == value) {
            return;
        }
        input_ = value;
        emit inputChanged();
    }

    Q_INVOKABLE void compute() {
        const auto result = session_.simplify(input_.toStdString());
        if (!result.ok) {
            status_ = QString::fromStdString(result.error);
            text_result_.clear();
            latex_result_.clear();
            ascii_result_.clear();
            numeric_result_.clear();
            emit outputChanged();
            return;
        }

        status_ = QStringLiteral("ok");
        text_result_ = QString::fromStdString(result.text);
        latex_result_ = QString::fromStdString(result.latex);
        ascii_result_ = QString::fromStdString(result.ascii);
        numeric_result_ = result.numeric_value.has_value()
            ? QString::number(*result.numeric_value, 'g', 16)
            : QStringLiteral("not numeric");
        emit outputChanged();
    }

    Q_INVOKABLE void sample2d(const QString& variable, double xMin, double xMax) {
        auto sampled = session_.sample_2d(input_.toStdString(), variable.toStdString(), xMin, xMax);
        if (sampled.is_error()) {
            status_ = QString::fromStdString(sampled.error().message);
            emit outputChanged();
            return;
        }

        QJsonArray points;
        for (const auto& point : sampled.value()) {
            QJsonObject item;
            item.insert(QStringLiteral("x"), point.x);
            item.insert(QStringLiteral("y"), point.y);
            points.append(item);
        }

        plot_json_ = QString::fromUtf8(QJsonDocument(points).toJson(QJsonDocument::Compact));
        emit plotChanged();
    }

signals:
    void inputChanged();
    void outputChanged();
    void plotChanged();

private:
    cas::gui::CasGuiSession session_;
    QString input_{QStringLiteral("sin(x)")};
    QString status_{QStringLiteral("ready")};
    QString text_result_;
    QString latex_result_;
    QString ascii_result_;
    QString numeric_result_;
    QString plot_json_{QStringLiteral("[]")};
};

} // namespace

int main(int argc, char* argv[]) {
    QGuiApplication app(argc, argv);

    GuiBackend backend;
    QQmlApplicationEngine engine;
    engine.rootContext()->setContextProperty(QStringLiteral("backend"), &backend);
    engine.loadFromModule(QStringLiteral("CASGui"), QStringLiteral("Main"));

    if (engine.rootObjects().isEmpty()) {
        return 1;
    }

    return QGuiApplication::exec();
}

#include "QtGuiMain.moc"

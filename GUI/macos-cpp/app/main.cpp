// main.cpp — bootstrap dell'applicazione.

#include <QGuiApplication>
#include <QFont>
#include <QQmlApplicationEngine>
#include <QQuickStyle>
#include <QIcon>
#include "viewmodels/AppCore.h"

int main(int argc, char* argv[]) {
    QGuiApplication app(argc, argv);
    app.setOrganizationName("CASCalc");
    app.setOrganizationDomain("cascalc.app");
    app.setApplicationName("CAS Calculator");
    app.setApplicationVersion("0.1.0");
    app.setWindowIcon(QIcon(":/resources/icon.icns"));
    app.setFont(QFont("Helvetica Neue", 13));

    QQuickStyle::setStyle("macOS");

    QQmlApplicationEngine engine;

    // AppCore è registrato come singleton via QML_SINGLETON in AppCore.h
    engine.load("qrc:/qt/qml/CAS/macos-cpp/qml/Main.qml");
    if (engine.rootObjects().isEmpty()) return -1;

    return app.exec();
}

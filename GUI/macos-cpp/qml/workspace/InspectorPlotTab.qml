import QtQuick
import QtQuick.Layouts
import QtQuick.Controls
import CAS

Item {
    ScrollView {
        anchors.fill: parent
        contentWidth: parent.width

        ColumnLayout {
            width: parent.width
            spacing: 12
            anchors.margins: 12

            Label {
                text: "Plot Settings"
                font.bold: true
                font.family: AppCore.theme.fontUI
            }

            GridLayout {
                columns: 2
                columnSpacing: 8
                rowSpacing: 8
                Layout.fillWidth: true

                Label { text: "Variabile"; color: AppCore.theme.textMuted; font.family: AppCore.theme.fontUI }
                TextField {
                    Layout.fillWidth: true
                    text: AppCore.plot.variable
                    font.family: AppCore.theme.fontUI
                    onTextEdited: AppCore.plot.variable = text
                }

                Label { text: "xRange"; color: AppCore.theme.textMuted; font.family: AppCore.theme.fontUI }
                RowLayout {
                    Layout.fillWidth: true
                    TextField {
                        Layout.fillWidth: true
                        text: AppCore.plot.xMin.toFixed(2)
                        font.family: AppCore.theme.fontUI
                        onEditingFinished: AppCore.plot.setXMin(Number.fromLocaleString(Qt.locale(), text))
                    }
                    Label { text: ":" }
                    TextField {
                        Layout.fillWidth: true
                        text: AppCore.plot.xMax.toFixed(2)
                        font.family: AppCore.theme.fontUI
                        onEditingFinished: AppCore.plot.setXMax(Number.fromLocaleString(Qt.locale(), text))
                    }
                }
            }

            RowLayout {
                Layout.fillWidth: true
                spacing: 6
                Button { 
                    text: "Clear All"
                    font.pixelSize: 10
                    onClicked: AppCore.plot.clearSeries()
                }
                Item { Layout.fillWidth: true }
                Button { text: "Reset View"; font.pixelSize: 10; onClicked: AppCore.resetPlotRange() }
            }

            Rectangle {
                Layout.fillWidth: true
                Layout.preferredHeight: 240
                color: AppCore.theme.bgCard
                radius: AppCore.theme.radiusM
                border.color: AppCore.theme.stroke
                clip: true

                PlotCanvas {
                    id: canvas
                    anchors.fill: parent
                    anchors.margins: 12
                    visible: AppCore.plot.hasPoints
                    seriesList: AppCore.plot.seriesList
                    xMin: AppCore.plot.xMin
                    xMax: AppCore.plot.xMax
                    yMin: AppCore.plot.yMin
                    yMax: AppCore.plot.yMax
                }

                MouseArea {
                    anchors.fill: canvas
                    visible: AppCore.plot.hasPoints
                    hoverEnabled: true
                    property real lastX: 0
                    onPressed: (mouse) => { lastX = mouse.x }
                    onPositionChanged: (mouse) => {
                        if (pressed) {
                            let dx = mouse.x - lastX
                            if (Math.abs(dx) > 0) {
                                let span = AppCore.plot.xMax - AppCore.plot.xMin
                                let shift = (dx / canvas.width) * span
                                AppCore.plot.setRange(AppCore.plot.xMin - shift, AppCore.plot.xMax - shift)
                                lastX = mouse.x
                            }
                        }
                    }
                    onWheel: (wheel) => {
                        let span = AppCore.plot.xMax - AppCore.plot.xMin
                        let factor = wheel.angleDelta.y > 0 ? 0.85 : 1.25
                        let mouseFrac = wheel.x / canvas.width
                        let mouseVal = AppCore.plot.xMin + mouseFrac * span
                        let newSpan = span * factor
                        AppCore.plot.setRange(mouseVal - mouseFrac * newSpan, mouseVal + (1 - mouseFrac) * newSpan)
                    }
                }

                Label {
                    anchors.centerIn: parent
                    visible: !AppCore.plot.hasPoints
                    text: AppCore.plot.status
                    color: AppCore.theme.textFaint
                    font.italic: true
                }
            }

            Label {
                text: "Series"
                font.bold: true
                font.family: AppCore.theme.fontUI
                visible: AppCore.plot.hasPoints
            }

            ListView {
                Layout.fillWidth: true
                Layout.preferredHeight: contentHeight
                model: AppCore.plot.seriesList
                interactive: false
                visible: AppCore.plot.hasPoints
                delegate: RowLayout {
                    width: parent.width
                    spacing: 8
                    Rectangle {
                        width: 12; height: 12; radius: 6
                        color: modelData.color
                    }
                    Label {
                        text: modelData.name
                        font.family: AppCore.theme.fontMono
                        font.pixelSize: 11
                        Layout.fillWidth: true
                        elide: Text.ElideRight
                    }
                }
            }

            Item { Layout.fillHeight: true }
        }
    }
}

import QtQuick
import QtQuick.Layouts
import QtQuick.Controls
import CAS

Item {
    ScrollView {
        anchors.fill: parent

        ColumnLayout {
            width: parent.width
            spacing: 10

            Label {
                visible: AppCore.selectedCellSteps.length === 0
                text: "Nessun passo disponibile."
                color: AppCore.theme.textFaint
            }

            Repeater {
                model: AppCore.selectedCellSteps

                delegate: Rectangle {
                    required property int index
                    required property var modelData

                    Layout.fillWidth: true
                    color: AppCore.theme.bgCard
                    border.color: AppCore.theme.stroke
                    radius: AppCore.theme.radiusM
                    implicitHeight: stepLayout.implicitHeight + 20

                    ColumnLayout {
                        id: stepLayout
                        anchors.fill: parent
                        anchors.margins: 10
                        spacing: 6

                        Label {
                            text: "#" + (index + 1)
                                  + " · rule " + modelData.ruleId
                                  + " · depth " + modelData.depth
                            color: AppCore.theme.textMuted
                            font.family: AppCore.theme.fontMono
                        }

                        MathView {
                            visible: modelData.rootLatex.length > 0
                            latex: modelData.rootLatex
                            fontSize: 14
                            color: AppCore.theme.textFaint
                            block: true
                        }

                        MathView {
                            visible: modelData.beforeLatex.length > 0
                            latex: modelData.beforeLatex
                            fontSize: 15
                            color: AppCore.theme.textFaint
                            block: true
                        }

                        MathView {
                            visible: modelData.afterLatex.length > 0
                            latex: modelData.afterLatex
                            fontSize: 17
                            color: AppCore.theme.text
                            block: true
                        }
                    }
                }
            }
        }
    }
}

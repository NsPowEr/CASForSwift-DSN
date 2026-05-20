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

            TextField {
                id: variableSearchField
                objectName: "variableSearchField"
                Layout.fillWidth: true
                placeholderText: "Filter variables by name or value"
                text: AppCore.variableSearchQuery
                font.family: AppCore.theme.fontUI
                onTextChanged: AppCore.variableSearchQuery = text
            }

            Label {
                text: "Session Context"
                font.bold: true
                font.family: AppCore.theme.fontUI
                color: AppCore.theme.textMuted
            }

            Rectangle {
                Layout.fillWidth: true
                implicitHeight: varColumn.implicitHeight + 24
                color: AppCore.theme.bgCard
                radius: AppCore.theme.radiusM
                border.color: AppCore.theme.stroke

                ColumnLayout {
                    id: varColumn
                    anchors.fill: parent
                    anchors.margins: 12
                    spacing: 8

                    Repeater {
                        model: AppCore.filteredVariables

                        delegate: ColumnLayout {
                            required property var modelData
                            Layout.fillWidth: true
                            spacing: 4

                            RowLayout {
                                spacing: 8
                                Rectangle {
                                    width: 8; height: 8; radius: 4
                                    color: AppCore.theme.accent
                                }
                                Label {
                                    text: modelData.name
                                    elide: Text.ElideRight
                                    font.bold: true
                                    font.family: AppCore.theme.fontMono
                                    color: AppCore.theme.text
                                    Layout.maximumWidth: 140
                                }
                                Label { text: ":"; color: AppCore.theme.textFaint }
                                Label {
                                    text: {
                                        let val = modelData.value
                                        if (val.includes("Matrix")) return "Matrix"
                                        if (val.includes("Poly")) return "Polynomial"
                                        return "Expression"
                                    }
                                    font.pixelSize: 10
                                    font.family: AppCore.theme.fontUI
                                    color: AppCore.theme.textFaint
                                }
                                Item { Layout.fillWidth: true }
                                Button {
                                    text: "Insert"
                                    font.pixelSize: 9
                                    flat: true
                                    onClicked: AppCore.insertVariableName(modelData.name)
                                }
                                Button {
                                    text: "Copy Name"
                                    font.pixelSize: 9
                                    flat: true
                                    onClicked: AppCore.copyVariableName(modelData.name)
                                }
                                Button {
                                    text: "Copy Value"
                                    font.pixelSize: 9
                                    flat: true
                                    onClicked: AppCore.copyVariableValue(modelData.name)
                                }
                            }

                            Rectangle {
                                Layout.fillWidth: true
                                height: mathView.implicitHeight + 16
                                color: Qt.rgba(1, 1, 1, 0.03)
                                radius: 4
                                clip: true

                                MathView {
                                    id: mathView
                                    anchors.centerIn: parent
                                    width: parent.width - 24
                                    latex: modelData.value
                                    fontSize: 14
                                    color: AppCore.theme.text
                                    block: false
                                }
                            }
                            
                            Rectangle {
                                Layout.fillWidth: true
                                height: 1
                                color: AppCore.theme.stroke
                                opacity: 0.5
                                visible: index < AppCore.variables.length - 1
                            }
                        }
                    }

                    Label {
                        visible: AppCore.filteredVariables.length === 0
                        text: AppCore.variables.length === 0
                              ? "No active definitions"
                              : "No variables match the current filter"
                        color: AppCore.theme.textFaint
                        font.italic: true
                        Layout.alignment: Qt.AlignCenter
                    }
                }
            }

            Label {
                text: "Built-in Functions"
                font.bold: true
                font.family: AppCore.theme.fontUI
                color: AppCore.theme.textMuted
                Layout.topMargin: 8
            }

            Flow {
                Layout.fillWidth: true
                spacing: 6
                Repeater {
                    model: AppCore.functions
                    delegate: Rectangle {
                        required property var modelData
                        width: fnLabel.implicitWidth + 16
                        height: 22
                        radius: 11
                        color: AppCore.theme.bgCard
                        border.color: AppCore.theme.stroke

                        Label {
                            id: fnLabel
                            anchors.centerIn: parent
                            text: modelData
                            font.pixelSize: 10
                            font.family: AppCore.theme.fontMono
                            color: AppCore.theme.textMuted
                        }
                    }
                }
            }

            Item { Layout.fillHeight: true }
        }
    }
}

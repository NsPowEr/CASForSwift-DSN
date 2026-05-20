import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import CAS

Rectangle {
    id: root
    color: AppCore.theme.bgCard
    radius: AppCore.theme.radiusL
    border.color: AppCore.theme.stroke
    visible: AppCore.palette.visible

    ColumnLayout {
        anchors.fill: parent
        anchors.margins: 12

        TextField {
            id: searchField
            Layout.fillWidth: true
            placeholderText: "Cerca comando..."
            font.family: AppCore.theme.fontUI
            focus: root.visible
            text: AppCore.palette.query
            onTextChanged: AppCore.palette.query = text
            onAccepted: {
                const commandId = AppCore.palette.firstMatchingCommandId()
                if (commandId.length > 0 && AppCore.isCommandEnabled(commandId)) {
                    AppCore.invokeCommand(commandId)
                } else {
                    AppCore.palette.visible = false
                }
            }
        }

        ListView {
            Layout.fillWidth: true
            Layout.fillHeight: true
            model: AppCore.palette.filteredCommands
            delegate: ItemDelegate {
                width: ListView.view ? ListView.view.width : 0
                font.family: AppCore.theme.fontUI
                enabled: AppCore.isCommandEnabled(modelData.id)
                height: commandLayout.implicitHeight + 14
                background: Rectangle {
                    radius: 8
                    color: parent.hovered
                        ? Qt.rgba(AppCore.theme.accent.r, AppCore.theme.accent.g, AppCore.theme.accent.b, 0.08)
                        : "transparent"
                    border.color: parent.highlighted
                        ? Qt.rgba(AppCore.theme.accent.r, AppCore.theme.accent.g, AppCore.theme.accent.b, 0.35)
                        : "transparent"
                }

                contentItem: ColumnLayout {
                    id: commandLayout
                    spacing: 3

                    RowLayout {
                        spacing: 8

                        Rectangle {
                            radius: 5
                            color: Qt.rgba(AppCore.theme.accent.r, AppCore.theme.accent.g, AppCore.theme.accent.b, 0.12)
                            implicitWidth: categoryLabel.implicitWidth + 12
                            implicitHeight: 18

                            Label {
                                id: categoryLabel
                                anchors.centerIn: parent
                                text: modelData.category
                                font.family: AppCore.theme.fontUI
                                font.pixelSize: 10
                                color: AppCore.theme.accent
                            }
                        }

                        Label {
                            Layout.fillWidth: true
                            text: AppCore.commandDisplayLabel(modelData.id)
                            font.family: AppCore.theme.fontUI
                            color: enabled ? AppCore.theme.text : AppCore.theme.textMuted
                            elide: Text.ElideRight
                        }

                        Label {
                            text: modelData.shortcut
                            font.family: AppCore.theme.fontMono
                            color: AppCore.theme.textFaint
                        }
                    }

                    Label {
                        visible: !enabled
                        text: AppCore.commandDisabledReason(modelData.id)
                        font.family: AppCore.theme.fontUI
                        font.pixelSize: 11
                        color: "#c96f1a"
                        wrapMode: Text.Wrap
                    }
                }

                onClicked: {
                    if (enabled) {
                        AppCore.invokeCommand(modelData.id)
                    }
                }
            }
        }
    }

    Keys.onEscapePressed: AppCore.palette.visible = false
}

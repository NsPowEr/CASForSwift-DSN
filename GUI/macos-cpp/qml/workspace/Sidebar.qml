import QtQuick
import QtQuick.Layouts
import QtQuick.Controls
import CAS

Rectangle {
    id: root
    color: AppCore.theme.bgRaised
    border.color: AppCore.theme.stroke
    border.width: 1
    property bool syncingTitle: false

    function refreshTitleEditor() {
        const idx = AppCore.sessions.activeIndex
        const items = AppCore.sessions.items
        syncingTitle = true
        titleField.text = (idx >= 0 && idx < items.length) ? items[idx].title : ""
        syncingTitle = false
    }

    Connections {
        target: AppCore.sessions
        function onActiveIndexChanged() { root.refreshTitleEditor() }
        function onTitlesChanged() { root.refreshTitleEditor() }
    }

    Component.onCompleted: refreshTitleEditor()

    ColumnLayout {
        anchors.fill: parent
        anchors.margins: 12
        spacing: 16

        Label {
            text: "Sessioni"
            font.bold: true
            font.family: AppCore.theme.fontUI
            color: AppCore.theme.textMuted

            Button {
                anchors.right: parent.right
                anchors.verticalCenter: parent.verticalCenter
                text: "Reset"
                font.pixelSize: 10
                flat: true
                onClicked: AppCore.reloadWorkspaceFromDisk()
            }
        }

        ListView {
            id: sessionList
            Layout.fillWidth: true
            Layout.fillHeight: true
            model: AppCore.sessions.items
            spacing: 4
            delegate: ItemDelegate {
                required property var modelData
                required property int index
                width: sessionList.width
                height: sessionCard.implicitHeight + 12
                font.family: AppCore.theme.fontUI
                highlighted: index === AppCore.sessions.activeIndex
                onClicked: AppCore.openSession(index)

                contentItem: ColumnLayout {
                    id: sessionCard
                    spacing: 2

                    RowLayout {
                        spacing: 8

                        Label {
                            text: modelData.title
                            font.family: AppCore.theme.fontUI
                            font.bold: highlighted
                            color: AppCore.theme.text
                            elide: Text.ElideRight
                            Layout.fillWidth: true
                        }

                        Rectangle {
                            visible: modelData.active
                            radius: 5
                            color: Qt.rgba(AppCore.theme.accent.r, AppCore.theme.accent.g, AppCore.theme.accent.b, 0.14)
                            implicitWidth: activeLabel.implicitWidth + 10
                            implicitHeight: 18

                            Label {
                                id: activeLabel
                                anchors.centerIn: parent
                                text: "active"
                                font.family: AppCore.theme.fontMono
                                font.pixelSize: 9
                                color: AppCore.theme.accent
                            }
                        }

                        Rectangle {
                            visible: modelData.selectedStatusLabel === "error"
                            radius: 5
                            color: Qt.rgba(AppCore.theme.error.r, AppCore.theme.error.g, AppCore.theme.error.b, 0.14)
                            implicitWidth: errorLabel.implicitWidth + 10
                            implicitHeight: 18

                            Label {
                                id: errorLabel
                                anchors.centerIn: parent
                                text: "error"
                                font.family: AppCore.theme.fontMono
                                font.pixelSize: 9
                                color: AppCore.theme.error
                            }
                        }
                    }

                    Label {
                        text: modelData.selectedPreview
                        font.family: AppCore.theme.fontMono
                        font.pixelSize: 11
                        color: AppCore.theme.textMuted
                        elide: Text.ElideRight
                    }

                    Label {
                        text: modelData.summary + " · cella " + modelData.selectedCellNumber + " · " + modelData.selectedStatusLabel
                        font.family: AppCore.theme.fontUI
                        font.pixelSize: 11
                        color: AppCore.theme.textFaint
                        elide: Text.ElideRight
                    }
                }

                background: Rectangle {
                    color: highlighted ? AppCore.theme.accentSoft : "transparent"
                    radius: AppCore.theme.radiusS
                }

                TapHandler {
                    acceptedButtons: Qt.RightButton
                    onTapped: sessionMenu.popup(index)
                }

                Menu {
                    id: sessionMenu
                    property int targetIndex: -1
                    function popup(idx) { targetIndex = idx; open() }

                    MenuItem {
                        text: "Rinomina..."
                        onTriggered: {
                            AppCore.openSession(sessionMenu.targetIndex)
                            titleField.forceActiveFocus()
                            titleField.selectAll()
                        }
                    }
                    MenuItem {
                        text: "Duplica"
                        onTriggered: AppCore.duplicateSession(sessionMenu.targetIndex)
                    }
                    MenuItem {
                        text: "Sposta su"
                        enabled: sessionMenu.targetIndex > 0
                        onTriggered: AppCore.moveSessionUp(sessionMenu.targetIndex)
                    }
                    MenuItem {
                        text: "Sposta giu"
                        enabled: sessionMenu.targetIndex < AppCore.sessions.items.length - 1
                        onTriggered: AppCore.moveSessionDown(sessionMenu.targetIndex)
                    }
                    MenuSeparator {}
                    MenuItem {
                        text: "Elimina"
                        enabled: AppCore.canDeleteSession
                        onTriggered: AppCore.deleteSession(sessionMenu.targetIndex)
                    }
                }
            }
        }

        Label {
            text: "Sessione attiva"
            color: AppCore.theme.textFaint
            font.family: AppCore.theme.fontUI
            font.pixelSize: 11
        }

        Label {
            visible: AppCore.sessions.activeIndex >= 0 && AppCore.sessions.activeIndex < AppCore.sessions.items.length
            text: {
                const idx = AppCore.sessions.activeIndex
                const item = AppCore.sessions.items[idx]
                return item.summary + " · cella " + item.selectedCellNumber + " · " + item.selectedStatusLabel
            }
            color: AppCore.theme.textMuted
            font.family: AppCore.theme.fontUI
            font.pixelSize: 11
            wrapMode: Text.Wrap
            Layout.fillWidth: true
        }

        Label {
            visible: AppCore.sessions.activeIndex >= 0 && AppCore.sessions.activeIndex < AppCore.sessions.items.length
            text: AppCore.sessions.items[AppCore.sessions.activeIndex].selectedPreview
            color: AppCore.theme.textFaint
            font.family: AppCore.theme.fontMono
            font.pixelSize: 11
            wrapMode: Text.Wrap
            Layout.fillWidth: true
        }

        TextField {
            id: titleField
            Layout.fillWidth: true
            placeholderText: "Titolo sessione"
            font.family: AppCore.theme.fontUI
            onEditingFinished: {
                if (!root.syncingTitle) {
                    AppCore.renameActiveSession(text)
                }
            }
        }

        RowLayout {
            Layout.fillWidth: true
            spacing: 8

            Button {
                text: "Su"
                font.family: AppCore.theme.fontUI
                enabled: AppCore.canMoveSessionUp
                onClicked: AppCore.moveSessionUp(AppCore.sessions.activeIndex)
            }

            Button {
                text: "Giu"
                font.family: AppCore.theme.fontUI
                enabled: AppCore.canMoveSessionDown
                onClicked: AppCore.moveSessionDown(AppCore.sessions.activeIndex)
            }

            Button {
                Layout.fillWidth: true
                text: "Duplica"
                font.family: AppCore.theme.fontUI
                onClicked: AppCore.duplicateSession(AppCore.sessions.activeIndex)
            }

            Button {
                Layout.fillWidth: true
                text: "Elimina"
                font.family: AppCore.theme.fontUI
                enabled: AppCore.canDeleteSession
                onClicked: AppCore.deleteSession(AppCore.sessions.activeIndex)
            }
        }

        RowLayout {
            Layout.fillWidth: true
            spacing: 8

            Button {
                Layout.fillWidth: true
                text: AppCore.undoActionLabel
                font.family: AppCore.theme.fontUI
                enabled: AppCore.canUndo
                onClicked: AppCore.undo()
            }

            Button {
                Layout.fillWidth: true
                text: AppCore.redoActionLabel
                font.family: AppCore.theme.fontUI
                enabled: AppCore.canRedo
                onClicked: AppCore.redo()
            }
        }

        Button {
            Layout.fillWidth: true
            text: "+ Nuova Sessione"
            font.family: AppCore.theme.fontUI
            onClicked: AppCore.newSession()
        }
    }
}

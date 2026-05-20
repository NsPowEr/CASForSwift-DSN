import QtQuick
import QtQuick.Layouts
import QtQuick.Controls
import CAS

Rectangle {
    id: root
    color: AppCore.theme.bg

    property alias sessionTitle: titleLabel.text

    ColumnLayout {
        anchors.fill: parent
        spacing: 0

        // ── Recovery / Warning Banner
        Rectangle {
            id: recoveryBanner
            Layout.fillWidth: true
            height: AppCore.storageNotification !== "" ? 44 : 0
            color: AppCore.storageSeverity >= 3 ? "#fdecea" : "#fff8e1"
            clip: true
            visible: height > 0
            
            Behavior on height { NumberAnimation { duration: 200 } }

            RowLayout {
                anchors.fill: parent
                anchors.leftMargin: 12
                anchors.rightMargin: 12
                spacing: 8

                Label {
                    text: "⚠️"
                    visible: AppCore.storageSeverity >= 2
                }

                Label {
                    Layout.fillWidth: true
                    text: AppCore.storageNotification
                    font.family: AppCore.theme.fontUI
                    font.pixelSize: 12
                    color: "#333"
                    elide: Text.ElideRight
                }

                Button {
                    text: "Chiudi"
                    flat: true
                    onClicked: AppCore.clearStorageNotification()
                }
            }

            Rectangle {
                anchors.bottom: parent.bottom
                width: parent.width; height: 1
                color: AppCore.storageSeverity >= 3 ? "#f5c6cb" : "#ffeeba"
            }
        }

        Rectangle {
            Layout.fillWidth: true
            height: 40
            color: "transparent"
            Label {
                id: titleLabel
                anchors.centerIn: parent
                text: AppCore.notebook.sessionTitle
                font.family: AppCore.theme.fontUI
                font.pixelSize: 14
                color: AppCore.theme.textMuted
            }
        }

        ListView {
            id: cellsList
            objectName: "cellsList"
            Layout.fillWidth: true
            Layout.fillHeight: true
            model: AppCore.notebook.cells
            spacing: 12
            clip: true
            focus: true
            boundsBehavior: Flickable.StopAtBounds

            onCurrentIndexChanged: {
                if (currentIndex !== AppCore.notebook.selectedIndex) {
                    AppCore.notebook.selectedIndex = currentIndex
                }
            }
            onContentYChanged: {
                if (Math.abs(AppCore.notebookScrollY - contentY) > 0.5) {
                    AppCore.notebookScrollY = contentY
                }
            }
            Connections {
                target: AppCore.notebook
                function onSelectedIndexChanged() {
                    if (cellsList.currentIndex !== AppCore.notebook.selectedIndex) {
                        cellsList.currentIndex = AppCore.notebook.selectedIndex
                    }
                    if (AppCore.notebook.selectedIndex >= 0) {
                        cellsList.positionViewAtIndex(AppCore.notebook.selectedIndex, ListView.Contain)
                    }
                    Qt.callLater(function() { cellsList.contentY = AppCore.notebookScrollY })
                }
                function onCellsChanged() {
                    Qt.callLater(function() { cellsList.contentY = AppCore.notebookScrollY })
                }
            }

            Connections {
                target: AppCore
                function onLayoutStateChanged() {
                    Qt.callLater(function() { cellsList.contentY = AppCore.notebookScrollY })
                }
            }
            
            delegate: CellView {
                required property var modelData
                width: cellsList.width - 24
                x: 12
                cell: modelData
            }

            footer: Item {
                width: parent.width
                height: 100
                RowLayout {
                    anchors.centerIn: parent
                    spacing: 10

                    Button {
                        text: "Sopra"
                        font.family: AppCore.theme.fontUI
                        onClicked: AppCore.insertCellAbove()
                    }
                    Button {
                        text: "Sotto"
                        font.family: AppCore.theme.fontUI
                        onClicked: AppCore.insertCellBelow()
                    }
                    Button {
                        text: "Aggiungi"
                        font.family: AppCore.theme.fontUI
                        onClicked: AppCore.addEmptyCell()
                    }
                }
            }
        }
    }
}

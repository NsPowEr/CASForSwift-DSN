import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import CAS

Rectangle {
    id: root
    width: 240
    height: Math.min(listView.contentHeight + 8, 300)
    color: AppCore.theme.bgCard
    border.color: AppCore.theme.stroke
    radius: AppCore.theme.radiusM
    visible: vm.visible
    z: 100

    property alias prefix: vm.prefix
    signal selected(string text, string insertText, string type)

    AutocompleteVM {
        id: vm
    }

    // Sync library from AppCore
    Connections {
        target: AppCore
        function onVariablesChanged() { vm.updateLibrary(AppCore.functions, AppCore.variables) }
    }
    Component.onCompleted: vm.updateLibrary(AppCore.functions, AppCore.variables)

    ListView {
        id: listView
        anchors.fill: parent
        anchors.margins: 4
        model: vm.suggestions
        clip: true
        currentIndex: 0

        delegate: ItemDelegate {
            width: listView.width
            height: 32
            contentItem: RowLayout {
                Label {
                    text: modelData.text
                    font.family: AppCore.theme.fontUI
                    color: AppCore.theme.text
                    Layout.fillWidth: true
                }
                Label {
                    text: modelData.detail || ""
                    font.family: AppCore.theme.fontUI
                    color: AppCore.theme.textMuted
                    font.pixelSize: 11
                    elide: Text.ElideRight
                }
                Label {
                    text: {
                        if (modelData.type === "function") return "ƒ"
                        if (modelData.type === "template") return "T"
                        return "x"
                    }
                    font.family: AppCore.theme.fontUI
                    color: AppCore.theme.textMuted
                    font.italic: true
                }
            }
            onClicked: {
                root.selected(modelData.text, modelData.insertText, modelData.type)
                vm.prefix = ""
            }
            highlighted: ListView.isCurrentItem
        }

        ScrollIndicator.vertical: ScrollIndicator {}
    }

    function moveUp() {
        if (listView.currentIndex > 0) listView.currentIndex--
    }

    function moveDown() {
        if (listView.currentIndex < listView.count - 1) listView.currentIndex++
    }

    function acceptCurrent() {
        if (vm.visible && listView.currentIndex >= 0 && listView.currentIndex < listView.count) {
            let item = vm.suggestions[listView.currentIndex]
            root.selected(item.text, item.insertText, item.type)
            vm.prefix = ""
            return true
        }
        return false
    }
}

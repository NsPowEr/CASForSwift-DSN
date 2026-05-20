import QtQuick
import QtQuick.Layouts
import QtQuick.Controls
import CAS

Rectangle {
    id: root
    color: AppCore.theme.bgRaised
    border.color: AppCore.theme.stroke

    function tabIndex(name) {
        const tabs = AppCore.inspectorTabs
        for (let i = 0; i < tabs.length; ++i) {
            if (tabs[i] === name) {
                return i
            }
        }
        return 0
    }

    ColumnLayout {
        anchors.fill: parent
        anchors.margins: 16
        spacing: 14

        Label {
            text: "Inspector"
            font.family: AppCore.theme.fontUI
            font.pixelSize: 18
            font.bold: true
            color: AppCore.theme.text
        }

        Rectangle {
            Layout.fillWidth: true
            visible: AppCore.helpContent.length > 0
            color: Qt.rgba(AppCore.theme.accent.r, AppCore.theme.accent.g, AppCore.theme.accent.b, 0.1)
            border.color: AppCore.theme.accentSoft
            radius: AppCore.theme.radiusM
            implicitHeight: helpText.implicitHeight + 24

            Label {
                id: helpText
                anchors.fill: parent
                anchors.margins: 12
                text: AppCore.helpContent
                textFormat: Text.RichText
                font.family: AppCore.theme.fontUI
                color: AppCore.theme.text
                wrapMode: Text.Wrap
            }
        }

        Rectangle {
            Layout.fillWidth: true
            color: AppCore.selectedCellState === "error"
                ? Qt.rgba(0.85, 0.20, 0.20, 0.10)
                : (AppCore.selectedCellState === "restored"
                    ? Qt.rgba(AppCore.theme.accent.r, AppCore.theme.accent.g, AppCore.theme.accent.b, 0.10)
                    : Qt.rgba(1, 1, 1, 0.03))
            border.color: AppCore.selectedCellState === "error" ? "#d9534f" : AppCore.theme.stroke
            radius: AppCore.theme.radiusM
            implicitHeight: stateColumn.implicitHeight + 20

            ColumnLayout {
                id: stateColumn
                anchors.fill: parent
                anchors.margins: 10
                spacing: 4

                Label {
                    text: AppCore.selectedCellInput.length > 0 ? AppCore.selectedCellInput : "Nessuna cella selezionata"
                    font.family: AppCore.theme.fontUI
                    color: AppCore.theme.text
                    wrapMode: Text.WrapAnywhere
                }

                Label {
                    text: AppCore.selectedCellState === "no-selection"
                        ? "Stato: nessuna cella"
                        : ("Stato: " + AppCore.selectedCellState + " · " + AppCore.selectedCellStatus)
                    font.family: AppCore.theme.fontMono
                    color: AppCore.theme.textMuted
                    wrapMode: Text.WrapAnywhere
                }

                Label {
                    visible: AppCore.selectedCellMeta.length > 0
                    text: AppCore.selectedCellMeta
                    font.family: AppCore.theme.fontUI
                    color: AppCore.theme.textFaint
                    wrapMode: Text.WrapAnywhere
                }
            }
        }

        TabBar {
            id: inspectorTabBar
            objectName: "inspectorTabBar"
            Layout.fillWidth: true
            currentIndex: root.tabIndex(AppCore.selectedInspectorTab)

            onCurrentIndexChanged: {
                const tabs = AppCore.inspectorTabs
                if (currentIndex >= 0 && currentIndex < tabs.length) {
                    AppCore.selectedInspectorTab = tabs[currentIndex]
                }
            }

            Repeater {
                model: AppCore.inspectorTabs
                delegate: TabButton {
                    required property var modelData
                    text: modelData
                    font.family: AppCore.theme.fontUI
                }
            }
        }

        StackLayout {
            id: inspectorStack
            objectName: "inspectorStack"
            Layout.fillWidth: true
            Layout.fillHeight: true
            currentIndex: root.tabIndex(AppCore.selectedInspectorTab)

            InspectorResultTab { objectName: "inspectorResultTab" }
            InspectorRepresentationsTab { objectName: "inspectorRepresentationsTab" }
            InspectorStepsTab { objectName: "inspectorStepsTab" }
            InspectorVariablesTab { objectName: "inspectorVariablesTab" }
            InspectorPlotTab { objectName: "inspectorPlotTab" }
        }
    }
}

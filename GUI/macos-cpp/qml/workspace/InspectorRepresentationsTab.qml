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
                visible: AppCore.selectedCellAlternatives.length === 0
                text: "Nessuna rappresentazione alternativa disponibile."
                color: AppCore.theme.textFaint
            }

            Repeater {
                model: AppCore.selectedCellAlternatives

                delegate: Rectangle {
                    required property var modelData
                    readonly property bool selectedRep: AppCore.selectedRepresentationId === modelData.id

                    Layout.fillWidth: true
                    color: selectedRep
                        ? Qt.rgba(AppCore.theme.accent.r, AppCore.theme.accent.g, AppCore.theme.accent.b, 0.12)
                        : AppCore.theme.bgCard
                    border.color: selectedRep ? AppCore.theme.accent : AppCore.theme.stroke
                    radius: AppCore.theme.radiusM
                    implicitHeight: representationLayout.implicitHeight + 20

                    TapHandler {
                        onTapped: AppCore.selectedRepresentationId = modelData.id
                    }

                    ColumnLayout {
                        id: representationLayout
                        anchors.fill: parent
                        anchors.margins: 10
                        spacing: 6

                        RowLayout {
                            Layout.fillWidth: true

                            Label {
                                text: modelData.label
                                font.family: AppCore.theme.fontUI
                                color: selectedRep ? AppCore.theme.accent : AppCore.theme.textMuted
                                font.bold: true
                            }

                            Item { Layout.fillWidth: true }

                            Label {
                                visible: selectedRep
                                text: "corrente"
                                font.family: AppCore.theme.fontMono
                                color: AppCore.theme.accent
                            }
                        }

                        Label {
                            text: modelData.value
                            color: AppCore.theme.text
                            wrapMode: Text.WrapAnywhere
                            font.family: AppCore.theme.fontMono
                        }
                    }
                }
            }
        }
    }
}

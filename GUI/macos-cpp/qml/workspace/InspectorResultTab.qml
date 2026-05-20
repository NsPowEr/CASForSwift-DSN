import QtQuick
import QtQuick.Layouts
import QtQuick.Controls
import CAS

Item {
    ScrollView {
        anchors.fill: parent

        ColumnLayout {
            width: parent.width
            spacing: 12

            Label {
                visible: !AppCore.selectedCellHasOutput
                text: AppCore.selectedCellState === "error"
                    ? "La cella ha prodotto un errore."
                    : "Nessun output disponibile per la cella selezionata."
                color: AppCore.theme.textFaint
                wrapMode: Text.Wrap
            }

            Rectangle {
                visible: AppCore.selectedCellHasOutput
                Layout.fillWidth: true
                color: AppCore.theme.bgCard
                border.color: AppCore.theme.stroke
                radius: AppCore.theme.radiusM
                implicitHeight: resultView.implicitHeight + 24

                MathView {
                    id: resultView
                    anchors.fill: parent
                    anchors.margins: 12
                    latex: AppCore.selectedCellOutput
                    fontSize: 18
                    color: AppCore.theme.text
                    block: true
                }
            }
        }
    }
}

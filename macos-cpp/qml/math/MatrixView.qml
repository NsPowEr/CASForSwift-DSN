import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import CAS

Rectangle {
    id: root
    color: AppCore.theme.bgCard
    radius: AppCore.theme.radiusM
    border.color: AppCore.theme.stroke
    clip: true

    property var model: null // MatrixTableModel*

    ColumnLayout {
        anchors.fill: parent
        anchors.margins: 12
        spacing: 12

        RowLayout {
            Label { text: "Dimensioni:"; color: AppCore.theme.textMuted }
            SpinBox { id: rows; value: 3; onValueChanged: root.model.resize(rows.value, cols.value) }
            Text { text: "×"; color: AppCore.theme.textFaint }
            SpinBox { id: cols; value: 3; onValueChanged: root.model.resize(rows.value, cols.value) }
            Item { Layout.fillWidth: true }
        }

        GridView {
            id: grid
            Layout.fillWidth: true
            Layout.fillHeight: true
            cellWidth: 60
            cellHeight: 40
            model: root.model
            
            // Note: QAbstractTableModel in GridView is tricky. 
            // In a real impl we'd use TableView (QtQuick.TableView)
            // But for scaffold a Repeater or simple TableView is better.
        }
        
        TableView {
            id: tableView
            Layout.fillWidth: true
            Layout.fillHeight: true
            model: root.model
            delegate: TextField {
                text: display
                onEditingFinished: display = text
                background: Rectangle {
                    color: AppCore.theme.bg
                    border.color: AppCore.theme.stroke
                }
            }
        }

        Button {
            Layout.fillWidth: true
            text: "Copia LaTeX"
            onClicked: {
                console.log(root.model.toLatex())
            }
        }
    }
}

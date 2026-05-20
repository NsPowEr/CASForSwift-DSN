import QtQuick
import QtQuick.Layouts
import QtQuick.Controls
import CAS

Item {
    id: root
    
    Rectangle {
        anchors.fill: parent
        color: "transparent"
        
        RowLayout {
            anchors.fill: parent
            anchors.leftMargin: 80 // Spazio per i traffic lights di macOS
            
            Label {
                text: "REAL CAS Calculator"
                font.bold: true
                font.family: AppCore.theme.fontUI
                color: AppCore.theme.text
            }
            
            Item { Layout.fillWidth: true }
            
            Button {
                text: "⌘K"
                flat: true
                font.family: AppCore.theme.fontUI
                onClicked: AppCore.palette.visible = true
            }
        }
    }
}

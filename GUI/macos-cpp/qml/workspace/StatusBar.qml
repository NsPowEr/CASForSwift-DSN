import QtQuick
import QtQuick.Layouts
import QtQuick.Controls
import CAS

Rectangle {
    id: root
    height: 28
    color: AppCore.theme.bgRaised
    border.color: AppCore.theme.stroke
    
    RowLayout {
        anchors.fill: parent
        anchors.leftMargin: 12
        anchors.rightMargin: 12
        
        Label {
            text: "Kernel: " + AppCore.kernelMode
            font.family: AppCore.theme.fontUI
            font.pixelSize: 11
            color: AppCore.theme.textFaint
        }

        Label {
            text: AppCore.storageStatus
            font.family: AppCore.theme.fontUI
            font.pixelSize: 11
            color: {
                switch (AppCore.storageSeverity) {
                    case 2: return "#e67e22" // Warning
                    case 3: return "#e74c3c" // Error
                    default: return AppCore.theme.textFaint
                }
            }
        }

        RowLayout {
            visible: AppCore.isKernelBusy
            spacing: 8
            
            BusyIndicator {
                implicitWidth: 16
                implicitHeight: 16
                running: AppCore.isKernelBusy
            }

            Button {
                text: "Stop"
                font.family: AppCore.theme.fontUI
                font.pixelSize: 10
                implicitHeight: 20
                flat: true
                onClicked: AppCore.interruptKernel()
            }
        }
        
        Item { Layout.fillWidth: true }
        
        Label {
            text: "v" + AppCore.kernelVersion
            font.family: AppCore.theme.fontUI
            font.pixelSize: 11
            color: AppCore.theme.textFaint
        }
    }
}

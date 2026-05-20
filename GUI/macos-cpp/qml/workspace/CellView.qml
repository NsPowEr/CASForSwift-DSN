import QtQuick
import QtQuick.Layouts
import QtQuick.Controls
import CAS

FocusScope {
    id: root
    required property var cell
    objectName: "cellView_" + cell.index

    readonly property bool selected: cell.active
    readonly property bool editing: cell.editing
    readonly property bool focusWithinCell: cell.focused
    readonly property bool hasVisualOutput: cell.hasOutput

    implicitHeight: frame.implicitHeight

    onActiveFocusChanged: {
        if (activeFocus && !input.activeFocus) {
            AppCore.focusCell(cell.index - 1)
        }
    }

    Rectangle {
        id: frame
        anchors.fill: parent
        color: cell.hasError
            ? Qt.rgba(0.70, 0.12, 0.12, 0.08)
            : (selected ? Qt.rgba(AppCore.theme.accent.r, AppCore.theme.accent.g, AppCore.theme.accent.b, 0.06)
                        : "transparent")
        border.color: cell.hasError
            ? "#d9534f"
            : (editing ? AppCore.theme.accent
                       : (focusWithinCell
                           ? Qt.rgba(AppCore.theme.accent.r, AppCore.theme.accent.g, AppCore.theme.accent.b, 0.60)
                           : (selected ? Qt.rgba(AppCore.theme.accent.r, AppCore.theme.accent.g, AppCore.theme.accent.b, 0.34)
                                       : AppCore.theme.stroke)))
        border.width: editing ? 2 : 1
        radius: AppCore.theme.radiusL
        implicitHeight: layout.implicitHeight + 24
    }

    TapHandler {
        acceptedButtons: Qt.LeftButton
        onTapped: {
            AppCore.notebook.selectedIndex = cell.index - 1
            AppCore.focusCell(cell.index - 1)
            root.forceActiveFocus()
        }
        onDoubleTapped: {
            AppCore.notebook.selectedIndex = cell.index - 1
            AppCore.beginEditingCell(cell.index - 1)
            input.forceActiveFocus()
        }
    }

    ColumnLayout {
        id: layout
        anchors.fill: frame
        anchors.margins: 14
        spacing: 8

        RowLayout {
            spacing: 7

            Rectangle {
                radius: 5
                color: selected ? AppCore.theme.accent : AppCore.theme.accentSoft
                implicitWidth: badge.implicitWidth + 16
                implicitHeight: 18

                Text {
                    id: badge
                    anchors.centerIn: parent
                    text: "In[" + cell.index + "]"
                    color: selected ? "#fff" : AppCore.theme.accent
                    font { family: AppCore.theme.fontMono; pointSize: 9; bold: true }
                }
            }

            Rectangle {
                radius: 5
                color: editing
                    ? Qt.rgba(AppCore.theme.accent.r, AppCore.theme.accent.g, AppCore.theme.accent.b, 0.18)
                    : (cell.hasError
                        ? Qt.rgba(0.85, 0.25, 0.25, 0.18)
                        : (hasVisualOutput
                            ? Qt.rgba(0.18, 0.50, 0.20, 0.14)
                            : Qt.rgba(1, 1, 1, 0.04)))
                border.color: editing ? AppCore.theme.accent : "transparent"
                implicitWidth: stateText.implicitWidth + 16
                implicitHeight: 18

                Text {
                    id: stateText
                    anchors.centerIn: parent
                    text: editing ? "editing" : (focusWithinCell ? "focused" : (selected ? "selected" : cell.statusText))
                    color: editing ? AppCore.theme.accent : AppCore.theme.textMuted
                    font { family: AppCore.theme.fontMono; pointSize: 9 }
                }
            }

            Item { Layout.fillWidth: true }

            RowLayout {
                spacing: 4

                ToolButton { text: "+↑"; enabled: AppCore.isCommandEnabled("insert_cell_above"); onClicked: AppCore.insertCellAbove() }
                ToolButton { text: "+↓"; enabled: AppCore.isCommandEnabled("insert_cell_below"); onClicked: AppCore.insertCellBelow() }
                ToolButton { text: "⇱"; enabled: AppCore.isCommandEnabled("merge_with_previous"); onClicked: AppCore.mergeWithPreviousCell() }
                ToolButton { text: "⇲"; enabled: AppCore.isCommandEnabled("merge_with_next"); onClicked: AppCore.mergeWithNextCell() }
                ToolButton { text: "⫶"; enabled: AppCore.isCommandEnabled("split_cell"); onClicked: AppCore.splitCurrentCell(input.cursorPosition) }
                ToolButton { text: "↑"; enabled: AppCore.isCommandEnabled("move_cell_up"); onClicked: AppCore.moveCurrentCellUp() }
                ToolButton { text: "↓"; enabled: AppCore.isCommandEnabled("move_cell_down"); onClicked: AppCore.moveCurrentCellDown() }
                ToolButton { text: "⧉"; enabled: AppCore.isCommandEnabled("duplicate_cell"); onClicked: AppCore.duplicateCurrentCell() }
                ToolButton { text: "✕"; enabled: AppCore.isCommandEnabled("delete_cell"); onClicked: AppCore.deleteCurrentCell() }
            }

            Text {
                text: "⇧↵ run · ⌘↵ next · ⌃↑↓ select"
                color: AppCore.theme.textFaint
                font { family: AppCore.theme.fontMono; pointSize: 9 }
            }
        }

        Rectangle {
            Layout.fillWidth: true
            color: cell.hasError ? Qt.rgba(1, 0, 0, 0.04) : AppCore.theme.bgRaised
            border.color: cell.hasError
                ? "#ff4444"
                : (editing
                    ? AppCore.theme.accent
                    : (selected
                        ? Qt.rgba(AppCore.theme.accent.r, AppCore.theme.accent.g, AppCore.theme.accent.b, 0.4)
                        : AppCore.theme.stroke))
            border.width: editing ? 2 : 1
            radius: AppCore.theme.radiusM
            implicitHeight: input.implicitHeight + 20

            MathInput {
                id: input
                objectName: "cellInput"
                anchors.fill: parent
                anchors.margins: 10
                latex: cell.inputLatex
                focus: cell.editing

                onActiveFocusChanged: {
                    if (activeFocus) {
                        AppCore.notebook.selectedIndex = cell.index - 1
                        AppCore.beginEditingCell(cell.index - 1)
                    } else if (cell.editing) {
                        AppCore.endEditingCell(cell.index - 1)
                    }
                }
                onTextChanged: cell.inputLatex = text
                onSubmitted: {
                    AppCore.notebook.selectedIndex = cell.index - 1
                    AppCore.executeCurrentCell()
                }
                onRunAndAdvanceRequested: {
                    AppCore.notebook.selectedIndex = cell.index - 1
                    AppCore.runAndAdvance()
                }
            }
        }

        // ── Live Preview (Ghost)
        Rectangle {
            id: livePreview
            visible: editing && input.text.length > 0
            Layout.fillWidth: true
            Layout.topMargin: -4
            color: "transparent"
            border.color: Qt.rgba(AppCore.theme.accent.r, AppCore.theme.accent.g, AppCore.theme.accent.b, 0.25)
            border.width: 1
            radius: AppCore.theme.radiusM
            implicitHeight: previewView.implicitHeight + 16
            
            // Dashed border effect
            Rectangle { anchors.fill: parent; radius: parent.radius; color: "transparent"; border.color: parent.border.color; border.width: 1; visible: false }

            Label {
                anchors.left: parent.left
                anchors.top: parent.top
                anchors.margins: 4
                text: "PREVIEW"
                font.family: AppCore.theme.fontMono
                font.pixelSize: 8
                font.bold: true
                color: AppCore.theme.accent
                opacity: 0.6
            }

            MathView {
                id: previewView
                anchors.centerIn: parent
                width: parent.width - 24
                latex: input.text
                fontSize: 18
                block: true
                color: AppCore.theme.textMuted
            }
        }

        Loader {
            visible: cell.hasOutput
            Layout.fillWidth: true
            sourceComponent: outputComp
        }

        Component {
            id: outputComp

            ColumnLayout {
                spacing: 4

                RowLayout {
                    Rectangle {
                        radius: 5
                        color: AppCore.theme.bgCard
                        implicitWidth: outBadge.implicitWidth + 16
                        implicitHeight: 18

                        Text {
                            id: outBadge
                            anchors.centerIn: parent
                            text: "Out[" + cell.index + "]"
                            color: AppCore.theme.textMuted
                            font { family: AppCore.theme.fontMono; pointSize: 9; bold: true }
                        }
                    }

                    Text {
                        text: cell.outputMeta
                        color: AppCore.theme.textFaint
                        font { family: AppCore.theme.fontMono; pointSize: 9 }
                    }

                    Rectangle {
                        visible: cell.restoredFromWorkspace
                        radius: 5
                        color: Qt.rgba(AppCore.theme.accent.r, AppCore.theme.accent.g, AppCore.theme.accent.b, 0.14)
                        implicitWidth: restoredLabel.implicitWidth + 14
                        implicitHeight: 18

                        Text {
                            id: restoredLabel
                            anchors.centerIn: parent
                            text: "restored"
                            color: AppCore.theme.accent
                            font { family: AppCore.theme.fontMono; pointSize: 8; bold: true }
                        }
                    }
                }

                Rectangle {
                    Layout.fillWidth: true
                    color: AppCore.theme.bg
                    border.color: AppCore.theme.stroke
                    radius: AppCore.theme.radiusM
                    implicitHeight: outView.implicitHeight + 24

                    MathView {
                        id: outView
                        anchors.fill: parent
                        anchors.margins: 12
                        latex: cell.outputLatex
                        fontSize: 22
                        block: true
                    }
                }

                RowLayout {
                    visible: cell.alternatives.length > 0
                    spacing: 6

                    Repeater {
                        model: cell.alternatives

                        delegate: Rectangle {
                            required property var modelData
                            readonly property bool highlighted: AppCore.selectedRepresentationId === modelData.id

                            radius: 7
                            color: highlighted
                                ? Qt.rgba(AppCore.theme.accent.r, AppCore.theme.accent.g, AppCore.theme.accent.b, 0.14)
                                : AppCore.theme.bgCard
                            border.color: highlighted ? AppCore.theme.accent : AppCore.theme.stroke
                            implicitWidth: altText.implicitWidth + 20
                            implicitHeight: 24

                            Text {
                                id: altText
                                anchors.centerIn: parent
                                text: modelData.label
                                color: highlighted ? AppCore.theme.accent : AppCore.theme.textMuted
                                font { family: AppCore.theme.fontMono; pointSize: 9 }
                            }

                            TapHandler {
                                onTapped: AppCore.selectedRepresentationId = modelData.id
                            }
                        }
                    }
                }
            }
        }
    }
}

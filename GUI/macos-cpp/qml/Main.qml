// Main.qml — root window con TitleBar nativa + body 3-colonne.
// Riproduce la struttura di CASMacWorkspace dal prototipo HTML.

import QtQuick
import QtQuick.Controls
import QtQuick.Controls.macOS
import QtQuick.Layouts
import QtQuick.Window
import Qt.labs.platform as Platform
import CAS

ApplicationWindow {
    id: win
    objectName: "mainWindow"
    width: 1280; height: 820
    minimumWidth: 1024; minimumHeight: 640
    title: "CAS Calculator — " + AppCore.notebook.sessionTitle
    visible: true
    font.family: AppCore.theme.fontUI

    Platform.MenuBar {
        Platform.Menu {
            title: "File"
            Platform.MenuItem { text: "Nuova Sessione"; shortcut: "Meta+N"; onTriggered: AppCore.newSession() }
            Platform.MenuItem { text: "Esporta come HTML..."; onTriggered: exportDialog.open() }
            Platform.MenuSeparator {}
            Platform.MenuItem { text: "Reset Workspace"; onTriggered: AppCore.clearWorkspace() }
            Platform.MenuSeparator {}
            Platform.MenuItem { text: "Esci"; role: Platform.MenuItem.QuitRole; onTriggered: Qt.quit() }
        }
        Platform.Menu {
            title: "Modifica"
            Platform.MenuItem { text: AppCore.undoActionLabel; enabled: AppCore.canUndo; shortcut: "Meta+Z"; onTriggered: AppCore.undo() }
            Platform.MenuItem { text: AppCore.redoActionLabel; enabled: AppCore.canRedo; shortcut: "Meta+Shift+Z"; onTriggered: AppCore.redo() }
            Platform.MenuSeparator {}
            Platform.MenuItem { text: "Inserisci Cella Sopra"; shortcut: "Meta+Alt+Up"; onTriggered: AppCore.insertCellAbove() }
            Platform.MenuItem { text: "Inserisci Cella Sotto"; shortcut: "Meta+Alt+Down"; onTriggered: AppCore.insertCellBelow() }
            Platform.MenuItem { text: "Duplica Cella"; shortcut: "Meta+Shift+D"; onTriggered: AppCore.duplicateCurrentCell() }
            Platform.MenuItem { text: "Elimina Cella"; enabled: AppCore.canDeleteCell; shortcut: "Meta+Delete"; onTriggered: AppCore.deleteCurrentCell() }
            Platform.MenuSeparator {}
            Platform.MenuItem { text: "Unisci con Precedente"; shortcut: "Meta+Alt+Left"; onTriggered: AppCore.mergeWithPreviousCell() }
            Platform.MenuItem { text: "Unisci con Successiva"; shortcut: "Meta+Alt+Right"; onTriggered: AppCore.mergeWithNextCell() }
        }
        Platform.Menu {
            title: "Sessione"
            Platform.MenuItem { text: "Duplica Sessione"; shortcut: "Meta+D"; onTriggered: AppCore.duplicateSession() }
            Platform.MenuItem { text: "Elimina Sessione"; enabled: AppCore.canDeleteSession; shortcut: "Meta+Backspace"; onTriggered: AppCore.deleteSession() }
            Platform.MenuSeparator {}
            Platform.MenuItem { text: "Prossima Sessione"; shortcut: "Ctrl+Tab"; onTriggered: AppCore.sessions.activeIndex = (AppCore.sessions.activeIndex + 1) % AppCore.sessions.titles.length }
            Platform.MenuItem { text: "Sessione Precedente"; shortcut: "Ctrl+Shift+Tab"; onTriggered: AppCore.sessions.activeIndex = (AppCore.sessions.activeIndex - 1 + AppCore.sessions.titles.length) % AppCore.sessions.titles.length }
        }
        Platform.Menu {
            title: "Kernel"
            Platform.MenuItem { text: "Esegui Cella"; shortcut: "Shift+Return"; onTriggered: AppCore.executeCurrentCell() }
            Platform.MenuItem { text: "Esegui e Avanza"; shortcut: "Meta+Return"; onTriggered: AppCore.runAndAdvance() }
            Platform.MenuSeparator {}
            Platform.MenuItem { text: "Interrompi"; enabled: AppCore.isKernelBusy; onTriggered: AppCore.interruptKernel() }
            Platform.MenuSeparator {}
            Platform.MenuItem { text: "Ricarica da Disco"; shortcut: "Meta+Shift+R"; onTriggered: AppCore.reloadWorkspaceFromDisk() }
        }
        Platform.Menu {
            title: "Finestra"
            Platform.MenuItem { text: "Command Palette"; shortcut: "Meta+K"; onTriggered: AppCore.palette.toggle() }
            Platform.MenuItem { text: AppCore.sidebarVisible ? "Nascondi Sidebar" : "Mostra Sidebar"; shortcut: "Meta+Alt+S"; onTriggered: AppCore.toggleSidebar() }
            Platform.MenuItem { text: AppCore.inspectorVisible ? "Nascondi Inspector" : "Mostra Inspector"; shortcut: "Meta+Alt+I"; onTriggered: AppCore.toggleInspector() }
            Platform.MenuItem { text: "Ripristina Layout"; shortcut: "Meta+Alt+0"; onTriggered: AppCore.resetLayoutState() }
            Platform.MenuSeparator {}
            Platform.MenuItem { text: "Plot Tab"; shortcut: "Meta+Shift+P"; onTriggered: { AppCore.selectedInspectorTab = "Plot"; AppCore.plotSelectedCell() } }
        }
    }

    Platform.FileDialog {
        id: exportDialog
        title: "Esporta Workspace come HTML"
        fileMode: Platform.FileDialog.SaveFile
        nameFilters: ["HTML files (*.html)"]
        onAccepted: AppCore.exportToHtml(file.toString().replace("file://", ""))
    }

    background: Rectangle { color: AppCore.theme.bg }
    color: AppCore.theme.bg

    // ── TitleBar unificata (overlay, hidesTitleBar)
    TitleBar {
        id: titleBar
        objectName: "titleBar"
        anchors { top: parent.top; left: parent.left; right: parent.right }
        height: 44
    }

    // ── Body: Sidebar | Notebook | Inspector
    SplitView {
        id: workspaceSplit
        objectName: "workspaceSplit"
        anchors { top: titleBar.bottom; left: parent.left; right: parent.right; bottom: statusBar.top }

        Sidebar {
            id: sidebar
            objectName: "sidebar"
            visible: AppCore.sidebarVisible
            SplitView.preferredWidth: AppCore.sidebarVisible ? AppCore.sidebarWidth : 0
            SplitView.minimumWidth: 180
            SplitView.maximumWidth: 420
            SplitView.fillHeight: true
            onWidthChanged: if (visible && width > 0) AppCore.sidebarWidth = width
        }

        Notebook {
            id: notebook
            objectName: "notebook"
            SplitView.fillWidth: true
            SplitView.fillHeight: true
        }

        Inspector {
            id: inspector
            objectName: "inspector"
            visible: AppCore.inspectorVisible
            SplitView.preferredWidth: AppCore.inspectorVisible ? AppCore.inspectorWidth : 0
            SplitView.minimumWidth: 260
            SplitView.maximumWidth: 520
            SplitView.fillHeight: true
            onWidthChanged: if (visible && width > 0) AppCore.inspectorWidth = width
        }
    }

    StatusBar {
        id: statusBar
        objectName: "statusBar"
        anchors { left: parent.left; right: parent.right; bottom: parent.bottom }
    }

    // ── Command palette (⌘K) flottante
    CommandPalette {
        id: palette
        objectName: "commandPalette"
        anchors.horizontalCenter: parent.horizontalCenter
        y: 120; width: 520
        visible: AppCore.palette.visible
    }

    Shortcut { objectName: "shortcutPalette"; sequence: "Meta+K"; onActivated: AppCore.palette.toggle() }
    Shortcut { objectName: "shortcutRunMeta"; sequence: "Meta+Return"; onActivated: AppCore.runAndAdvance() }
    Shortcut { objectName: "shortcutRunShift"; sequence: "Shift+Return"; onActivated: AppCore.executeCurrentCell() }
    Shortcut { objectName: "shortcutPlot"; sequence: "Meta+Shift+P"; onActivated: AppCore.plotSelectedCell() }
    Shortcut { objectName: "shortcutReloadWorkspace"; sequence: "Meta+Shift+R"; onActivated: AppCore.reloadWorkspaceFromDisk() }
    Shortcut { objectName: "shortcutToggleSidebar"; sequence: "Meta+Alt+S"; onActivated: AppCore.toggleSidebar() }
    Shortcut { objectName: "shortcutToggleInspector"; sequence: "Meta+Alt+I"; onActivated: AppCore.toggleInspector() }
    Shortcut { objectName: "shortcutResetLayout"; sequence: "Meta+Alt+0"; onActivated: AppCore.resetLayoutState() }
    Shortcut { objectName: "shortcutNewSession"; sequence: "Meta+N"; onActivated: AppCore.newSession() }
    Shortcut { objectName: "shortcutDuplicateSession"; sequence: "Meta+D"; onActivated: AppCore.duplicateSession() }
    Shortcut { objectName: "shortcutDuplicateCell"; sequence: "Meta+Shift+D"; onActivated: AppCore.duplicateCurrentCell() }
    Shortcut { objectName: "shortcutInsertCellAbove"; sequence: "Meta+Alt+Up"; onActivated: AppCore.insertCellAbove() }
    Shortcut { objectName: "shortcutInsertCellBelow"; sequence: "Meta+Alt+Down"; onActivated: AppCore.insertCellBelow() }
    Shortcut { objectName: "shortcutDeleteSession"; sequence: "Meta+Backspace"; onActivated: AppCore.deleteSession() }
    Shortcut { objectName: "shortcutDeleteCell"; sequence: "Meta+Delete"; onActivated: AppCore.deleteCurrentCell() }
    Shortcut { objectName: "shortcutSplitCell"; sequence: "Meta+Shift+\\"; onActivated: AppCore.splitCurrentCell() }
    Shortcut { objectName: "shortcutMergePrev"; sequence: "Meta+Alt+Left"; onActivated: AppCore.mergeWithPreviousCell() }
    Shortcut { objectName: "shortcutMergeNext"; sequence: "Meta+Alt+Right"; onActivated: AppCore.mergeWithNextCell() }
    Shortcut { objectName: "shortcutSelectPrev"; sequence: "Ctrl+Up"; onActivated: AppCore.selectPreviousCell() }
    Shortcut { objectName: "shortcutSelectNext"; sequence: "Ctrl+Down"; onActivated: AppCore.selectNextCell() }
    Shortcut { objectName: "shortcutUndo"; sequence: "Meta+Z"; onActivated: AppCore.undo() }
    Shortcut { objectName: "shortcutRedo"; sequence: "Meta+Shift+Z"; onActivated: AppCore.redo() }
    Shortcut { objectName: "shortcutMoveCellUp"; sequence: "Alt+Up"; onActivated: AppCore.moveCurrentCellUp() }
    Shortcut { objectName: "shortcutMoveCellDown"; sequence: "Alt+Down"; onActivated: AppCore.moveCurrentCellDown() }
}

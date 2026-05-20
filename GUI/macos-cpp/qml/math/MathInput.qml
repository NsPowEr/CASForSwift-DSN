import QtQuick
import QtQuick.Controls
import CAS

FocusScope {
    id: root
    property alias text: input.text
    property alias latex: input.text
    property alias cursorPosition: input.cursorPosition
    property alias selectionStart: input.selectionStart
    property alias selectionEnd: input.selectionEnd
    property alias selectedText: input.selectedText
    signal submitted()
    signal runAndAdvanceRequested()

    function tokenStartForPosition(pos) {
        let txt = input.text
        let start = pos - 1
        while (start >= 0 && /[a-zA-Z0-9_\\]/.test(txt[start])) {
            start--
        }
        return start + 1
    }

    function moveCursorIntoTemplate(start, insertion, hadSelection) {
        let placeholders = insertion.match(/\{\}/g)
        if (placeholders && placeholders.length > 0) {
            let firstIdx = insertion.indexOf("{}")
            if (hadSelection && placeholders.length > 1) {
                let secondIdx = insertion.indexOf("{}", firstIdx + 2)
                input.cursorPosition = start + secondIdx + 1
                return
            }
            input.cursorPosition = start + firstIdx + 1
            return
        }
        if (insertion.endsWith("()")) {
            input.cursorPosition = start + insertion.length - 1
            return
        }
        input.cursorPosition = start + insertion.length
    }

    function buildInsertion(insertText, type, selected) {
        let insertion = insertText || ""
        if (selected.length === 0) {
            return insertion
        }

        if (type === "function" && insertion.endsWith("()")) {
            return insertion.slice(0, insertion.length - 1) + selected + ")"
        }
        if (type === "template") {
            if (insertion.indexOf("\\frac{}{}") === 0) {
                return "\\frac{" + selected + "}{}"
            }
            if (insertion.indexOf("{}") !== -1) {
                return insertion.replace("{}", "{" + selected + "}")
            }
        }
        return insertion
    }

    TextArea {
        id: input
        anchors.fill: parent
        font.family: AppCore.theme.fontMono
        font.pointSize: 14
        color: AppCore.theme.text
        background: null
        focus: true

        MathHighlighter {
            target: input
        }

        Keys.onPressed: (event) => {
            if (autocomplete.visible) {
                if (event.key === Qt.Key_Up) {
                    autocomplete.moveUp()
                    event.accepted = true
                    return
                } else if (event.key === Qt.Key_Down) {
                    autocomplete.moveDown()
                    event.accepted = true
                    return
                } else if (event.key === Qt.Key_Tab || event.key === Qt.Key_Return) {
                    if (autocomplete.acceptCurrent()) {
                        event.accepted = true
                        return
                    }
                } else if (event.key === Qt.Key_Escape) {
                    autocomplete.visible = false
                    event.accepted = true
                    return
                }
            }

            if (event.key === Qt.Key_Return && (event.modifiers & Qt.ShiftModifier)) {
                submitted()
                event.accepted = true
            } else if (event.key === Qt.Key_Return && (event.modifiers & Qt.MetaModifier)) {
                runAndAdvanceRequested()
                event.accepted = true
            }
        }

        onCursorPositionChanged: {
            updateAutocomplete()
        }

        onTextChanged: {
            updateAutocomplete()
        }

        function updateAutocomplete() {
            let pos = cursorPosition
            let start = root.tokenStartForPosition(pos)
            let prefix = text.substring(start, pos)
            if (prefix.length >= 1 && (/^[a-zA-Z]/.test(prefix) || /^\\[a-zA-Z]/.test(prefix))) {
                autocomplete.prefix = prefix
                let cursorRect = input.cursorRectangle
                autocomplete.x = cursorRect.x
                autocomplete.y = cursorRect.y + cursorRect.height + 4
            } else {
                autocomplete.visible = false
            }
        }
    }

    AutocompletePopup {
        id: autocomplete
        onSelected: (selectedText, insertText, type) => {
            let pos = input.cursorPosition
            let start = root.tokenStartForPosition(pos)
            let rangeStart = Math.min(input.selectionStart, input.selectionEnd)
            let rangeEnd = Math.max(input.selectionStart, input.selectionEnd)
            let activeSelection = rangeEnd > rangeStart ? input.selectedText : ""
            let replaceStart = activeSelection.length > 0 ? rangeStart : start
            let replaceEnd = activeSelection.length > 0 ? rangeEnd : pos
            let insertion = root.buildInsertion(insertText || selectedText, type, activeSelection)
            input.remove(replaceStart, replaceEnd)
            input.insert(replaceStart, insertion)

            root.moveCursorIntoTemplate(replaceStart, insertion, activeSelection.length > 0)
            input.focus = true
        }
    }
}

import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

ApplicationWindow {
    id: root
    width: 1180
    height: 760
    visible: true
    title: "REAL CAS GUI Lab"

    color: "#101114"

    RowLayout {
        anchors.fill: parent
        spacing: 0

        Rectangle {
            Layout.preferredWidth: 260
            Layout.fillHeight: true
            color: "#181A20"

            ColumnLayout {
                anchors.fill: parent
                anchors.margins: 18
                spacing: 14

                Label {
                    text: "CAS Lab"
                    color: "#F5F7FA"
                    font.pixelSize: 26
                    font.bold: true
                }

                Label {
                    text: "Manual tester staccabile"
                    color: "#9DA5B4"
                    wrapMode: Text.WordWrap
                }

                Rectangle {
                    Layout.fillWidth: true
                    height: 1
                    color: "#2A2D36"
                }

                Label {
                    text: "Operazioni v1"
                    color: "#F5F7FA"
                    font.bold: true
                }

                Label {
                    text: "Simplify\nLaTeX\nASCII 2D\nPlot 2D"
                    color: "#C9D1DF"
                    lineHeight: 1.25
                }

                Item { Layout.fillHeight: true }

                Label {
                    text: "Core untouched: GUI -> adapter -> cas_core"
                    color: "#7E8796"
                    wrapMode: Text.WordWrap
                }
            }
        }

        ColumnLayout {
            Layout.fillWidth: true
            Layout.fillHeight: true
            spacing: 12

            Rectangle {
                Layout.fillWidth: true
                Layout.preferredHeight: 112
                color: "#12141A"

                RowLayout {
                    anchors.fill: parent
                    anchors.margins: 18
                    spacing: 12

                    TextArea {
                        Layout.fillWidth: true
                        Layout.fillHeight: true
                        text: backend.input
                        color: "#F5F7FA"
                        selectedTextColor: "#111318"
                        selectionColor: "#A6D189"
                        placeholderText: "Inserisci una formula: integrate(sin(x)^2, x), diff(x*sin(x), x), sin(x)"
                        onTextChanged: backend.input = text
                        background: Rectangle {
                            color: "#1C1F28"
                            radius: 14
                            border.color: "#303543"
                        }
                    }

                    Button {
                        text: "Calcola"
                        onClicked: backend.compute()
                    }

                    Button {
                        text: "Plot 2D"
                        onClicked: {
                            backend.compute()
                            backend.sample2d("x", -6.28, 6.28)
                        }
                    }
                }
            }

            SplitView {
                Layout.fillWidth: true
                Layout.fillHeight: true
                orientation: Qt.Horizontal

                ScrollView {
                    SplitView.preferredWidth: 520

                    ColumnLayout {
                        width: parent.width
                        spacing: 12

                        OutputCard {
                            title: "Status"
                            body: backend.status
                        }

                        OutputCard {
                            title: "Text"
                            body: backend.textResult
                        }

                        OutputCard {
                            title: "LaTeX"
                            body: backend.latexResult
                        }

                        OutputCard {
                            title: "Numeric"
                            body: backend.numericResult
                        }

                        OutputCard {
                            title: "ASCII 2D"
                            body: backend.asciiResult
                            mono: true
                        }
                    }
                }

                Rectangle {
                    SplitView.fillWidth: true
                    color: "#0D0F14"

                    ColumnLayout {
                        anchors.fill: parent
                        anchors.margins: 18

                        Label {
                            text: "Plot 2D"
                            color: "#F5F7FA"
                            font.pixelSize: 20
                            font.bold: true
                        }

                        Canvas {
                            id: plotCanvas
                            Layout.fillWidth: true
                            Layout.fillHeight: true
                            onPaint: {
                                const ctx = getContext("2d")
                                ctx.reset()
                                ctx.fillStyle = "#11141B"
                                ctx.fillRect(0, 0, width, height)

                                const points = JSON.parse(backend.plotJson)
                                if (points.length < 2) {
                                    return
                                }

                                let minY = points[0].y
                                let maxY = points[0].y
                                for (const p of points) {
                                    minY = Math.min(minY, p.y)
                                    maxY = Math.max(maxY, p.y)
                                }
                                const spanY = Math.max(maxY - minY, 1e-9)
                                const minX = points[0].x
                                const maxX = points[points.length - 1].x
                                const spanX = Math.max(maxX - minX, 1e-9)

                                ctx.strokeStyle = "#2A3040"
                                ctx.lineWidth = 1
                                ctx.beginPath()
                                ctx.moveTo(0, height / 2)
                                ctx.lineTo(width, height / 2)
                                ctx.moveTo(width / 2, 0)
                                ctx.lineTo(width / 2, height)
                                ctx.stroke()

                                ctx.strokeStyle = "#A6D189"
                                ctx.lineWidth = 2
                                ctx.beginPath()
                                for (let i = 0; i < points.length; ++i) {
                                    const x = ((points[i].x - minX) / spanX) * width
                                    const y = height - ((points[i].y - minY) / spanY) * height
                                    if (i === 0) {
                                        ctx.moveTo(x, y)
                                    } else {
                                        ctx.lineTo(x, y)
                                    }
                                }
                                ctx.stroke()
                            }

                            Connections {
                                target: backend
                                function onPlotChanged() { plotCanvas.requestPaint() }
                            }
                        }
                    }
                }
            }
        }
    }

    component OutputCard: Rectangle {
        required property string title
        required property string body
        property bool mono: false

        width: parent.width
        implicitHeight: content.implicitHeight + 28
        color: "#191C24"
        radius: 16
        border.color: "#2A2E3A"

        ColumnLayout {
            id: content
            anchors.left: parent.left
            anchors.right: parent.right
            anchors.margins: 14
            anchors.verticalCenter: parent.verticalCenter
            spacing: 8

            Label {
                text: title
                color: "#9DA5B4"
                font.bold: true
            }

            Text {
                Layout.fillWidth: true
                text: body.length > 0 ? body : "..."
                color: "#F5F7FA"
                wrapMode: Text.Wrap
                font.family: mono ? "Menlo" : "Helvetica Neue"
            }
        }
    }
}

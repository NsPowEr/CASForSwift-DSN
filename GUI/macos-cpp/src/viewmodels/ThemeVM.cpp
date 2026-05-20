// ThemeVM.cpp — porting della funzione palette() dal prototipo HTML.

#include "ThemeVM.h"
#include <QSettings>

ThemeVM::ThemeVM(QObject* parent) : QObject(parent) { load(); }

void ThemeVM::setDark(bool v) {
    if (m_dark == v) return;
    m_dark = v; emit changed(); persist();
}

void ThemeVM::setAccent(QColor c) {
    if (m_accent == c) return;
    m_accent = c; emit changed(); persist();
}

// Esattamente le derivazioni del prototipo:
//   bg:        dark ? #0B0B0E : #FFFFFF
//   bgRaised:  dark ? #17171C : #F7F7F9
//   bgCard:    dark ? #1F1F26 : #FFFFFF
//   stroke:    dark ? rgba(255,255,255,0.07) : rgba(0,0,0,0.07)
//   text:      dark ? #ECECEF : #1A1A1F
//   textMuted: dark ? #9D9DA8 : #6B6B75
//   textFaint: dark ? #5C5C66 : #9C9CA6
//   accentSoft: accent + alpha 0x22 (light) or 0x33 (dark)

QColor ThemeVM::bg() const        { return m_dark ? QColor("#0B0B0E") : QColor("#FFFFFF"); }
QColor ThemeVM::bgRaised() const  { return m_dark ? QColor("#17171C") : QColor("#F7F7F9"); }
QColor ThemeVM::bgCard() const    { return m_dark ? QColor("#1F1F26") : QColor("#FFFFFF"); }
QColor ThemeVM::stroke() const    { return m_dark ? QColor(255,255,255, 18) : QColor(0,0,0, 18); }
QColor ThemeVM::text() const      { return m_dark ? QColor("#ECECEF") : QColor("#1A1A1F"); }
QColor ThemeVM::textMuted() const { return m_dark ? QColor("#9D9DA8") : QColor("#6B6B75"); }
QColor ThemeVM::textFaint() const { return m_dark ? QColor("#5C5C66") : QColor("#9C9CA6"); }

QColor ThemeVM::accentSoft() const {
    QColor c = m_accent;
    c.setAlpha(m_dark ? 51 : 34); // 0x33 / 0x22
    return c;
}

void ThemeVM::persist() {
    QSettings s; s.setValue("theme/dark", m_dark);
    s.setValue("theme/accent", m_accent.name(QColor::HexArgb));
}

void ThemeVM::load() {
    QSettings s;
    if (s.contains("theme/dark"))   m_dark = s.value("theme/dark").toBool();
    if (s.contains("theme/accent")) m_accent = QColor(s.value("theme/accent").toString());
}

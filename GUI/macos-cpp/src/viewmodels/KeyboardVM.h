#pragma once
#include <QObject>

namespace cas::gui {

class KeyboardVM : public QObject {
    Q_OBJECT
public:
    explicit KeyboardVM(QObject* parent = nullptr) : QObject(parent) {}
    // Placeholder for math keyboard logic
};

} // namespace cas::gui

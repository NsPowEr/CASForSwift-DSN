#pragma once

#include <functional>

class QWindow;

namespace cas::platform {

void applyMacOSChrome(QWindow* window);
void observeAppearance(std::function<void(bool dark)> callback);

} // namespace cas::platform

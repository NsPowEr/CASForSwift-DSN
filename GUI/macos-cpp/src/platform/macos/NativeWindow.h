// NativeWindow.h — pure C++ interface to macOS-specific window chrome.
// No Objective-C types or AppKit headers here: consumers stay Qt-only.
//
// Two namespaces coexist for backward compatibility:
//   - cas::platform   : original API (applyMacOSChrome, observeAppearance)
//   - cas::gui::native: granular API requested by the task spec
//
// This header is safe to include from any .cpp or .h file.
// The implementation lives in NativeWindow.mm (Objective-C++).

#pragma once

#include <functional>

// Forward-declare QWindow so consumers do not pull in Qt headers.
class QWindow;

// ─── Legacy API (originally in cas::platform) ────────────────────────────────
// Kept for backward compatibility with existing callers (e.g. AppCore).

namespace cas::platform {

/// Apply the full macOS chrome in one shot:
///   • titlebarAppearsTransparent + NSWindowTitleHidden
///   • NSWindowStyleMaskFullSizeContentView
///   • NSVisualEffectView behind the entire window
///   • Traffic-light buttons repositioned to 12 px / 16 px
void applyMacOSChrome(QWindow* window);

/// Register a callback that is invoked every time the system appearance
/// switches between dark and light mode.  The bool argument is true when
/// the new appearance is Dark.
void observeAppearance(std::function<void(bool dark)> callback);

} // namespace cas::platform

// ─── Granular API (cas::gui::native) ─────────────────────────────────────────
// Each function is intentionally narrow so callers can compose them.

namespace cas::gui::native {

/// Attach an NSVisualEffectView behind the Qt content view of @p window.
/// Material: NSVisualEffectMaterialHUDWindow (frosted, high-contrast).
/// Blending mode: NSVisualEffectBlendingModeBehindWindow.
/// The effect view is auto-resized together with the window.
///
/// Must be called after the window is shown (i.e. the NSWindow handle
/// is already available from the Qt platform).
void apply_vibrancy(QWindow* window);

/// Reposition the three traffic-light buttons (close / minimise / zoom)
/// so that their top-left corner is at (@p x, @p y) in the NSWindow's
/// coordinate system (origin at bottom-left, points not pixels).
///
/// Typical values matching a 44-pt unified title bar: x=12, y=16.
void set_traffic_lights_pos(QWindow* window, double x, double y);

/// Configure the NSWindow for a unified, borderless title-bar look:
///   • titlebarAppearsTransparent = YES
///   • titleVisibility = NSWindowTitleHidden
///   • styleMask |= NSWindowStyleMaskFullSizeContentView
///   • movableByWindowBackground = YES
///
/// This is the minimal subset needed for a "frameless" macOS window
/// where Qt content extends behind the traffic lights.
void set_title_bar_style(QWindow* window);

} // namespace cas::gui::native

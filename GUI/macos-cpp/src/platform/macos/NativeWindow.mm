// NativeWindow.mm — Objective-C++ implementation of macOS native window chrome.
//
// Responsibilities:
//   1. Expose native NSWindow from a QWindow handle via Qt's platform interface.
//   2. Apply vibrancy (NSVisualEffectView) behind the Qt content view.
//   3. Reposition traffic-light buttons to a custom (x, y) coordinate.
//   4. Set the unified, full-size-content-view title-bar style.
//   5. Observe system appearance changes (dark/light) and fire a C++ callback.
//
// Objective-C ARC is enabled via the -fobjc-arc compile flag set in CMake.
// Do NOT mix manual retain/release calls here.
//
// Compile guard: this file is only compiled on Apple platforms (see CMakeLists).

#include "NativeWindow.h"

// ── Objective-C / Cocoa ───────────────────────────────────────────────────────
#import <AppKit/AppKit.h>
#import <Cocoa/Cocoa.h>

// ── Qt ────────────────────────────────────────────────────────────────────────
#include <QWindow>
#include <QGuiApplication>
#include <QtGui/qpa/qplatformnativeinterface.h>

// ─────────────────────────────────────────────────────────────────────────────
// Internal helper: retrieve the NSWindow backing a QWindow.
// Returns nil if the window is not yet shown or the platform is not Cocoa.
// ─────────────────────────────────────────────────────────────────────────────
namespace {

[[nodiscard]] NSWindow* ns_window_for(QWindow* qw) {
    if (!qw) return nil;
    auto* pni = QGuiApplication::platformNativeInterface();
    if (!pni) return nil;
    return reinterpret_cast<NSWindow*>(
        pni->nativeResourceForWindow("nswindow", qw));
}

} // anonymous namespace

// ─────────────────────────────────────────────────────────────────────────────
// Legacy API  —  cas::platform
// ─────────────────────────────────────────────────────────────────────────────
namespace cas::platform {

void applyMacOSChrome(QWindow* qw) {
    NSWindow* ns = ns_window_for(qw);
    if (!ns) return;

    // Unified title bar: content extends behind traffic lights.
    ns.titlebarAppearsTransparent = YES;
    ns.titleVisibility            = NSWindowTitleHidden;
    ns.styleMask                 |= NSWindowStyleMaskFullSizeContentView;
    ns.movableByWindowBackground  = YES;

    // Traffic lights: 12 pt from left edge, 16 pt from bottom of title bar
    // (NSWindow coordinates have origin at the bottom-left of the screen).
    const CGFloat kX = 12.0;
    const CGFloat kY = 16.0;
    const CGFloat kSpacing = 20.0;
    NSButton* buttons[] = {
        [ns standardWindowButton:NSWindowCloseButton],
        [ns standardWindowButton:NSWindowMiniaturizeButton],
        [ns standardWindowButton:NSWindowZoomButton]
    };
    for (NSUInteger i = 0; i < 3; ++i) {
        if (!buttons[i]) continue;
        NSRect f = buttons[i].frame;
        f.origin.x = kX + static_cast<CGFloat>(i) * kSpacing;
        f.origin.y = kY;
        buttons[i].frame = f;
    }

    // Full-window vibrancy: renders behind the Qt content view.
    NSVisualEffectView* fx =
        [[NSVisualEffectView alloc] initWithFrame:ns.contentView.bounds];
    fx.material      = NSVisualEffectMaterialUnderWindowBackground;
    fx.blendingMode  = NSVisualEffectBlendingModeBehindWindow;
    fx.state         = NSVisualEffectStateActive;
    fx.autoresizingMask = NSViewWidthSizable | NSViewHeightSizable;
    [ns.contentView addSubview:fx positioned:NSWindowBelow relativeTo:nil];
}

void observeAppearance(std::function<void(bool dark)> cb) {
    // The static local ensures we register only once per process lifetime.
    static id observer =
        [[NSDistributedNotificationCenter defaultCenter]
            addObserverForName:@"AppleInterfaceThemeChangedNotification"
                        object:nil
                         queue:[NSOperationQueue mainQueue]
                    usingBlock:^(NSNotification*) {
                NSString* style = [[NSUserDefaults standardUserDefaults]
                    stringForKey:@"AppleInterfaceStyle"];
                const bool dark = [style isEqualToString:@"Dark"];
                cb(dark);
            }];
    (void)observer; // suppress unused-variable warning
}

} // namespace cas::platform

// ─────────────────────────────────────────────────────────────────────────────
// Granular API  —  cas::gui::native
// ─────────────────────────────────────────────────────────────────────────────
namespace cas::gui::native {

// ── apply_vibrancy ────────────────────────────────────────────────────────────
// Inserts an NSVisualEffectView as the bottom-most subview of the NSWindow's
// content view.  The material chosen (HUDWindow) gives a frosted, opaque-ish
// background suitable for both light and dark appearances.
//
// Guards:
//   • No-op if ns_window_for() returns nil (window not yet on screen).
//   • No-op if a VisualEffectView with tag 9001 already exists (idempotent).
void apply_vibrancy(QWindow* qw) {
    NSWindow* ns = ns_window_for(qw);
    if (!ns) return;

    // Idempotency guard: tag 9001 is our sentinel.
    const NSInteger kTag = 9001;
    if ([ns.contentView viewWithTag:kTag]) return;

    NSVisualEffectView* fx =
        [[NSVisualEffectView alloc] initWithFrame:ns.contentView.bounds];
    fx.tag          = kTag;
    fx.material     = NSVisualEffectMaterialHUDWindow;
    fx.blendingMode = NSVisualEffectBlendingModeBehindWindow;
    fx.state        = NSVisualEffectStateActive;
    fx.autoresizingMask = NSViewWidthSizable | NSViewHeightSizable;

    [ns.contentView addSubview:fx positioned:NSWindowBelow relativeTo:nil];
}

// ── set_traffic_lights_pos ────────────────────────────────────────────────────
// Repositions the three standard window control buttons to the given (x, y)
// position in NSWindow's flipped coordinate system (origin = bottom-left).
// The three buttons are placed horizontally with a 20-pt spacing.
//
// Typical call:  set_traffic_lights_pos(window, 12.0, 16.0)
void set_traffic_lights_pos(QWindow* qw, double x, double y) {
    NSWindow* ns = ns_window_for(qw);
    if (!ns) return;

    const CGFloat spacing = 20.0;
    NSButton* buttons[] = {
        [ns standardWindowButton:NSWindowCloseButton],
        [ns standardWindowButton:NSWindowMiniaturizeButton],
        [ns standardWindowButton:NSWindowZoomButton]
    };
    for (NSUInteger i = 0; i < 3; ++i) {
        NSButton* btn = buttons[i];
        if (!btn) continue;
        NSRect f = btn.frame;
        f.origin.x = static_cast<CGFloat>(x) + static_cast<CGFloat>(i) * spacing;
        f.origin.y = static_cast<CGFloat>(y);
        btn.frame  = f;
    }
}

// ── set_title_bar_style ───────────────────────────────────────────────────────
// Configures the NSWindow for a unified, borderless title bar:
//   • titlebarAppearsTransparent → Qt content paints behind traffic lights.
//   • titleVisibility hidden     → no document title visible in the bar.
//   • FullSizeContentView mask   → content rect spans the entire window frame.
//   • movableByWindowBackground  → user can drag the window from any blank area.
void set_title_bar_style(QWindow* qw) {
    NSWindow* ns = ns_window_for(qw);
    if (!ns) return;

    ns.titlebarAppearsTransparent = YES;
    ns.titleVisibility            = NSWindowTitleHidden;
    ns.styleMask                 |= NSWindowStyleMaskFullSizeContentView;
    ns.movableByWindowBackground  = YES;
}

} // namespace cas::gui::native

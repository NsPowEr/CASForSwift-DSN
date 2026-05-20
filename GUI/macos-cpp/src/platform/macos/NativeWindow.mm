// NativeWindow.mm — bridging Cocoa per features non coperte da Qt.
// Si chiama da AppCore::AppCore() dopo che la finestra è creata.
//
//  - Vibrancy (NSVisualEffectView) sotto la sidebar e la title bar
//  - Traffic-light positioning custom (per allinearli alla title bar unificata)
//  - Toggle full-size content view
//  - Notification when system appearance changes (auto switch dark/light)

#include "NativeWindow.h"
#import <AppKit/AppKit.h>
#import <Cocoa/Cocoa.h>
#include <QWindow>
#include <QGuiApplication>
#include <QtGui/qpa/qplatformnativeinterface.h>

namespace cas::platform {

NSWindow* nsWindowFor(QWindow* qw) {
    return reinterpret_cast<NSWindow*>(
        QGuiApplication::platformNativeInterface()->nativeResourceForWindow("nswindow", qw));
}

void applyMacOSChrome(QWindow* qw) {
    NSWindow* ns = nsWindowFor(qw);
    if (!ns) return;

    // Title bar unificata + content sotto la barra
    ns.titlebarAppearsTransparent = YES;
    ns.titleVisibility = NSWindowTitleHidden;
    ns.styleMask |= NSWindowStyleMaskFullSizeContentView;
    ns.movableByWindowBackground = YES;

    // Traffic lights al posto giusto (12 px dall'angolo, centrati sui 44 px della titlebar)
    NSButton* close = [ns standardWindowButton:NSWindowCloseButton];
    NSButton* mini  = [ns standardWindowButton:NSWindowMiniaturizeButton];
    NSButton* zoom  = [ns standardWindowButton:NSWindowZoomButton];
    for (NSButton* b in @[close, mini, zoom]) {
        NSRect f = b.frame; f.origin.y = 16; b.frame = f;
    }

    // Vibrancy material come background della finestra
    NSVisualEffectView* fx = [[NSVisualEffectView alloc] initWithFrame:ns.contentView.bounds];
    fx.material  = NSVisualEffectMaterialUnderWindowBackground;
    fx.blendingMode = NSVisualEffectBlendingModeBehindWindow;
    fx.state = NSVisualEffectStateActive;
    fx.autoresizingMask = NSViewWidthSizable | NSViewHeightSizable;
    [ns.contentView addSubview:fx positioned:NSWindowBelow relativeTo:nil];
}

void observeAppearance(std::function<void(bool dark)> cb) {
    static id observer = [[NSDistributedNotificationCenter defaultCenter]
        addObserverForName:@"AppleInterfaceThemeChangedNotification"
                    object:nil queue:[NSOperationQueue mainQueue]
                usingBlock:^(NSNotification*) {
            NSString* style = [[NSUserDefaults standardUserDefaults]
                stringForKey:@"AppleInterfaceStyle"];
            const bool dark = [style isEqualToString:@"Dark"];
            cb(dark);
        }];
    (void)observer;
}

} // namespace

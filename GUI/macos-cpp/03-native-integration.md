# macOS C++ — Native Integration

Solo le feature non coperte da Qt richiedono Objective-C++ (`.mm`). Tutto il resto resta Qt puro.

## Cosa risolviamo in `platform/macos/`

| Feature | File | API Apple |
|---|---|---|
| Title bar unificata + vibrancy | `NativeWindow.mm` | `NSWindowStyleMaskFullSizeContentView`, `NSVisualEffectView` |
| Traffic-light positioning | `NativeWindow.mm` | `[NSWindow standardWindowButton:]` |
| Auto dark/light follow | `NativeWindow.mm` | `AppleInterfaceThemeChangedNotification` |
| Touch Bar (Intel Mac) | `TouchBar.mm` | `NSTouchBar`, `NSTouchBarItem` |
| Services menu | `ServicesProvider.mm` | `NSServicesProvider` |
| Quick Look preview celle | `QuickLook.mm` | `QLPreviewPanel` |
| Dock badge (job in coda) | `DockBadge.mm` | `NSDockTile` |
| Continuity Camera + Apple Pencil scribble | `ContinuityInput.mm` | `NSWritingTools`, `NKView` |
| Spotlight/Core Spotlight indexing | `SpotlightIndex.mm` | `CSSearchableIndex` |

## Pattern raccomandato

C++ side dichiara una funzione libera in un namespace `cas::platform`. L'implementazione `.mm` include AppKit. Un `.h` puro C++ è esposto al resto del codice — niente Qt né Cocoa nei consumer.

```cpp
// NativeWindow.h
namespace cas::platform {
    void applyMacOSChrome(QWindow* w);
    void observeAppearance(std::function<void(bool dark)>);
}
```

## Build CMake

```cmake
add_library(cas_platform_macos STATIC
    NativeWindow.mm NativeWindow.h
    TouchBar.mm     TouchBar.h
    # ...
)
set_source_files_properties(NativeWindow.mm TouchBar.mm
    PROPERTIES COMPILE_FLAGS "-fobjc-arc -fobjc-weak")
target_link_libraries(cas_platform_macos PUBLIC
    Qt6::Gui
    "-framework AppKit"
    "-framework Cocoa"
    "-framework CoreServices"
    "-framework QuartzCore")
```

## Code signing & notarization

`Info.plist.in` minimi:
- `LSApplicationCategoryType` = `public.app-category.education` (o `productivity`)
- `NSHighResolutionCapable` = YES
- `LSMinimumSystemVersion` = 12.0
- `NSSupportsAutomaticGraphicsSwitching` = YES

Entitlements:
- `com.apple.security.app-sandbox` = YES
- `com.apple.security.files.user-selected.read-write` = YES (per import/export PDF)
- `com.apple.security.network.client` = YES (kernel remoto)

CI script: `xcrun notarytool submit ... --wait` post-build.

## Accessibility

Qt fornisce `QAccessible` automaticamente. Per la matematica (cruciale): override di `MathRenderer::accessibilityValue()` con conversione **MathSpeak** dell'AST. Esempio: `\frac{x+1}{x-1}` → "fraction with numerator x plus 1 and denominator x minus 1".

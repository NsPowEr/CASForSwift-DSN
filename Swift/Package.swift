// swift-tools-version: 5.9
import PackageDescription

let package = Package(
    name: "CASEngine",
    platforms: [
        .iOS(.v16),
        .macOS(.v13)
    ],
    products: [
        .library(name: "CASEngine", targets: ["CASEngine"]),
    ],
    targets: [
        .target(
            name: "CCASEngine",
            path: "Sources/CCASEngine",
            exclude: [
                "src/ui/main.cpp"
            ],
            sources: [
                "src/algebra",
                "src/ast",
                "src/calculus",
                "src/capi",
                "src/formatter",
                "src/foundation",
                "src/lexer",
                "src/linalg",
                "src/numeric",
                "src/numtheory",
                "src/parser",
                "src/symbolic"
            ],
            publicHeadersPath: "include",
            cxxSettings: [
                .headerSearchPath("include"),
                .headerSearchPath("src"),
                .define("INTERNAL_CAS_BUILD")
            ]
        ),
        .target(
            name: "CASEngine",
            dependencies: ["CCASEngine"],
            path: "Sources/CASEngine"
        )
    ],
    cxxLanguageStandard: .cxx17
)

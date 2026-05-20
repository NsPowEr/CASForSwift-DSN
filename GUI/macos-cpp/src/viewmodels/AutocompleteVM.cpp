#include "AutocompleteVM.h"
#include <QVariantMap>
#include <algorithm>

namespace cas::gui {

AutocompleteVM::AutocompleteVM(QObject* parent) : QObject(parent) {
}

void AutocompleteVM::setPrefix(const QString& p) {
    if (m_prefix == p) return;
    m_prefix = p;
    emit prefixChanged();
    refreshSuggestions();
}

void AutocompleteVM::updateLibrary(const QStringList& functions, const QVariantList& variables) {
    m_library.clear();
    
    // Add functions
    for (const auto& fn : functions) {
        m_library.push_back({fn, "function", fn + "()", "Function"});
    }
    
    // Add variables
    for (const auto& v : variables) {
        QVariantMap vm = v.toMap();
        QString name = vm["name"].toString();
        const QString value = vm["value"].toString();
        m_library.push_back({name, "variable", name, value.isEmpty() ? "Variable" : value});
    }
    
    // Add common LaTeX templates and math snippets
    m_library.push_back({"frac", "template", "\\frac{}{}", "Fraction"});
    m_library.push_back({"sqrt", "template", "\\sqrt{}", "Square root"});
    m_library.push_back({"sum", "template", "\\sum_{}^{}", "Summation"});
    m_library.push_back({"int", "template", "\\int_{}^{} ", "Integral"});
    m_library.push_back({"lim", "template", "\\lim_{x \\to }", "Limit"});
    m_library.push_back({"alpha", "template", "\\alpha", "Greek letter"});
    m_library.push_back({"beta", "template", "\\beta", "Greek letter"});
    m_library.push_back({"theta", "template", "\\theta", "Greek letter"});
    m_library.push_back({"pi", "template", "\\pi", "Greek letter"});
    
    refreshSuggestions();
}

void AutocompleteVM::refreshSuggestions() {
    m_suggestions.clear();
    if (m_prefix.isEmpty()) {
        emit suggestionsChanged();
        emit visibleChanged();
        return;
    }

    struct ScoredItem {
        LibraryItem item;
        int score;
    };
    std::vector<ScoredItem> matches;
    const QString lookup_prefix = normalizedLookupKey(m_prefix);

    for (const auto& item : m_library) {
        const QString lookup_target = normalizedLookupKey(item.text);
        int score = fuzzyScore(lookup_prefix, lookup_target);
        if (score > 0) {
            matches.push_back({item, score});
        }
    }

    // Sort by score descending, then by length ascending
    std::sort(matches.begin(), matches.end(), [](const auto& a, const auto& b) {
        if (a.score != b.score) return a.score > b.score;
        if (a.item.type != b.item.type) {
            if (a.item.type == "function") return true;
            if (b.item.type == "function") return false;
            if (a.item.type == "variable") return true;
            if (b.item.type == "variable") return false;
        }
        return a.item.text.length() < b.item.text.length();
    });

    // Take top 8
    int limit = std::min<int>(8, matches.size());
    for (int i = 0; i < limit; ++i) {
        QVariantMap m;
        m["text"] = matches[i].item.text;
        m["type"] = matches[i].item.type;
        m["insertText"] = matches[i].item.insertText;
        m["detail"] = matches[i].item.detail;
        m_suggestions.append(m);
    }

    emit suggestionsChanged();
    emit visibleChanged();
}

int AutocompleteVM::fuzzyScore(const QString& prefix, const QString& target) const {
    if (prefix.isEmpty()) return 0;
    if (target == prefix) return 200;
    if (target.startsWith(prefix)) return 150 - (target.length() - prefix.length());
    if (target.contains(prefix)) return 90 - static_cast<int>(target.indexOf(prefix));
    
    // Simple subsequence check
    int pIdx = 0;
    int tIdx = 0;
    int matches = 0;
    while (pIdx < prefix.length() && tIdx < target.length()) {
        if (prefix[pIdx] == target[tIdx]) {
            pIdx++;
            matches++;
        }
        tIdx++;
    }
    
    if (pIdx == prefix.length()) {
        return 50 + matches; // Found all chars in order
    }
    
    return 0;
}

QString AutocompleteVM::normalizedLookupKey(const QString& value) const {
    QString normalized = value.trimmed().toLower();
    while (normalized.startsWith('\\')) {
        normalized.remove(0, 1);
    }
    return normalized;
}

} // namespace cas::gui

#include <gtest/gtest.h>

#include <QVariantMap>

#include "viewmodels/AutocompleteVM.h"

using cas::gui::AutocompleteVM;

TEST(GuiAutocomplete, PrefixAndLatexMatchingStayUseful) {
    AutocompleteVM vm;
    QVariantList variables;
    variables.append(QVariantMap{{"name", "sigma"}, {"value", "3"}});
    variables.append(QVariantMap{{"name", "sampleVar"}, {"value", "x+1"}});

    vm.updateLibrary(QStringList{"sin", "solve", "simplify"}, variables);

    vm.setPrefix("si");
    ASSERT_FALSE(vm.suggestions().isEmpty());
    const QVariantMap first_prefix = vm.suggestions().front().toMap();
    EXPECT_EQ(first_prefix.value("text").toString(), QString("sin"));
    EXPECT_EQ(first_prefix.value("type").toString(), QString("function"));
    EXPECT_EQ(first_prefix.value("insertText").toString(), QString("sin()"));

    vm.setPrefix("\\sq");
    ASSERT_FALSE(vm.suggestions().isEmpty());
    const QVariantMap first_latex = vm.suggestions().front().toMap();
    EXPECT_EQ(first_latex.value("text").toString(), QString("sqrt"));
    EXPECT_EQ(first_latex.value("type").toString(), QString("template"));
    EXPECT_EQ(first_latex.value("insertText").toString(), QString("\\sqrt{}"));
}

TEST(GuiAutocomplete, SuggestionsExposeDetailForVariablesAndTemplates) {
    AutocompleteVM vm;
    QVariantList variables;
    variables.append(QVariantMap{{"name", "sampleVar"}, {"value", "42"}});

    vm.updateLibrary(QStringList{"abs"}, variables);

    vm.setPrefix("sam");
    ASSERT_FALSE(vm.suggestions().isEmpty());
    const QVariantMap variable = vm.suggestions().front().toMap();
    EXPECT_EQ(variable.value("text").toString(), QString("sampleVar"));
    EXPECT_EQ(variable.value("type").toString(), QString("variable"));
    EXPECT_EQ(variable.value("detail").toString(), QString("42"));
    EXPECT_EQ(variable.value("insertText").toString(), QString("sampleVar"));

    vm.setPrefix("fra");
    ASSERT_FALSE(vm.suggestions().isEmpty());
    const QVariantMap templ = vm.suggestions().front().toMap();
    EXPECT_EQ(templ.value("text").toString(), QString("frac"));
    EXPECT_EQ(templ.value("type").toString(), QString("template"));
    EXPECT_EQ(templ.value("detail").toString(), QString("Fraction"));
    EXPECT_EQ(templ.value("insertText").toString(), QString("\\frac{}{}"));
}

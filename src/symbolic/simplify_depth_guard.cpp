#include "simplify_impl.hpp"

namespace cas::symbolic::detail {

thread_local int simplification_depth = 0;

DepthGuard::DepthGuard() { ++simplification_depth; }
DepthGuard::~DepthGuard() { --simplification_depth; }
bool DepthGuard::exceeded() const { return simplification_depth > MAX_SIMPLIFICATION_DEPTH; }

} // namespace cas::symbolic::detail

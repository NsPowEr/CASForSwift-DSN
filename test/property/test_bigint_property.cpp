// F0.4 — BigInt Property Tests: Knuth D, Toom-3.
// Identity: (a*b)/b == a  and  a/b*b + a%b == a.

#include "cas/bigint.hpp"
#include <gtest/gtest.h>
#include <rapidcheck/gtest.h>
#include <vector>
#include <cstdint>

using namespace cas;

namespace {

RC_GTEST_PROP(BigIntPropertyTest, KnuthD_DivisionIdentity, 
              (const std::vector<std::uint32_t>& alimbs, 
               const std::vector<std::uint32_t>& blimbs)) {
    BigInt a = BigInt::from_limbs_le(alimbs);
    BigInt b = BigInt::from_limbs_le(blimbs);
    
    if (b.is_zero()) return;
    
    BigInt q = a / b;
    BigInt r = a % b;
    
    // a = q*b + r
    BigInt check = q * b + r;
    RC_ASSERT(check.decimal() == a.decimal());
    RC_ASSERT(r.abs() < b.abs());
}

RC_GTEST_PROP(BigIntPropertyTest, MultiLimb_KaratsubaToom3_Roundtrip, 
              (const std::vector<std::uint32_t>& alimbs, 
               const std::vector<std::uint32_t>& blimbs)) {
    BigInt a = BigInt::from_limbs_le(alimbs);
    BigInt b = BigInt::from_limbs_le(blimbs);
    
    if (b.is_zero() || a.is_zero()) return;
    
    BigInt prod = a * b;
    RC_ASSERT((prod / b).decimal() == a.decimal());
}

} // namespace

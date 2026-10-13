#pragma once
// ============================================================================
// scompute.hpp — 核心数值：从两个 53 位均匀数协同生成 S = -ln(u)
// ----------------------------------------------------------------------------
// 把两个 53 位数拼成 106 位整数：
//   U = u_hi_bits * 2^53 + u_lo_bits            ∈ [0, 2^106)
// 归一化均匀数：
//   u = (U + 1) * 2^-106                        ∈ (0, 1]
//   -ln u = 106·ln2 - ln(U + 1)
//
// 不能直接算 U（106 位整数 long double 存不下），分两种情况：
//   情况 A：u_hi_bits > 0
//     U + 1 = u_hi_bits·2^53 · (1 + (u_lo_bits+1)/(u_hi_bits·2^53))
//     S = 53·ln2 - ln(u_hi_bits) - ln(1 + δ)
//     其中 δ = (u_lo_bits+1)/(u_hi_bits·2^53)，用 log1pl 精确计算
//   情况 B：u_hi_bits == 0
//     U + 1 = u_lo_bits + 1
//     -ln u = 106·ln2 - ln(u_lo_bits + 1)
// ============================================================================

#include <cmath>
#include <cstdint>

namespace rs {

// 80 位精度 ln2（hex float 字面量，避免十进制舍入）
inline constexpr long double LN2 = 0x1.62e42fefa39efp-1L;

// 精算 S：需要两个 53 位块
inline long double S_from_two(std::uint64_t u_hi_bits,
                              std::uint64_t u_lo_bits)
{
    if (u_hi_bits == 0) {
        long double lo = (long double)(u_lo_bits + 1);   // [1, 2^53]
        return 106.0L * LN2 - std::logl(lo);
    } else {
        long double hi = (long double)u_hi_bits;
        long double S  = 53.0L * LN2 - std::logl(hi);
        long double delta = ((long double)u_lo_bits + 1.0L) / hi * 0x1p-53L;
        return S - std::log1pl(delta);
    }
}

// 粗估 S：只看高位块，结果是保守上界（>= 真实 S），
// 用于快速拒绝：若上界都低于 1σ 阈值，真实 S 必然更低，可直接丢弃。
inline long double S_from_hi(std::uint64_t u_hi_bits)
{
    if (u_hi_bits == 0) return 106.0L * LN2;   // 真实 S ≤ 106·ln2（由 u_lo 决定）
    return 53.0L * LN2 - std::logl((long double)u_hi_bits);
}

} // namespace rs

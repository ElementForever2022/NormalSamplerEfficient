#pragma once
// ============================================================================
// scompute.hpp — 核心数值：从两个 53 位均匀数协同生成"双边尾坐标"
// ----------------------------------------------------------------------------
// 均匀数（106 位）：
//   U = u_hi_bits * 2^53 + u_lo_bits            ∈ [0, 2^106)
//   u = (U + 1) * 2^-106                        ∈ (0, 1]
//
// 由 u 反解标准正态分位 x = Phi_bar^{-1}(u) = Phi^{-1}(1-u) 可知 x ~ N(0,1)
// 是**带符号**的：u <= 1/2 给出右尾 x >= 0，u > 1/2 给出左尾 x < 0。
// 分档必须按 |x| 双边来做，于是定义统一坐标
//
//   y = -ln(min(u, 1-u)) >= ln2
//
// 它同时覆盖两条尾：
//   |x| >= kσ  ⟺  Phi_bar(|x|) <= Phi_bar(k)  ⟺  min(u,1-u) <= Phi_bar(k)
//              ⟺  y >= SH[k] = -ln(Phi_bar(k))
// 即左右尾共用同一张阈值表 SH，这是本次分档修正的关键。
//
// 两条尾各自的 S 值（S = -ln(u)，与旧版一致，便于日志对照）：
//   右尾：S      = -ln(u)      = 106·ln2 - ln(U + 1)
//   左尾：y_L    = -ln(1-u)    = 106·ln2 - ln(c),  c = 2^106 - U - 1
//
// 注意：S_from_two 只在右尾（u <= 1/2）有意义，S_left_from_two 只在左尾
// 有意义；在另一条尾上调用的结果是 ~1e-15 量级的舍入残渣（73 量级相减相消），
// 不要拿去比较阈值。统一入口用 Y_from_bits() / Y_upper_from_hi()。
//
// 不能直接算 U / c（106 位整数 long double 存不下），分两种情况（以右尾为例）：
//   情况 A：u_hi_bits > 0
//     U + 1 = u_hi_bits·2^53 · (1 + (u_lo_bits+1)/(u_hi_bits·2^53))
//     S = 53·ln2 - ln(u_hi_bits) - ln(1 + δ)
//     其中 δ = (u_lo_bits+1)/(u_hi_bits·2^53)，用 log1p 精确计算
//   情况 B：u_hi_bits == 0
//     U + 1 = u_lo_bits + 1
//     S = 106·ln2 - ln(u_lo_bits + 1)
// 左尾的 c = c_hi·2^53 + c_lo（c_hi = 2^53-1-u_hi_bits，c_lo = 2^53-1-u_lo_bits，
// 无进位），与情况 A/B 同构，只是不含 "+1"。
// ============================================================================

#include <cmath>
#include <cstdint>

namespace rs {

// ln2，64 位尾数的 hex float 字面量（80 位 long double 下与 std::log(2.0L) 同值；
// 旧的 0x1.62e42fefa39efp-1L 只带 53 位，会给 106·ln2 引入 ~2.4e-15 的绝对误差）
inline constexpr long double LN2 = 0x1.62e42fefa39ef357p-1L;

// u <= 1/2 ⟺ U + 1 <= 2^105 ⟺ u_hi_bits < 2^52（u_hi_bits == 2^52 时 U+1 > 2^105）
inline constexpr std::uint64_t RIGHT_TAIL_LIMIT = std::uint64_t(1) << 52;

// 53 位块全 1，用于取补（c_hi / c_lo 的基数 2^53 - 1）
inline constexpr std::uint64_t MASK53 = (std::uint64_t(1) << 53) - 1;

// 106 位坐标的上界常量：106·ln2，也是 c = 1（u 达到最大可表示值）时的坐标
inline constexpr long double Y_MAX = 106.0L * LN2;

// 快速拒绝的安全余量：粗估上界自身带 ~1e-18 舍入误差，为了严格保证
// "上界确实低于阈值才拒绝"，拒绝判据用 bound < SH[1] - REJECT_MARGIN。
// 代价：只有落在 1σ 边界 1e-15 邻域内的样本改走精算路径，可忽略不计。
inline constexpr long double REJECT_MARGIN = 1e-15L;

// u 是否落在右半（x >= 0）：u <= 1/2
inline bool is_right_tail(std::uint64_t u_hi_bits)
{
    return u_hi_bits < RIGHT_TAIL_LIMIT;
}

// 精算右尾坐标 S = -ln(u)：需要两个 53 位块（仅当 u <= 1/2 时使用）
// 注：用 std::log/std::log1p 的 long double 重载（C++ 标准保证），
// 不用 std::logl/std::log1pl（部分 libstdc++ 不在 std 内提供 C 包装）。
inline long double S_from_two(std::uint64_t u_hi_bits,
                              std::uint64_t u_lo_bits)
{
    if (u_hi_bits == 0) {
        long double lo = (long double)(u_lo_bits + 1);   // [1, 2^53]
        return 106.0L * LN2 - std::log(lo);
    } else {
        long double hi = (long double)u_hi_bits;
        long double S  = 53.0L * LN2 - std::log(hi);
        long double delta = ((long double)u_lo_bits + 1.0L) / hi * 0x1p-53L;
        return S - std::log1p(delta);
    }
}

// 精算左尾坐标 y_L = -ln(1-u)：需要两个 53 位块（仅当 u > 1/2 时使用）
// 1-u = c·2^-106, c = (2^53-1-u_hi_bits)·2^53 + (2^53-1-u_lo_bits)
// 直接对补数取对数，避免 1-u 在 u→1 时的灾难性相消。
// 特例 u == 1（c == 0，概率 2^-106，实际不可达）：|x| = ∞，坐标封顶到 Y_MAX。
inline long double S_left_from_two(std::uint64_t u_hi_bits,
                                   std::uint64_t u_lo_bits)
{
    std::uint64_t c_hi = MASK53 - u_hi_bits;
    std::uint64_t c_lo = MASK53 - u_lo_bits;
    if (c_hi == 0) {
        if (c_lo == 0) return Y_MAX;                     // u == 1：封顶
        return 106.0L * LN2 - std::log((long double)c_lo);
    }
    long double chi   = (long double)c_hi;
    long double delta = ((long double)c_lo / chi) * 0x1p-53L;
    return 53.0L * LN2 - std::log(chi) - std::log1p(delta);
}

// 粗估右尾坐标 S 的上界（>= 真实 S），只看高位块。
inline long double S_from_hi(std::uint64_t u_hi_bits)
{
    if (u_hi_bits == 0) return 106.0L * LN2;   // 真实 S <= 106·ln2（由 u_lo 决定）
    return 53.0L * LN2 - std::log((long double)u_hi_bits);
}

// 粗估左尾坐标 y_L 的上界（>= 真实 y_L），只看高位块：
//   c >= c_hi·2^53  ⟹  y_L = 106·ln2 - ln c <= 53·ln2 - ln(c_hi)
// c_hi == 0（即 u 只差最后 2^-53 就到 1）时用全局上界 Y_MAX 兜底。
inline long double S_left_from_hi(std::uint64_t u_hi_bits)
{
    std::uint64_t c_hi = MASK53 - u_hi_bits;
    if (c_hi == 0) return Y_MAX;
    return 53.0L * LN2 - std::log((long double)c_hi);
}

// ---- 统一入口：自动按 u 选左右尾，避免误用单尾函数 ----

// 双边尾坐标 y = -ln(min(u, 1-u)) = max(-ln u, -ln(1-u)) >= ln2
inline long double Y_from_bits(std::uint64_t u_hi_bits,
                               std::uint64_t u_lo_bits)
{
    return is_right_tail(u_hi_bits) ? S_from_two(u_hi_bits, u_lo_bits)
                                    : S_left_from_two(u_hi_bits, u_lo_bits);
}

// 双边尾坐标的上界（>= 真实 y），只看高位块，用于快速拒绝
inline long double Y_upper_from_hi(std::uint64_t u_hi_bits)
{
    return is_right_tail(u_hi_bits) ? S_from_hi(u_hi_bits)
                                    : S_left_from_hi(u_hi_bits);
}

} // namespace rs

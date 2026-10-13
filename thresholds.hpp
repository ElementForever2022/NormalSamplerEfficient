#pragma once
// ============================================================================
// thresholds.hpp — σ 档位阈值表（双边：±kσ 各自一个阈值）
// ----------------------------------------------------------------------------
// 对每个档位 k = 0..NBINS：
//   TH[k] : σ 值 k
//   SH[k] : 右尾阈值，SH[k] = -ln(Phi_bar(k))，Phi_bar(k) = erfc(k/√2)/2
//           S = -ln(u) >= SH[k]  ⟺  x >= +kσ
//   LH[k] : 左尾阈值，LH[k] = -ln(Phi(k)) = -ln(1 - Phi_bar(k))
//           S = -ln(u) <= LH[k]  ⟺  x <= -kσ
//
// 用统一双边坐标 y = -ln(min(u, 1-u))（恒有 y >= ln2 = SH[0] = LH[0]）表述：
//   |x| >= kσ  ⟺  y >= SH[k]        ← 左右尾共用 SH，分档只与 SH 比较
//   x >= +kσ   ⟺  S >= SH[k]        ← 右尾（S = -ln u）
//   x <= -kσ   ⟺  S <= LH[k]        ← 左尾
// 所以 SH[k] 与 LH[k] 就是档位 k 的"上下界"：低于 SH[k]（且高于 LH[k]）
// 才属于 |x| < kσ，不会把左尾的样本算进右尾的区间。
//
// SH/LH 由 init_thresholds() 在启动时用 long double + erfc 精确计算，不硬编码。
// ============================================================================

#include "config.hpp"

namespace rs {

extern double      TH[NBINS + 1];   // [0] = 0σ（占位，非阈值）
extern long double SH[NBINS + 1];   // 右尾阈值 = -ln(Phi_bar(k))，SH[0] = ln2
extern long double LH[NBINS + 1];   // 左尾阈值 = -ln(Phi(k))，LH[0] = ln2

// 给定 x > 0，求 S = -ln(Phi_bar(x))（单调递增）
long double s_from_x(long double x);

// 给定双边坐标 y >= ln2，反解 |x| >= 0（二分求解，只在记录路径上调用）
long double x_from_s(long double y);

// 由双边坐标 y 取档位：返回满足 y >= SH[k] 的最大 k，即 |x| 落在 [k, k+1)σ
// （取"最大 k"这一步同时施加了上界 SH[k+1]：若 y < SH[k+1] 就不会取到 k+1，
//   所以档位判定天然是"下界 SH[k] <= y < SH[k+1] 上界"的双边判定）
inline int ring_of_y(long double y)
{
    for (int i = NBINS; i >= 1; --i)
        if (y >= SH[i]) return i;
    return 0;   // y >= SH[0] = ln2 恒成立，故 0 档（|x| < 1σ）兜底
}

// 启动时调用一次，填充 TH / SH / LH
void init_thresholds();

} // namespace rs

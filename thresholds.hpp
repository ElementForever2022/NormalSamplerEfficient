#pragma once
// ============================================================================
// thresholds.hpp — sigma 档位阈值表
// ----------------------------------------------------------------------------
// TH[i] : 第 i 档的 sigma 值（i = 1..NBINS）
// SH[i] : 对应的 S 阈值，S = -ln(Phi_bar(x))
//         Phi_bar(x) = P(N(0,1) > x) = erfc(x/√2)/2
//
// SH 由 init_thresholds() 在启动时用 long double + erfc 精确计算，不硬编码。
// ============================================================================

#include "config.hpp"

namespace rs {

extern double      TH[NBINS + 1];   // [0] 占位
extern long double SH[NBINS + 1];   // [0] 占位

// 给定 x > 0，求 S = -ln(Phi_bar(x))
long double s_from_x(long double x);

// 给定 S，反解 x（牛顿迭代，仅用于日志记录，不在主循环热路径上）
long double x_from_s(long double S);

// 启动时调用一次，填充 TH / SH
void init_thresholds();

} // namespace rs

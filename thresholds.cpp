// thresholds.cpp — 阈值表实现（双边 ±kσ）
#include "thresholds.hpp"

#include <cmath>

#include "scompute.hpp"

namespace rs {

double      TH[NBINS + 1] = {};
long double SH[NBINS + 1] = {};
long double LH[NBINS + 1] = {};

// 注：用 std::log/std::erfc/std::sqrt/std::acos 的 long double 重载
// （C++ 标准保证），不用 std::logl 等 C 包装（部分 libstdc++ 不在 std 内提供）。
long double s_from_x(long double x)
{
    long double p = 0.5L * std::erfc(x / std::sqrt(2.0L));
    return -std::log(p);
}

// 反解：给定双边坐标 y >= ln2，求 |x|（s_from_x 单调递增）
// y = ln2 ⟺ |x| = 0；y < ln2 无解，返回 0。
// 只在"刷新纪录"时调用，不在主循环热路径上，故用二分求稳（4 步牛顿在
// x→0 附近会因渐近导数 x + 1/x 失真而失败）。
long double x_from_s(long double y)
{
    if (y <= LN2) return 0.0L;

    long double hi = 1.0L;
    while (s_from_x(hi) < y) hi *= 2.0L;      // 倍增出上界：s_from_x(hi) >= y
    long double lo = 0.0L;                    // s_from_x(0) = ln2 < y

    for (int i = 0; i < 100; ++i) {
        long double mid = 0.5L * (lo + hi);
        if (s_from_x(mid) < y) lo = mid;
        else                   hi = mid;
    }
    return 0.5L * (lo + hi);
}

void init_thresholds()
{
    for (int i = 0; i <= NBINS; ++i) {
        TH[i] = (double)i;
        // 右尾阈值；i = 0 时为 -ln(1/2) = ln2，即双边坐标的下界
        SH[i] = s_from_x((long double)i);
        // 左尾阈值 = -ln(Phi(i)) = -ln(1 - Phi_bar(i))
        // 用 log1p 处理 1 - e^{-SH} 在 i 较大时的相消误差
        LH[i] = -std::log1p(-std::exp(-SH[i]));
    }
}

} // namespace rs

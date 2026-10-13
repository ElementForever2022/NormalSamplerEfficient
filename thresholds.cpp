// thresholds.cpp — 阈值表实现
#include "thresholds.hpp"

#include <cmath>

namespace rs {

double      TH[NBINS + 1] = {};
long double SH[NBINS + 1] = {};

long double s_from_x(long double x)
{
    long double p = 0.5L * std::erfcl(x / std::sqrtl(2.0L));
    return -std::logl(p);
}

long double x_from_s(long double S)
{
    long double t   = 2.0L * S;
    long double l2pi = std::logl(2.0L * std::acosl(-1.0L));
    long double x   = std::sqrtl(t - std::logl(t) - l2pi + 1.0L);
    if (x < 1.0L) x = 1.0L;
    for (int i = 0; i < 4; ++i) {
        long double f  = s_from_x(x) - S;
        long double df = x + 1.0L / x;
        x -= f / df;
        if (x < 1.0L) x = 1.0L;
    }
    return x;
}

void init_thresholds()
{
    for (int i = 0; i <= NBINS; ++i) {
        TH[i] = (double)i;
        SH[i] = (i == 0) ? 0.0L : s_from_x((long double)i);
    }
}

} // namespace rs

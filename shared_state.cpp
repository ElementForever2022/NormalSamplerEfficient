// shared_state.cpp — 共享状态定义与工具函数
#include "shared_state.hpp"

#include <ctime>

namespace rs {

std::mutex              g_mutex;
std::condition_variable g_cv;
std::atomic<double>     g_best_S{0.0};
Record                  g_best;
std::atomic<bool>       g_stop{false};
Record                  g_final;
std::atomic<bool>       g_reached{false};
std::atomic<bool>       g_has_new_record{false};

std::atomic<SampleCount> g_total_samples{0};

std::mutex  g_bin_mutex;
SampleCount g_tail[NBINS + 1] = {0};
SampleCount g_ring[NBINS + 1] = {0};

void flush_bins(LocalBins& lb)
{
    std::lock_guard<std::mutex> lk(g_bin_mutex);
    for (int i = 0; i <= NBINS; ++i) {
        if (lb.tail[i]) { g_tail[i] += lb.tail[i]; lb.tail[i] = 0; }
        if (lb.ring[i]) { g_ring[i] += lb.ring[i]; lb.ring[i] = 0; }
    }
}

std::string format_gmt(std::chrono::system_clock::time_point tp)
{
    std::time_t t = std::chrono::system_clock::to_time_t(tp);
    std::tm tm_buf{};
#if defined(_WIN32)
    gmtime_s(&tm_buf, &t);
#else
    gmtime_r(&t, &tm_buf);
#endif
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", &tm_buf);
    return std::string("GMT ") + buf;
}

} // namespace rs

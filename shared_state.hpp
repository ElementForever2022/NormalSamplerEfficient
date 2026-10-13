#pragma once
// ============================================================================
// shared_state.hpp — 线程间共享状态
// ----------------------------------------------------------------------------
// 所有全局可变状态集中在这里，排查并发问题时只需看这一个文件：
//   g_mutex / g_best       —— 最稀有纪录（保护它们用 g_mutex）
//   g_cv                   —— 通知记录线程"有新纪录 / 要停止"
//   g_best_S               —— g_best.S 的原子镜像，采样线程无锁快筛
//   g_stop / g_reached     —— 停止标志 / 是否达标
//   g_total_samples        —— 全局总样本数（批量原子更新）
//   g_ring                 —— 全局环形分箱计数（保护它用 g_bin_mutex）
//
// 锁顺序约定（避免死锁）：g_mutex -> g_bin_mutex，不可反向。
// ============================================================================

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <mutex>
#include <string>

#include "config.hpp"

namespace rs {

// 一条稀有纪录
struct Record {
    long double S = 0.0L;
    long double x = 0.0L;
    SampleCount sample_index = 0;
    std::chrono::system_clock::time_point when;
    bool        valid = false;
};

// ---- 纪录相关（g_mutex 保护）----
extern std::mutex              g_mutex;
extern std::condition_variable g_cv;
extern std::atomic<double>     g_best_S;    // double 而非 long double：避免 16 字节原子
extern Record                  g_best;
extern std::atomic<bool>       g_stop;
extern Record                  g_final;
extern std::atomic<bool>       g_reached;
extern std::atomic<bool>       g_has_new_record;

// ---- 计数相关 ----
extern std::atomic<SampleCount> g_total_samples;

// ---- 分箱计数（g_bin_mutex 保护）----
extern std::mutex  g_bin_mutex;
extern SampleCount g_ring[NBINS + 1];

// 线程本地分箱，攒够 FLUSH_INTERVAL 再并入全局，降低锁竞争
struct LocalBins {
    SampleCount ring[NBINS + 1] = {0};
};

// 把本地分箱并入全局并清零
void flush_bins(LocalBins& lb);

// 时间点格式化为 "GMT YYYY-MM-DD HH:MM:SS"
std::string format_gmt(std::chrono::system_clock::time_point tp);

} // namespace rs

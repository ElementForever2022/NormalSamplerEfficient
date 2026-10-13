#pragma once
// ============================================================================
// config.hpp — 全局配置常量
// ----------------------------------------------------------------------------
// 修改目标 sigma、档位数、记录间隔等都只需改这一个文件。
// ============================================================================

#include <cstdint>

namespace rs {

// 目标：找到 S = -ln(u) >= SH[TARGET_SIGMA] 的样本即停止
// 等价于 |x| >= TARGET_SIGMA σ 的标准正态尾部事件
static constexpr double TARGET_SIGMA       = 10.0;

// 记录线程每隔多少秒写一次"每分钟"快照（新纪录会立即写，不受此限制）
static constexpr int RECORD_INTERVAL_SEC   = 60;

// 统计档位数：1σ, 2σ, ..., NBINS σ
static constexpr int NBINS                 = 12;

// 默认线程数，0 = 取 hardware_concurrency()
static constexpr int NTHREADS_DEFAULT      = 0;

// 样本计数类型，必须 64 位（总样本数可达 10^23 量级）
using SampleCount = std::uint64_t;
static_assert(sizeof(SampleCount) >= 8, "样本计数必须至少 64 位");

// 每累计多少个样本把本地分箱刷到全局（降低锁竞争）
static constexpr SampleCount FLUSH_INTERVAL = 1ull << 16;

// 每累计多少个样本更新一次全局总样本数（原子计数，减少 cache 抖动）
static constexpr SampleCount TOTAL_BATCH    = 1024;

} // namespace rs

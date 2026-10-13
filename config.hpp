#pragma once
// ============================================================================
// config.hpp — 全局配置常量
// ----------------------------------------------------------------------------
// 修改目标 sigma、档位数、记录间隔等都只需改这一个文件。
// ============================================================================

#include <cstdint>

namespace rs {

// 目标：找到双边 |x| >= TARGET_SIGMA σ 的样本即停止
// 等价于双边尾概率 2·Phi_bar(TARGET_SIGMA)；分档/纪录也全部按 |x| 双边口径。
static constexpr double TARGET_SIGMA       = 10.0;

// 目标必须是整数档位（阈值表按 1σ, 2σ, ... 离散给出）
static_assert(TARGET_SIGMA == (double)(int)TARGET_SIGMA,
              "TARGET_SIGMA 必须是整数（阈值表只离散到整数 σ）");

// 记录线程每隔多少秒写一次"每分钟"快照（新纪录会立即写，不受此限制）
static constexpr int RECORD_INTERVAL_SEC   = 60;

// 统计档位数：|x| ∈ [0,1)σ, [1,2)σ, ..., [NBINS, ∞)σ（双边，左右尾合计）
static constexpr int NBINS                 = 12;

// 最后一档（[NBINS, ∞)）要能覆盖目标 σ，否则达标样本无处分档
static_assert(NBINS >= (int)TARGET_SIGMA, "NBINS 必须 >= TARGET_SIGMA");

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

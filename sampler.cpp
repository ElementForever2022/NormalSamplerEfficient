// sampler.cpp — 采样线程实现
#include "sampler.hpp"

#include "config.hpp"
#include "rng.hpp"
#include "scompute.hpp"
#include "shared_state.hpp"
#include "thresholds.hpp"

namespace rs {

void sample_loop(std::uint64_t seed)
{
    Xoshiro256pp rng(seed);

    LocalBins   lb;
    SampleCount idx        = 0;
    SampleCount last_flush = 0;

    const long double SH1 = SH[1];                          // 1σ 阈值：快速拒绝线
    const long double SHT = SH[NBINS >= 10 ? 10 : NBINS];  // 目标阈值（默认 10σ）

    while (!g_stop.load(std::memory_order_relaxed)) {
        std::uint64_t u_hi_bits = rng.next() >> 11;   // 53 位
        std::uint64_t u_lo_bits = rng.next() >> 11;   // 53 位

        ++idx;
        if ((idx % TOTAL_BATCH) == 0)
            g_total_samples.fetch_add(TOTAL_BATCH, std::memory_order_relaxed);

        // ---- 级 1：主项粗估（保守上界）----
        long double S_est = S_from_hi(u_hi_bits);

        // ---- 快速出口：< 1σ 直接丢弃（上界都不到 1σ，真实值必然不到）----
        if (S_est < SH1) {
            lb.ring[0] += 1;
            if (idx - last_flush >= FLUSH_INTERVAL) {
                flush_bins(lb);
                last_flush = idx;
            }
            continue;
        }

        // ---- 级 2：精算 ----
        long double S_full = S_from_two(u_hi_bits, u_lo_bits);

        // ---- 分档：找满足 S >= SH[i] 的最大 i ----
        int ring = 0;
        for (int i = NBINS; i >= 1; --i) {
            if (S_full >= SH[i]) { ring = i; break; }
        }
        lb.ring[ring] += 1;
        for (int j = ring; j >= 1; --j) lb.tail[j] += 1;

        // ---- 更新最稀有纪录（双重检查：先无锁读原子镜像，再加锁写）----
        if (S_full > (long double)g_best_S.load(std::memory_order_relaxed)) {
            bool is_new = false;
            {
                std::lock_guard<std::mutex> lk(g_mutex);
                if (!g_best.valid || S_full > g_best.S) {
                    g_best.S            = S_full;
                    g_best.x            = x_from_s(S_full);
                    g_best.sample_index = idx;
                    g_best.when         = std::chrono::system_clock::now();
                    g_best.valid        = true;
                    g_best_S.store((double)S_full, std::memory_order_relaxed);
                    // 在锁内置位，防止记录线程在 wait 判定与真正入睡之间错过通知
                    g_has_new_record.store(true, std::memory_order_relaxed);
                    is_new = true;
                }
            }
            if (is_new) {
                g_cv.notify_all();
            }
        }

        // ---- 达标停止 ----
        if (S_full >= SHT) {
            {
                std::lock_guard<std::mutex> lk(g_mutex);
                g_final = g_best;
                g_reached.store(true, std::memory_order_relaxed);
                g_stop.store(true, std::memory_order_relaxed);
            }
            flush_bins(lb);
            g_total_samples.fetch_add(idx % TOTAL_BATCH,
                                      std::memory_order_relaxed);
            g_cv.notify_all();
            return;
        }

        if (idx - last_flush >= FLUSH_INTERVAL) {
            flush_bins(lb);
            last_flush = idx;
        }
    }

    // 被别的线程停机：把自己攒下的计数刷干净再退出
    flush_bins(lb);
    SampleCount rem = idx % TOTAL_BATCH;
    if (rem) g_total_samples.fetch_add(rem, std::memory_order_relaxed);
}

} // namespace rs

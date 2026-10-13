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

    const long double SH1 = SH[1];                    // 1σ 阈值：快速拒绝线
    const long double SHT = SH[(int)TARGET_SIGMA > NBINS ? NBINS
                                                         : (int)TARGET_SIGMA];

    while (!g_stop.load(std::memory_order_relaxed)) {
        std::uint64_t u_hi_bits = rng.next() >> 11;   // 53 位
        std::uint64_t u_lo_bits = rng.next() >> 11;   // 53 位

        ++idx;
        if ((idx % TOTAL_BATCH) == 0)
            g_total_samples.fetch_add(TOTAL_BATCH, std::memory_order_relaxed);

        // ---- 级 0：定左右尾 ----
        // x = Phi_bar^{-1}(u) ~ N(0,1) 是带符号的：u <= 1/2 → 右尾 x >= 0，
        // u > 1/2 → 左尾 x < 0。分档按 |x| 双边，两条尾都要保留。
        const bool right = is_right_tail(u_hi_bits);

        // ---- 级 1：粗估双边坐标的上界，上界不到 1σ 才丢弃 ----
        // 留 REJECT_MARGIN 余量，保证"上界 < 1σ 阈值"严格成立，不会误杀。
        if (Y_upper_from_hi(u_hi_bits) < SH1 - REJECT_MARGIN) {
            lb.ring[0] += 1;
            if (idx - last_flush >= FLUSH_INTERVAL) {
                flush_bins(lb);
                last_flush = idx;
            }
            continue;
        }

        // ---- 级 2：精算双边坐标 y = -ln(min(u, 1-u)) ----
        // 右尾 y = -ln(u)，左尾 y = -ln(1-u)（直接对补数取对数，无相消误差）。
        const long double y = Y_from_bits(u_hi_bits, u_lo_bits);

        // ---- 分档：|x| ∈ [k, k+1)σ，上下界即 SH[k] / SH[k+1]（±kσ、±(k+1)σ）----
        const int ring = ring_of_y(y);
        lb.ring[ring] += 1;

        // ---- 更新最稀有纪录（双重检查：先无锁读原子镜像，再加锁写）----
        if (y > (long double)g_best_y.load(std::memory_order_relaxed)) {
            bool is_new = false;
            {
                std::lock_guard<std::mutex> lk(g_mutex);
                if (!g_best.valid || y > g_best.y) {
                    g_best.y            = y;
                    g_best.x            = right ? x_from_s(y) : -x_from_s(y);
                    g_best.sample_index = idx;
                    g_best.when         = std::chrono::system_clock::now();
                    g_best.valid        = true;
                    g_best_y.store((double)y, std::memory_order_relaxed);
                    // 在锁内置位，防止记录线程在 wait 判定与真正入睡之间错过通知
                    g_has_new_record.store(true, std::memory_order_relaxed);
                    is_new = true;
                }
            }
            if (is_new) {
                g_cv.notify_all();
            }
        }

        // ---- 达标停止：|x| >= TARGET_SIGMA σ（双边，左右尾任一方向都算）----
        if (y >= SHT) {
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

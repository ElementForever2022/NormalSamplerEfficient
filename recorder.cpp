// recorder.cpp — 记录线程实现
#include "recorder.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "config.hpp"
#include "scompute.hpp"
#include "shared_state.hpp"
#include "thresholds.hpp"

namespace rs {

void record_loop(const std::string& filename)
{
    std::ofstream ofs(filename, std::ios::out | std::ios::app);
    if (!ofs) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return;
    }

    const long double LN10 = std::log(10.0L);

    // ---- 内部辅助：把当前全局分箱写入文件 ----
    // 紧凑格式：每类分布压缩到 2 行，bin 间用 " | " 分隔。
    // 分档是**双边** |x|：bin k 表示 |x| ∈ [k, k+1)σ（左尾 x<=-k 与右尾 x>=k 合计），
    // 判定上下界即 ±kσ / ±(k+1)σ 对应的 SH 阈值（见 thresholds.hpp）。
    auto write_bins = [&](const char* indent) {
        SampleCount ring[NBINS + 1];
        {
            std::lock_guard<std::mutex> lk(g_bin_mutex);
            for (int i = 0; i <= NBINS; ++i) ring[i] = g_ring[i];
        }

        const int half = (NBINS + 1) / 2;   // 每行放一半 bin

        // 环形分布（观测计数）：2 行
        ofs << indent << "环形分布 (|x| 落在 [k, k+1)σ, 双边):";
        for (int i = 0; i <= NBINS; ++i) {
            if (i == half) ofs << "\n" << indent << "                ";
            else if (i > 0) ofs << " | ";
            ofs << "[" << std::setprecision(2) << TH[i] << ",";
            if (i < NBINS) ofs << std::setprecision(2) << TH[i + 1];
            else           ofs << "∞";
            ofs << ")σ:" << static_cast<unsigned long long>(ring[i]);
        }
        ofs << "\n";

        // 各档理论概率的 log10（稀有度）：2 行，与上面对齐
        // 双边：P(|x| ∈ [k,k+1)) = 2(Phi_bar(k) - Phi_bar(k+1))，末档 = 2·Phi_bar(NBINS)
        ofs << indent << "稀有度 log10 P: ";
        for (int i = 0; i <= NBINS; ++i) {
            if (i == half) ofs << "\n" << indent << "                ";
            else if (i > 0) ofs << " | ";
            long double p_hi = std::exp(-SH[i]);   // Phi_bar(TH[i]) = P(x >= TH[i])
            long double p = (i < NBINS) ? 2.0L * (p_hi - std::exp(-SH[i + 1]))
                                        : 2.0L * p_hi;
            ofs << std::setprecision(3) << (double)(std::log(p) / LN10);
        }
        ofs << "\n";
    };

    auto write_snapshot = [&](const Record& r, const char* tag) {
        SampleCount total = g_total_samples.load(std::memory_order_relaxed);
        // 纪录时间 = 样本产生时刻（r.when）；写入时间 = 本条日志落盘时刻。
        // 新纪录两者近似相同；"每分钟"快照两者相差可达一个记录间隔。
        auto now = std::chrono::system_clock::now();
        ofs << "[" << tag << "] ";
        if (!r.valid) {
            ofs << "写入时间: " << format_gmt(now)
                << " | 总抽样数: " << static_cast<unsigned long long>(total)
                << " (尚无有效样本)\n";
        } else {
            ofs << "纪录时间: " << format_gmt(r.when)
                << " | 写入时间: " << format_gmt(now)
                << " | y: " << std::setprecision(10) << (double)r.y
                << " | x: " << std::setprecision(10) << (double)r.x
                << " | 稀有度 log10 P: "
                // 双边：log10(2·Phi_bar(y)) = (ln2 - y)/ln10
                << std::setprecision(4) << (double)((LN2 - r.y) / LN10)
                << " | 该样本是第几个: "
                << static_cast<unsigned long long>(r.sample_index)
                << " | 当前总抽样数: "
                << static_cast<unsigned long long>(total)
                << "\n";
        }
        ofs << "  分级统计:\n";
        write_bins("    ");
        ofs.flush();
    };

    // ---- 文件头 ----
    ofs << "# ===== 记录开始 " << format_gmt(std::chrono::system_clock::now())
        << " (目标 |x| >= " << TARGET_SIGMA << "σ, 双边) =====\n";
    ofs << "# 阈值表 (离线 C++ 计算, 双边 ±kσ):\n";
    ofs << "#   SH[k] = -ln(Phi_bar(k)) 右尾阈值: S = -ln(u) >= SH[k]  ⟺  x >= +kσ\n";
    ofs << "#   LH[k] = -ln(Phi(k))     左尾阈值: S = -ln(u) <= LH[k]  ⟺  x <= -kσ\n";
    ofs << "#   双边坐标 y = -ln(min(u,1-u)): |x| >= kσ ⟺ y >= SH[k]（左右尾共用）\n";
    for (int i = 1; i <= NBINS; ++i) {
        ofs << "#   " << i << "σ : SH = "
            << std::setprecision(12) << (double)SH[i]
            << " | LH = " << (double)LH[i] << "\n";
    }
    ofs.flush();

    // ---- 主循环：等新纪录或超时 ----
    while (!g_stop.load(std::memory_order_relaxed)) {
        bool is_new_record = false;
        {
            std::unique_lock<std::mutex> lk(g_mutex);
            g_cv.wait_for(
                lk, std::chrono::seconds(RECORD_INTERVAL_SEC),
                [] {
                    return g_stop.load(std::memory_order_relaxed) ||
                           g_has_new_record.load(std::memory_order_relaxed);
                });
            is_new_record = g_has_new_record.exchange(false);
        }
        if (g_stop.load(std::memory_order_relaxed)) break;

        Record snapshot;
        {
            std::lock_guard<std::mutex> g(g_mutex);
            snapshot = g_best;
        }

        if (is_new_record)
            write_snapshot(snapshot, "新纪录");
        else
            write_snapshot(snapshot, "每分钟");
    }

    // ---- 停机快照 ----
    Record final_rec;
    {
        std::lock_guard<std::mutex> lk(g_mutex);
        final_rec = g_final.valid ? g_final : g_best;
    }
    if (g_reached.load()) {
        ofs << "# ===== 达标！|x| >= " << TARGET_SIGMA << "σ =====\n";
        write_snapshot(final_rec, "最终/达标");
    } else {
        ofs << "# ===== 结束 =====\n";
        write_snapshot(final_rec, "最终");
    }
    ofs << "# ===== 记录结束 " << format_gmt(std::chrono::system_clock::now())
        << " =====\n";
    ofs.close();
}

} // namespace rs

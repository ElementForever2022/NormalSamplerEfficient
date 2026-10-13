// recorder.cpp — 记录线程实现
#include "recorder.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>

#include "config.hpp"
#include "shared_state.hpp"
#include "thresholds.hpp"

namespace rs {

namespace {

// 打印 2*n（n 为 SampleCount）。总样本可达 1e23 量级，2*n 可能溢出 64 位，
// 有 __int128 时用它做十进制转换，避免回绕。
void write_doubled(std::ostream& os, SampleCount n)
{
#if defined(__SIZEOF_INT128__)
    unsigned __int128 v = (unsigned __int128)n * 2u;
    if (v == 0) { os << '0'; return; }
    char buf[40];
    int  p = 40;
    while (v) { buf[--p] = char('0' + (int)(v % 10)); v /= 10; }
    os.write(buf + p, 40 - p);
#else
    os << (n <= (~SampleCount(0) / 2) ? n * 2 : n);
#endif
}

} // namespace

void record_loop(const std::string& filename)
{
    std::ofstream ofs(filename, std::ios::out | std::ios::app);
    if (!ofs) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return;
    }

    // ---- 内部辅助：把当前全局分箱写入文件 ----
    // 紧凑格式：每类分布压缩到 2 行，bin 间用 " | " 分隔。
    // 尾部累计每个 bin 同时给出单边 (>=kσ) 与双边 (±kσ = 2×单边，正态对称) 计数。
    auto write_bins = [&](const char* indent) {
        SampleCount tail[NBINS + 1];
        SampleCount ring[NBINS + 1];
        {
            std::lock_guard<std::mutex> lk(g_bin_mutex);
            for (int i = 0; i <= NBINS; ++i) {
                tail[i] = g_tail[i];
                ring[i] = g_ring[i];
            }
        }

        const int half = (NBINS + 1) / 2;   // 每行放一半 bin

        // 尾部累计：2 行
        ofs << indent << "尾部累计 (单边 | >=kσ ; 双边 | ±kσ):";
        for (int i = 1; i <= NBINS; ++i) {
            if (i == half + 1) ofs << "\n" << indent << "                ";
            else if (i > 1)    ofs << " | ";
            ofs << std::setprecision(2) << TH[i] << "σ:"
                << static_cast<unsigned long long>(tail[i])
                << "(±";
            write_doubled(ofs, tail[i]);
            ofs << ")";
        }
        ofs << "\n";

        // 环形分布：2 行
        ofs << indent << "环形分布 (|x| 落在 [k, k+1)σ):";
        for (int i = 0; i <= NBINS; ++i) {
            if (i == half) ofs << "\n" << indent << "                ";
            else if (i > 0) ofs << " | ";
            ofs << "[" << std::setprecision(2) << TH[i] << ",";
            if (i < NBINS) ofs << std::setprecision(2) << TH[i + 1];
            else           ofs << "∞";
            ofs << ")σ:" << static_cast<unsigned long long>(ring[i]);
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
                << " | S: " << std::setprecision(10) << (double)r.S
                << " | x: " << std::setprecision(10) << (double)r.x
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
        << " (目标 sigma=" << TARGET_SIGMA << ") =====\n";
    ofs << "# 阈值表 (离线 C++ 计算, SH = -ln(Phi_bar(x))):\n";
    for (int i = 1; i <= NBINS; ++i) {
        ofs << "#   " << i << "σ : SH = "
            << std::setprecision(12) << (double)SH[i] << "\n";
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
        ofs << "# ===== 达标！sigma >= " << TARGET_SIGMA << " =====\n";
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

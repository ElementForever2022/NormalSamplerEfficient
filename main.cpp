// ============================================================================
// main.cpp — 程序入口
// ----------------------------------------------------------------------------
// 编译: 见 Makefile
// 运行: ./rare_sample_efficient [seed] [outfile] [nthreads]
// 默认: seed=42, outfile=rare_records_efficient.txt, nthreads=自动
// ============================================================================

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

#include "config.hpp"
#include "recorder.hpp"
#include "sampler.hpp"
#include "shared_state.hpp"
#include "thresholds.hpp"

int main(int argc, char** argv)
{
    using namespace rs;

    std::uint64_t seed     = 42;
    std::string   outfile  = "rare_records_efficient.txt";
    int           nthreads = NTHREADS_DEFAULT;

    if (argc >= 2) seed     = std::strtoull(argv[1], nullptr, 10);
    if (argc >= 3) outfile  = argv[2];
    if (argc >= 4) nthreads = std::atoi(argv[3]);

    init_thresholds();

    if (nthreads <= 0) {
        unsigned hc = std::thread::hardware_concurrency();
        nthreads = (hc == 0) ? 4 : (int)hc;
    }

    std::cout << "种子: " << seed << "\n"
              << "输出文件: " << outfile << "\n"
              << "目标: |x| >= " << TARGET_SIGMA << "σ 时停止（双边）\n"
              << "随机源: xoshiro256++ (双块 106 位)\n"
              << "档位: [0,1)σ ~ [" << NBINS << ",∞)σ（双边 |x|）\n"
              << "线程数: " << nthreads << "\n"
              << "阈值表 (双边 ±kσ; SH = -ln(Phi_bar(k)) 右尾, LH = -ln(Phi(k)) 左尾):\n";
    for (int i = 1; i <= NBINS; ++i) {
        std::cout << "  " << i << "σ : SH = "
                  << std::setprecision(12) << (double)SH[i]
                  << " | LH = " << (double)LH[i] << "\n";
    }
    std::cout << std::endl;

    // 每个线程用不同的种子（seed + t·golden ratio），保证流独立
    std::vector<std::thread> workers;
    for (int t = 0; t < nthreads; ++t)
        workers.emplace_back(sample_loop,
                             seed + (std::uint64_t)t * 0x9E3779B97F4A7C15ULL);

    std::thread recorder(record_loop, outfile);

    for (auto& w : workers) w.join();
    g_stop.store(true, std::memory_order_relaxed);
    g_cv.notify_all();
    recorder.join();
    return 0;
}

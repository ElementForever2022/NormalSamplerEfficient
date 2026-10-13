# rare_sample_efficient — 高效稀有事件采样器

多线程随机采样程序：从均匀分布 `u ∈ (0, 1]` 生成 `S = -ln(u)`（指数分布），
映射到标准正态尾部事件 `|x| >= kσ`，持续寻找最稀有的样本，直到找到
`>= 10σ` 的事件为止。

## 快速开始

```bash
# Linux / macOS
make
./rare_sample_efficient [seed] [outfile] [nthreads]

# Windows (MSYS2/MinGW)
mingw32-make
rare_sample_efficient.exe [seed] [outfile] [nthreads]
```

参数（全部可选）：

| 参数     | 默认值                        | 说明                        |
|----------|-------------------------------|-----------------------------|
| seed     | 42                            | 全局随机种子                |
| outfile  | rare_records_efficient.txt    | 输出文件（追加写入）        |
| nthreads | 0 = 自动（hardware_concurrency） | 采样线程数              |

调试构建（带 AddressSanitizer / UBSan）：

```bash
make debug
```

## 数学原理

对 `u ~ U(0, 1]`，`S = -ln(u)` 服从标准指数分布：`P(S > s) = e^{-s}`。

令 `S = -ln(Φ̄(x))`，其中 `Φ̄(x) = P(N(0,1) > x) = erfc(x/√2)/2` 是标准正态
单侧尾概率，则 `S` 与 `x`（σ 档位）一一对应：

```
P(|x| >= kσ) = Φ̄(k) = e^{-SH[k]}
```

因此"找大 S 样本"等价于"找极端 σ 事件"。目标 10σ 对应
`SH[10] ≈ 53.23`，即约 `e^{-53.23} ≈ 7.6e-24` 的概率——这正是程序要
长时间运行才能碰到的稀有事件。

### 106 位均匀数

用 xoshiro256++ 产生两个 53 位整数拼成 106 位：

```
U = u_hi * 2^53 + u_lo ∈ [0, 2^106)
u = (U + 1) * 2^-106 ∈ (0, 1]
S = 106·ln2 - ln(U + 1)
```

`U = 0` 时 `u` 达到最小值 `2^-106`，对应最大 `S ≈ 73.5`（约 12σ）——
这是本采样器的分辨率上限。

### 快速拒绝

`S_from_hi(u_hi)` 只用高位块估算 S 的**上界**。若上界都低于 1σ 阈值，
真实 S 必然更低，直接丢弃（约 84% 的样本走此路径），避免完整对数计算。

## 文件结构

```
rand_sample/
├── config.hpp          全局配置常量（目标σ、档位数、记录间隔…）
├── thresholds.hpp/.cpp σ 档位阈值表 TH/SH，erfc 精确计算
├── rng.hpp             xoshiro256++ 伪随机数发生器
├── scompute.hpp        核心数值：S_from_two（精算）/ S_from_hi（粗估上界）
├── shared_state.hpp/.cpp 所有线程间共享状态与锁（排查并发问题看这里）
├── sampler.hpp/.cpp    采样线程主循环
├── recorder.hpp/.cpp   记录线程：新纪录/定时/最终快照写文件
├── main.cpp            入口：解析参数、启动线程
├── Makefile            构建脚本
└── README.md
```

### 模块依赖关系

```
main
 ├── config
 ├── thresholds  ←── config
 ├── sampler     ←── config, rng, scompute, shared_state, thresholds
 ├── recorder    ←── config, shared_state, thresholds
 └── shared_state ←── config
```

## 输出格式

```
# ===== 记录开始 GMT ... (目标 sigma=10) =====
# 阈值表 (SH = -ln(Phi_bar(x))):
#   1σ : SH = 1.84102164501
#   ...
[新纪录] 时间: GMT ... | S: 11.45 | x: 4.25 | 该样本是第几个: 8468 | 当前总抽样数: 20480
  分级统计:
    尾部累计 (|x| >= k sigma):
      1σ : 38107085        ← S >= SH[1] 的累计样本数
      2σ : 5464547
      ...
    环形分布 (|x| 落在 [k, k+1) sigma):
      [0, 1)σ : 202082355  ← 恰好落在该区间的样本数
      [1, 2)σ : 32642538
      ...
[每分钟] ...               ← 每 60 秒一条（新纪录会立即写）
# ===== 达标！sigma >= 10 =====
[最终/达标] ...
```

**校验方法**：`尾部累计(kσ) / 总抽样数` 应趋近 `Φ̄(k)`。
例如 1σ ≈ 0.1587、2σ ≈ 0.0228、3σ ≈ 0.00135。
若偏差过大，优先检查 `thresholds.cpp` 的 SH 表与 `scompute.hpp`。

## 并发设计（排查并发问题必读）

所有共享状态集中在 `shared_state.hpp/.cpp`：

| 状态              | 保护方式             | 说明                          |
|-------------------|----------------------|-------------------------------|
| `g_best`          | `g_mutex`            | 当前最稀有纪录                |
| `g_best_S`        | 原子 double          | `g_best.S` 的无锁镜像，热路径快筛 |
| `g_has_new_record`| 原子 bool + `g_mutex` | 在锁内置位，防止丢失唤醒      |
| `g_stop` / `g_reached` | 原子 bool      | 停机 / 达标标志               |
| `g_total_samples` | 原子计数             | 每 1024 样本批量累加          |
| `g_tail` / `g_ring` | `g_bin_mutex`      | 全局分箱计数                  |
| 线程本地 `LocalBins` | 无锁（线程私有）  | 每 65536 样本刷入全局         |

**锁顺序约定**：`g_mutex → g_bin_mutex`，不可反向（避免死锁）。

### 线程模型

```
main
 ├─ worker[0] ── sample_loop(seed + 0·φ)   采样 + 更新纪录
 ├─ worker[1] ── sample_loop(seed + 1·φ)   （每个线程独立 RNG 流）
 ├─ ...
 └─ recorder ── record_loop(outfile)        等 g_cv，写日志
```

停止流程：任一 worker 达到 10σ → 置 `g_stop` + `g_reached` → 所有 worker
退出循环 → `main` join 后再置 `g_stop` 唤醒 recorder → recorder 写最终快照。

## 常见问题排查

**Q: 编译报错 `erfcl` / `log1pl` 未声明？**
确保 `-std=c++17` 且 `<cmath>` 已包含。MSYS2/MinGW 需链接 libm（默认自动）。

**Q: 统计偏差大（尾部计数明显偏离 Φ̄(k)）？**
1. 检查启动时打印的阈值表是否与 README 中参考值一致；
2. `S_from_hi` 必须返回**上界**（否则快速拒绝会误杀稀有样本）——见
   `scompute.hpp` 注释；
3. `S_from_two` 中 `log1pl(delta)` 不可简化为 `delta`（hi 小时误差达 0.3）。

**Q: 输出文件里中文乱码？**
Windows 控制台默认 GBK，UTF-8 源文件中文在 `Get-Content` 下可能显示异常。
用 `Get-Content -Encoding UTF8` 查看，或在代码中改用英文输出。

**Q: 程序长时间不退出？**
10σ 事件概率约 `7.6e-24`，单线程每秒约 2000 万样本，理论期望需要极长时间。
这是设计使然（分布式跑批场景）。可用 Ctrl+C 手动停止，已写入的日志不丢失。

**Q: 想改目标 σ？**
改 `config.hpp` 的 `TARGET_SIGMA` 和 `NBINS`（NBINS 需 >= 目标 σ）。
`SH` 表由 `erfc` 启动时计算，无需手动更新。

**Q: 如何验证单线程结果可复现？**
固定 seed 且 `nthreads=1` 时，样本序列完全确定。多线程时各线程 RNG 流独立，
但调度顺序不定，日志中"该样本是第几个"等字段可能不同。

## 已修复的历史 bug（供对照排查）

1. **阈值表发散**：旧版 `s_from_x` 用渐近级数 `1 - 1/x² + 3/x⁴ - ...`，
   在 x=1、2 处发散，1σ 阈值算出负值。已改用 `erfc` 精确计算。
2. **快速拒绝方向错误**：`S_from_hi(0)` 旧版返回下界，破坏"上界才能安全
   拒绝"的不变量。已改为返回 `106·ln2` 上界。
3. **S 精算误差**：`S_from_two` 旧版用 `-δ` 近似 `-ln(1+δ)` 且漏了 `+1`，
   hi 小时误差达 0.3σ 档位。已改用 `log1pl` 精确计算。
4. **丢失唤醒**：`g_has_new_record` 旧版在锁外置位，记录线程可能错过通知。
   已移入 `g_mutex` 临界区。
5. **计数丢失**：线程退出时最后不足 1024 个样本未累加到
   `g_total_samples`。已补上 `idx % TOTAL_BATCH` 的尾数刷新。

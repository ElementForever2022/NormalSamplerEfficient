#pragma once
// ============================================================================
// rng.hpp — xoshiro256++ 伪随机数发生器
// ----------------------------------------------------------------------------
// 每个采样线程持有独立实例，用 splitmix64 从 64 位种子扩展出 256 位状态。
// next() 返回 64 位；>> 11 后得到 53 位均匀整数（可精确表示为 double 尾数）。
// ============================================================================

#include <cstdint>

namespace rs {

struct Xoshiro256pp {
    std::uint64_t s[4];

    explicit Xoshiro256pp(std::uint64_t seed)
    {
        for (int i = 0; i < 4; ++i) s[i] = splitmix64(seed);
    }

    std::uint64_t next()
    {
        const std::uint64_t result = rotl(s[0] + s[3], 23) + s[0];
        const std::uint64_t t = s[1] << 17;
        s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
        s[2] ^= t;
        s[3] = rotl(s[3], 45);
        return result;
    }

private:
    static std::uint64_t rotl(std::uint64_t x, int k)
    {
        return (x << k) | (x >> (64 - k));
    }

    static std::uint64_t splitmix64(std::uint64_t& x)
    {
        std::uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
};

} // namespace rs

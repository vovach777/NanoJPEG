// NanoJPEG++ (version 3.5) -- vovach777's JPEG Decoder based on NanoJPEG
// NanoJPEG -- KeyJ's Tiny Baseline JPEG Decoder
// version 1.3.5 (2016-11-14)
// Copyright (c) 2009-2016 Martin J. Fiedler <martin.fiedler@gmx.net>
// published under the terms of the MIT license
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

// NanoJPEG++ is a C++17 rewrite of NanoJPEG with some improvements.

#pragma once
#include <utility>
#include <vector>
#include <array>
#include <stdexcept>
#include <cstdint>
#include <cassert>
#include <algorithm>
#include "profiling.hpp"

namespace nanojpeg
{
typedef enum _nj_result {
    NJ_OK = 0,        // no error, decoding successful
    NJ_NO_JPEG,       // not a JPEG file
    NJ_UNSUPPORTED,   // unsupported format
    NJ_OUT_OF_MEM,    // out of memory
    NJ_INTERNAL_ERR,  // internal error
    NJ_SYNTAX_ERROR,  // syntax error
    __NJ_FINISHED,    // used internally, will never be reported
} nj_result_t;

class nj_exception : public std::domain_error {
    public:
        nj_result_t value;
        nj_exception(nj_result_t value) : std::domain_error("jpeg decoding error!"),value(value) {}
};

inline thread_local nj_result_t nj_error{NJ_OK};
#define njThrow(e) { if (nj_error == NJ_OK ) nj_error = e; return; }

#ifndef HAS_BUILTIN_CLZ
#if defined(__has_builtin)
#if __has_builtin(__builtin_clz)
#define HAS_BUILTIN_CLZ
#endif
#else
#if defined(_MSC_VER)
#include <intrin.h>
#pragma intrinsic(_BitScanReverse)
    inline int __builtin_clz(unsigned long mask)
    {
        unsigned long where;
        // Search from LSB to MSB for first set bit.
        // Returns zero if no set bit is found.
        if (_BitScanReverse(&where, mask))
            return static_cast<int>(31 - where);
        return 32; // Undefined Behavior.
    }
#define HAS_BUILTIN_CLZ
#else
    // inline int __builtin_clz(unsigned long mask)
    // {
    //     if (mask == 0)
    //         return 32;
    //     int where = 31;
    //     while ((mask & (1UL << where)) == 0)
    //         where--;
    //     return 32 - where;
    // }
#endif
#endif
#endif

    inline int leading_ones(int peek)
    {
        #ifdef HAS_BUILTIN_CLZ
        return __builtin_clz(~(peek << 16));
        #else
        return 0;
        #endif
    }



    constexpr inline uint8_t njClip(int x)
    {
        // return (x < 0) ? 0 : ((x > 0xFF) ? 0xFF : (uint8_t)x);
        return std::clamp(x, 0, 0xFF);
    }

    constexpr inline int njDecode16(const uint8_t *pos)
    {
        return (pos[0] << 8) | pos[1];
    }


template <typename Type, int SimdSize>
union alignas(32) SimdX{
    Type data[SimdSize]{};
    static constexpr auto simd_size = SimdSize;
    using simd_type = Type;

    constexpr inline auto operator-() const {
        SimdX res;
        for (int i = 0; i < SimdSize; ++i)
            res[i] = -data[i];
        return res;
    }

    constexpr inline auto& operator /= (Type v) {
        for (int i = 0; i < SimdSize; ++i)
            data[i] /= v;
        return *this;
    }

    constexpr inline auto operator / (Type v) const {
        SimdX res = *this;
        res /= v;
        return res;
    }

    constexpr inline auto& operator *= (Type v) {
        for (int i = 0; i < SimdSize; ++i)
            data[i] *= v;
        return *this;
    }

    constexpr inline auto operator * (Type v) const {
        SimdX res = *this;
        res *= v;
        return res;
    }

    constexpr inline auto& operator *= (SimdX v) {
        for (int i = 0; i < SimdSize; ++i)
            data[i] *= v[i];
        return *this;
    }

    constexpr inline auto operator * (SimdX v) const {
        SimdX res = *this;
        res *= v;
        return res;
    }

    constexpr inline auto& operator += (SimdX v) {
        for (int i = 0; i < SimdSize; ++i)
            data[i] += v[i];
        return *this;
    }

    constexpr auto operator[](int i) const {
        return data[i];
    }

    constexpr auto& operator[](int i) {
        return data[i];
    }

    constexpr inline auto operator + (SimdX v) const {
        auto res = *this;
        res += v;
        return res;
    }

    constexpr inline auto& operator -= (SimdX v) {
        for (int i = 0; i < SimdSize; ++i)
            data[i] -= v[i];
        return *this;
    }

    constexpr inline auto operator - (SimdX v) const {
        auto res = *this;
        res -= v;
        return res;
    }

    constexpr inline auto& operator += (Type v) {
        for (int i = 0; i < SimdSize; ++i)
            data[i] += v;
        return *this;
    }

    constexpr inline auto operator + (Type v) const {
        auto res = *this;
        res += v;
        return res;
    }

    constexpr inline auto& operator -= (Type v) {
        return operator+=(-v);
    }

    constexpr inline auto operator - (Type v) const {
        return operator+(-v);
    }

    constexpr friend inline SimdX operator * (Type v, const SimdX& rhs) {
            return rhs * v;
    }


    constexpr friend inline SimdX operator + (Type v, const SimdX& rhs) {
            return rhs + v;
    }

    constexpr inline auto& operator << (int v) {
        if (v > 0) {
            *this *= Type(1 << v);
        }
        return *this;
    }

    constexpr inline auto& operator >> (int v) {
        if (v > 0) {
            *this /= Type(1 << v);
        }
        return *this;
    }

    constexpr inline auto operator << (Type v) const {
        SimdX res = *this;
        res <<= v;
        return res;
    }
    constexpr inline auto operator >> (Type v) const {
        SimdX res = *this;
        res >>= v;
        return res;
    }
    constexpr inline int size() {
        return SimdSize;
    }
    inline void load(const Type * src, int skip)
    {
        #pragma omp simd
        for (int i = 0; i < SimdSize; ++i)
        {
            data[i] = src[i*skip];
        }
    }
    inline void store(Type * store)
    {
        #pragma omp simd
        for (int i = 0; i < SimdSize; ++i)
        {
            store[i] = data[i];
        }

    }
};

using float8 = SimdX<float,8>;

inline void idct8(float8 * v, unsigned char* out, int stride)
{
    /* Even part */

    float8 tmp10 = v[0] + v[4]; /* phase 3 */
    float8 tmp11 = v[0] - v[4];

    float8 tmp13 = v[2] + v[6];                                /* phases 5-3 */
    float8 tmp12 = (v[2] - v[6]) * 1.414213562f - tmp13; /* 2*c4 */

    float8 tmp0 = tmp10 + tmp13; /* phase 2 */
    float8 tmp3 = tmp10 - tmp13;
    float8 tmp1 = tmp11 + tmp12;
    float8 tmp2 = tmp11 - tmp12;

    /* Odd part */

    float8 z13 = v[5] + v[3]; /* phase 6 */
    float8 z10 = v[5] - v[3];
    float8 z11 = v[1] + v[7];
    float8 z12 = v[1] - v[7];

    float8 tmp7 = z11 + z13;                         /* phase 5 */
    tmp11 = (z11 - z13) * 1.414213562f; /* 2*c4 */

    float8 z5 = (z10 + z12) * 1.847759065f; /* 2*c2 */
    tmp10 = z5 - z12 * 1.082392200f;       /* 2*(c2-c6) */
    tmp12 = z5 - z10 * 2.613125930f;       /* 2*(c2+c6) */

    float8 tmp6 = tmp12 - tmp7; /* phase 2 */
    float8 tmp5 = tmp11 - tmp6;
    float8 tmp4 = tmp10 - tmp5;
    if ( out ) {

        auto v = {
                    tmp0 + tmp7 + 128.5f,
                    tmp1 + tmp6 + 128.5f,
                    tmp2 + tmp5 + 128.5f,
                    tmp3 + tmp4 + 128.5f,
                    tmp3 - tmp4 + 128.5f,
                    tmp2 - tmp5 + 128.5f,
                    tmp1 - tmp6 + 128.5f,
                    tmp0 - tmp7 + 128.5f,
                    };


        #pragma omp simd
        for (int i=0; i<8; ++i)
        {
            for (int j=0; j < 8; ++j)
            {
                out[j] = njClip( std::cbegin(v)[i][j]);
            }
            out += stride;
        }
    } else {
        v[0] = tmp0 + tmp7;
        v[1] = tmp1 + tmp6;
        v[2] = tmp2 + tmp5;
        v[3] = tmp3 + tmp4;
        v[4] = tmp3 - tmp4;
        v[5] = tmp2 - tmp5;
        v[6] = tmp1 - tmp6;
        v[7] = tmp0 - tmp7;
    }
}

#ifdef __AVX__
#include <x86intrin.h>
#endif
inline void transpose(float8 *a8)
{
#ifdef __AVX__
__m256* m8 = reinterpret_cast<__m256*>(a8);
__m256 &row0=m8[0];__m256 &row1=m8[1];__m256 &row2=m8[2];__m256 &row3=m8[3];__m256 &row4=m8[4];__m256 &row5=m8[5];__m256 &row6=m8[6];__m256 &row7=m8[7];
__m256 __t0, __t1, __t2, __t3, __t4, __t5, __t6, __t7;
__m256 __tt0, __tt1, __tt2, __tt3, __tt4, __tt5, __tt6, __tt7;
__t0 = _mm256_unpacklo_ps(row0, row1);
__t1 = _mm256_unpackhi_ps(row0, row1);
__t2 = _mm256_unpacklo_ps(row2, row3);
__t3 = _mm256_unpackhi_ps(row2, row3);
__t4 = _mm256_unpacklo_ps(row4, row5);
__t5 = _mm256_unpackhi_ps(row4, row5);
__t6 = _mm256_unpacklo_ps(row6, row7);
__t7 = _mm256_unpackhi_ps(row6, row7);
__tt0 = _mm256_shuffle_ps(__t0,__t2,_MM_SHUFFLE(1,0,1,0));
__tt1 = _mm256_shuffle_ps(__t0,__t2,_MM_SHUFFLE(3,2,3,2));
__tt2 = _mm256_shuffle_ps(__t1,__t3,_MM_SHUFFLE(1,0,1,0));
__tt3 = _mm256_shuffle_ps(__t1,__t3,_MM_SHUFFLE(3,2,3,2));
__tt4 = _mm256_shuffle_ps(__t4,__t6,_MM_SHUFFLE(1,0,1,0));
__tt5 = _mm256_shuffle_ps(__t4,__t6,_MM_SHUFFLE(3,2,3,2));
__tt6 = _mm256_shuffle_ps(__t5,__t7,_MM_SHUFFLE(1,0,1,0));
__tt7 = _mm256_shuffle_ps(__t5,__t7,_MM_SHUFFLE(3,2,3,2));
row0 = _mm256_permute2f128_ps(__tt0, __tt4, 0x20);
row1 = _mm256_permute2f128_ps(__tt1, __tt5, 0x20);
row2 = _mm256_permute2f128_ps(__tt2, __tt6, 0x20);
row3 = _mm256_permute2f128_ps(__tt3, __tt7, 0x20);
row4 = _mm256_permute2f128_ps(__tt0, __tt4, 0x31);
row5 = _mm256_permute2f128_ps(__tt1, __tt5, 0x31);
row6 = _mm256_permute2f128_ps(__tt2, __tt6, 0x31);
row7 = _mm256_permute2f128_ps(__tt3, __tt7, 0x31);

#else
    for(int i=0; i<8-1; i++)
        for(int j=i+1; j<8; j++)
        {
           std::swap(a8[i][j],a8[j][i]);
        }
#endif
}


inline void idct8x8(float8 * block, unsigned char*out, int stride)
{
    idct8(block,nullptr,0);
    transpose(block);
    idct8(block,out, stride);
}

    struct BitstreamContext
    {

        uint64_t bits = 0; // stores bits read from the buffer
        const uint8_t *buffer_end = nullptr;
        const uint8_t *ptr = nullptr; // pointer to the position inside a buffer
        int bits_valid = 0;           // number of bits left in bits field

        constexpr BitstreamContext(const uint8_t *buffer, int64_t buffer_size)
        {
            buffer_end = buffer + buffer_size;
            ptr = buffer;
            bits_valid = 0;
            bits = 0;
        }

        constexpr inline uint8_t readjbyte()
        {
            if (ptr < buffer_end)
            {
                const uint8_t value = *ptr++;
                if (value != 0xff)
                {
                    return value;
                }
                if (*ptr++ == 0)
                    return 0xff;
                ptr -= 2;
                buffer_end = ptr;
            }
            return 0;
        }

        constexpr void refill()
        {
            while (bits_valid <= 56)
            {
                bits |= uint64_t(readjbyte()) << (56 - bits_valid);
                bits_valid += 8;
            }
        }

        constexpr inline void skip(int n)
        {
            bits <<= n;
            bits_valid -= n;
        }

        constexpr inline uint32_t peek()
        {
            if (bits_valid<32)
                refill();
            return bits >> 32;
        }
    };

    template <int LOCKUP_SIZE>
    struct HuffCode
    {
        struct DHTItem
        {
            uint32_t mask;
            uint16_t code;
            uint16_t symbolbit;
        };
        std::array<DHTItem, 256> abc_dht{};
        std::array<uint8_t, 17> bitseek{};
        std::array<uint16_t, 1 << LOCKUP_SIZE> fast_lockup{};

        const uint8_t *dht{nullptr};
        int max_peek{0};

        void reset() {
            std::fill(abc_dht.begin(), abc_dht.end(), DHTItem{});
            std::fill(bitseek.begin(), bitseek.end(), 0);
            std::fill(fast_lockup.begin(), fast_lockup.end(), 0);
            dht = nullptr;
            max_peek = 0;
        }

        void build_lockup()
        {
            max_peek = 0;

            // Build the fast lookup table for the Huffman decoding.
            assert(dht != nullptr);
            const uint8_t *pCount = dht;
            const uint8_t *pSymbol = dht + 16;
            int remain = (1 << LOCKUP_SIZE), spread = (1 << LOCKUP_SIZE);

            int symbols_count = 0;
            for (int codelen = 1; codelen <= LOCKUP_SIZE; ++codelen)
            {
                spread >>= 1;
                int currcnt = *pCount++;
                if (!currcnt)
                    continue;
                symbols_count += currcnt;
                remain -= currcnt << (LOCKUP_SIZE - codelen);

                for (int i = 0; i < currcnt; ++i)
                {
                    auto code = (*pSymbol++) << 8 | codelen;
                    for (int j = spread; j; --j)
                    {
                        fast_lockup.at(max_peek++) = code;
                    }
                }
            }
        }

        int build_index()
        {
            const uint8_t *symbols = dht + 16;
            const uint8_t *counts = dht;
            int huffman_code = 0;

            int dht_size = 16;
            uint32_t seek = 0;
            int ones_pos = 0;
            int ones_seek = 0;
            for (int bitlen = 1; bitlen <= 16; ++bitlen)
            {
                int count = *counts++;
                dht_size += count;
                if (count > 0)
                {
                    for (int i = 0; i < count; ++i)
                    {
                        const int code = huffman_code++;
                        auto &item = abc_dht.at(seek);
                        item.symbolbit = ((*symbols++) << 8) | bitlen; // symbol and bitlen
                        item.code = code << (16 - bitlen);
                        item.mask = (0xFFFF << (16 - bitlen)) & 0xffff;
                        const int ones = leading_ones(item.code);

                        if (ones_pos < ones)
                        {
                            while (ones_pos < ones)
                            {
                                bitseek[ones_pos++] = std::min(0xff, ones_seek); // fill up the bitseek table
                            }
                            ones_seek = seek; // set the bitseek table
                        }
                        seek += 1;
                    }
                }
                huffman_code <<= 1;
            }

            while (ones_pos < 17)
            {
                bitseek[ones_pos++] = std::min(0xff, ones_seek); // fill up the bitseek table
            }
            return dht_size;
        }

        uint16_t find_slow(const int peek) const noexcept
        {

            for (auto it = abc_dht.cbegin() + bitseek[leading_ones(peek)]; it != abc_dht.cend(); ++it)
            {
                if (!((it->code ^ peek) & it->mask))
                {
                    return it->symbolbit; // symbol and bitlen
                }
            }
            return uint16_t{0}; // error
        }

        inline uint16_t find(const int peek16) const noexcept
        {
            static_assert(LOCKUP_SIZE >= 4 && LOCKUP_SIZE <= 16, "LOCKUP_SIZE must be between 4 and 16");
            static_assert(sizeof(DHTItem) == 8);

            const int peek = peek16 >> (16 - LOCKUP_SIZE);
            if (peek >= max_peek)
            {
                return find_slow(peek16);
            }
            const auto symbolbits = fast_lockup[peek];
            return symbolbits;
        }

        inline auto njGetVLC(BitstreamContext &bs) const
        {

            uint32_t peek = bs.peek();
            const auto symbolbit = find(peek >> 16);
            uint32_t codelen = symbolbit & 0xff;
            if (codelen == 0 || codelen > 16)
            {
                if (nj_error == NJ_OK)
                    nj_error = NJ_INTERNAL_ERR;
                return std::make_pair(0,uint8_t{});
            }

            uint8_t code = uint8_t(symbolbit >> 8);
            uint32_t bits = code & 0xfU; // number of extra bits

            if (!bits)
            {
                bs.skip(codelen); // skip the bits we just read
                return std::make_pair(0, code);
            }
            peek <<= codelen; // remove the bits we just read from the peek buffer
            int bitvalue = (peek >> (32U - bits));
            bs.skip(codelen + bits);

            if (bitvalue < (1 << (bits - 1)))
                bitvalue += ((-1) << bits) + 1;
            return std::make_pair( bitvalue, code);
        }
    };

    using HuffCodeDC = HuffCode<6>;
    using HuffCodeAC = HuffCode<10>;

    struct nj_component_t
    {
        int cid{};
        int ssx{}, ssy{};
        int width{}, height{};
        int chroma_w_log2{};
        int chroma_h_log2{};
        int stride{};
        int qtsel{};
        int actabsel{}, dctabsel{};
        int dcpred{};
        int size{};
        std::vector<uint8_t> pixels{};

    };

        inline void njDecodeBlock(profiling::StopWatch &profile, BitstreamContext &bs, int &dcpred, const  HuffCodeDC&dc, const  HuffCodeAC&ac,const float *qtab, uint8_t * out, int stride)
        {
            alignas(32) union Block {
                float data[64]{};
                float8 rows[8];
                constexpr float operator [](int i) const {
                    return data[i];
                }
                constexpr float& operator [](int i) {
                    return data[i];
                }
                constexpr operator float*(){
                    return &data[0];
                }
                constexpr operator float8*(){
                    return &rows[0];
                }
                constexpr Block() {}

            } block;

            // DC coef
            dcpred += std::get<0>( dc.njGetVLC(bs) );
            //transponsed!!!
            static const uint8_t ZZ[64] =
            {0,  8,  1,  2,  9, 16, 24, 17,
            10,  3,  4, 11, 18, 25, 32, 40,
            33, 26, 19, 12,  5,  6, 13, 20,
            27, 34, 41, 48, 56, 49, 42, 35,
            28, 21, 14,  7, 15, 22, 29, 36,
            43, 50, 57, 58, 51, 44, 37, 30,
            23, 31, 38, 45, 52, 59, 60, 53,
            46, 39, 47, 54, 61, 62, 55, 63};

            block[0] = (dcpred)*qtab[0]; // DC component scaling and quantization

            int coef{0};

            do
            {
                auto [value, code] = ac.njGetVLC(bs);
                if (!code)
                    break; // EOB
                if (!(code & 0x0F) && (code != 0xF0))
                     njThrow(NJ_SYNTAX_ERROR);
                coef += (code >> 4) + 1; // RLE jumps (block fills zeros)
                if (coef > 63)
                {
                    njThrow(NJ_SYNTAX_ERROR);
                }
                block[ZZ[coef]] = value * qtab[coef]; // DCT coefficients scaling and quantization
            } while (coef < 63);

            if (nj_error != NJ_OK)
                return;

            if (coef)
            {
                idct8x8(block, out , stride);
            }
            else
            { // only DC component
                auto value = njClip(block[0] + 128.5f);

                for (int i = 0; i < 8; ++i)
                {
                    std::fill(out, out + 8, value);
                    out += stride; // next line
                }
            }
        }

    struct nj_context_t
    {
        const uint8_t *pos{};
        int size{};
        int length{};
        int width{}, height{};
        int mbwidth{}, mbheight{};
        int mbsizex{}, mbsizey{};
        int ncomp{};
        profiling::StopWatch allocations_penalty{};
        std::vector<nj_component_t> comp{};
        //int qtused{}, qtavail{};
        std::vector<std::array<float, 64>> qtab{};
        int rstinterval{};
        std::vector<HuffCodeDC> huff_DC{};
        std::vector<HuffCodeAC> huff_AC{};
        bool is_ycck{false}; // YCCK color space

        inline void allocate_pixels()
        {
            for (auto &c : comp)
            {
                c.pixels.resize(c.size);
            }
        }

        inline int get_yuv_format() {
            if (!( ncomp == 3 &&
                    comp[0].chroma_h_log2==0 && comp[0].chroma_w_log2==0 &&
                    comp[1].chroma_h_log2== comp[2].chroma_h_log2 && comp[1].chroma_w_log2== comp[2].chroma_w_log2 ))
                return 0; // not yuv format
            switch ( comp[1].chroma_w_log2 << 4 | comp[1].chroma_h_log2 )
            {
                case 0x20:  return 411;
                case 0x11:  return 420;
                case 0x10:  return 422;
                case 0x00:  return 444;
                case 0x01:  return 440;
            }
            return 0; // not yuv format

        }

        inline void Decode()
        {
            nj_error = NJ_OK;

            auto sequence = {0xff, 0xd8};
            auto newpos = std::search(pos, pos + size, sequence.begin(), sequence.end());
            if ( newpos == pos + size)
                njThrow(NJ_NO_JPEG);
            Skip(newpos - pos + 2);

            while (nj_error==NJ_OK)
            {

                if ((size < 2) || (pos[0] != 0xFF))
                {
                    njThrow(NJ_SYNTAX_ERROR);
                }
                Skip(2);
                switch (pos[-1])
                {
                // case 0xD9:
                //     return;
                case 0xC0:
                    DecodeSOF<false, true>();
                    break;
                case 0xC4:
                    DecodeDHT();
                    break;
                case 0xDB:
                    DecodeDQT();
                    break;
                case 0xDD:
                    DecodeDRI();
                    break;
                case 0xDA:
                    DecodeScan();
                    if ( nj_error == NJ_OK )
                        nj_error = __NJ_FINISHED;
                    return;
                case 0xFE:
                    SkipMarker();
                    break;
                case 0xEE:
                    DecodeLength();
                    if (length == 12)
                    {
                        is_ycck = pos[11] == 2;
                    }
                    Skip(length);

                    break;
                default:
                    if ((pos[-1] & 0xF0) == 0xE0)
                        SkipMarker();
                    else
                        njThrow(NJ_UNSUPPORTED);
                }
            }
        }

        inline void Skip(int count)
        {
            pos += count;
            size -= count;
            length -= count;
            if (size < 0)
                njThrow(NJ_SYNTAX_ERROR);
        }

        inline void DecodeLength(void)
        {
            if (size < 2)
                njThrow(NJ_SYNTAX_ERROR);
            length = njDecode16(pos);
            if (length > size)
                njThrow(NJ_SYNTAX_ERROR);
            Skip(2);
        }

        inline void SkipMarker(void)
        {
            DecodeLength();
            Skip(length);
        }

        template <bool search_sof, bool allocate_memory>
        inline void DecodeSOF(void)
        {
            if constexpr (search_sof)
            {
                auto sequence = {0xff, 0xc0};
                auto newpos = std::search(pos, pos + size, sequence.begin(), sequence.end());
                Skip(newpos - pos + 2);
                if (size < 2)
                    njThrow(NJ_SYNTAX_ERROR);
            }
            int i, ssxmax = 0, ssymax = 0;
            DecodeLength();
            if (length < 9)
                njThrow(NJ_SYNTAX_ERROR);
            if (pos[0] != 8)
                njThrow(NJ_UNSUPPORTED);
            height = njDecode16(pos + 1);
            width = njDecode16(pos + 3);
            if (!width || !height)
                njThrow(NJ_SYNTAX_ERROR);
            ncomp = pos[5];
            Skip(6);
            switch (ncomp)
            {
            case 1:
            case 3:
            case 4: // grayscale, YCbCr, CMYK
                break;
            default:
                njThrow(NJ_UNSUPPORTED);
            }
            if (length < (ncomp * 3))
                njThrow(NJ_SYNTAX_ERROR);
            comp.resize(ncomp);
            for (auto &c : comp)
            {
                c.cid = pos[0];
                if (!(c.ssx = pos[1] >> 4))
                    njThrow(NJ_SYNTAX_ERROR);
                if (c.ssx & (c.ssx - 1))
                    njThrow(NJ_UNSUPPORTED); // non-power of two
                if (!(c.ssy = pos[1] & 15))
                    njThrow(NJ_SYNTAX_ERROR);
                if (c.ssy & (c.ssy - 1))
                    njThrow(NJ_UNSUPPORTED); // non-power of two
                if ((pos[2]) & 0xFC)
                    njThrow(NJ_SYNTAX_ERROR);
                c.qtsel = pos[2];
                Skip(3);
                //qtused |= 1 << c.qtsel;
                if (c.ssx > ssxmax)
                    ssxmax = c.ssx;
                if (c.ssy > ssymax)
                    ssymax = c.ssy;
            }
            if (ncomp == 1)
            {
                comp[0].ssx = comp[0].ssy = ssxmax = ssymax = 1;
            }
            mbsizex = ssxmax << 3;
            mbsizey = ssymax << 3;
            mbwidth = (width + mbsizex - 1) / mbsizex;
            mbheight = (height + mbsizey - 1) / mbsizey;


            for (auto &c : comp)
            {
                c.width = (width * c.ssx + ssxmax - 1) / ssxmax;
                c.height = (height * c.ssy + ssymax - 1) / ssymax;
                c.stride = mbwidth * c.ssx << 3;
                if (((c.width < 3) && (c.ssx != ssxmax)) || ((c.height < 3) && (c.ssy != ssymax)))
                    njThrow(NJ_UNSUPPORTED);
                c.size = c.stride * mbheight * c.ssy << 3;
                if constexpr (allocate_memory)
                {
                    c.pixels.resize(c.size);
                }

                int srcw = c.width;
                int srch = c.height;
                while (srcw < width)
                {
                    c.chroma_w_log2 += 1;
                    srcw <<= 1;
                }
                while (srch < height)
                {
                    c.chroma_h_log2 += 1;
                    srch <<= 1;
                }
            }
            Skip(length);
        }

        inline void DecodeDHT(void)
        {
            DecodeLength();
            huff_DC.resize(2);
            huff_AC.resize(2);
            while (length >= 17)
            {
                int i = pos[0];
                if (i & 0xEC)
                    njThrow(NJ_SYNTAX_ERROR);
                if (i & 0x02)
                    njThrow(NJ_UNSUPPORTED);
                int idx = i & 1;
                bool is_AC = i & 0x10;
                Skip(1);
                const uint8_t *dht = pos;
                if (is_AC)
                {
                    auto &tab = huff_AC[idx];
                    if (tab.dht)
                        tab.reset();
                    tab.dht = dht;
                    tab.build_lockup();
                    Skip(tab.build_index());
                }
                else
                {
                    auto &tab = huff_DC[idx];
                    if (tab.dht)
                        tab.reset();
                    tab.dht = dht;
                    tab.build_lockup();
                    Skip(tab.build_index());
                }
            }
            if (length)
                njThrow(NJ_SYNTAX_ERROR);
        }

        inline void DecodeDQT(void)
        {
            DecodeLength();

            static const uint8_t AANDctScaleFactor[64] =
                {128, 178, 178, 167, 246, 167, 151, 232,
                 232, 151, 128, 209, 219, 209, 128, 101,
                 178, 197, 197, 178, 101, 69, 139, 167,
                 177, 167, 139, 69, 35, 96, 131, 151,
                 151, 131, 96, 35, 49, 91, 118, 128,
                 118, 91, 49, 46, 81, 101, 101, 81,
                 46, 42, 69, 79, 69, 42, 35, 54,
                 54, 35, 28, 37, 28, 19, 19, 10};
            qtab.resize(4);
            while (length >= 65)
            {
                int i = pos[0];
                if (i & 0xFC)
                    njThrow(NJ_SYNTAX_ERROR);
                //qtavail |= 1 << i;
                auto p = pos + 1;
                auto scale = AANDctScaleFactor;
                for (auto &t : qtab.at(i))
                    t = (/*Integral Promotion */ (*p++) * (*scale++)) * float(1. / 1024.); // scale aan quantization table
                Skip(65);
            }
            if (length)
                njThrow(NJ_SYNTAX_ERROR);
        }

        inline void DecodeDRI(void)
        {
            DecodeLength();
            if (length < 2)
                njThrow(NJ_SYNTAX_ERROR);
            rstinterval = njDecode16(pos);
            Skip(length);
        }

        inline void DecodeScan()
        {
            DecodeLength();
            if (length < (4 + 2 * ncomp))
                njThrow(NJ_SYNTAX_ERROR);
            if (pos[0] != ncomp)
                njThrow(NJ_UNSUPPORTED);
            Skip(1);
            for (auto &c : comp)
            {
                if (pos[0] != c.cid)
                    njThrow(NJ_SYNTAX_ERROR);
                if (pos[1] & 0xEE)
                    njThrow(NJ_SYNTAX_ERROR);
                c.dctabsel = (pos[1] >> 4) & 1;
                c.actabsel = (pos[1] & 1);
                Skip(2);
            }

            if (pos[0] || (pos[1] != 63) || pos[2])
                njThrow(NJ_UNSUPPORTED);
            Skip(length);
            BitstreamContext bitstream(pos, size);
            int rstcount = rstinterval, nextrst = 0;

            for (int mby = 0; mby < mbheight; ++mby)
            {
                for (int mbx = 0; mbx < mbwidth; ++mbx)
                {


                    for (auto &c : comp)
                    {

                        for (int sby = 0; sby < c.ssy; ++sby)
                            for (int sbx = 0; sbx < c.ssx; ++sbx)
                            {
                                njDecodeBlock(allocations_penalty, bitstream, c.dcpred, huff_DC[c.dctabsel],huff_AC[c.actabsel],qtab[c.qtsel].data(), c.pixels.data() + (((mby * c.ssy + sby) * c.stride + mbx * c.ssx + sbx) << 3),c.stride);
                                if (nj_error != NJ_OK)
                                    return;
                            }
                    }
                    if (rstinterval && !(--rstcount) && ( mby != mbheight-1 ))
                    {
                        bitstream.refill();
                        Skip( bitstream.ptr - pos );
                        while (size >= 2 && njDecode16(pos) == 0xffff) Skip(1); //skips fill bytes
                        if (size < 2)
                            njThrow(NJ_SYNTAX_ERROR);
                        int  code = njDecode16(pos);
                        Skip(2);
                        // if (code == 0xFFD9) // EOI marker
                        // {
                        //     while ( size >= 2 &&  njDecode16(pos) == 0xFFD9 ) Skip(2);
                        //     return; // end of image
                        // }
                        if ( (code & 0xFFF8) != 0xFFD0)
                            njThrow(NJ_SYNTAX_ERROR);

                        if ((code & 7) != nextrst)
                            njThrow(NJ_SYNTAX_ERROR);

                        bitstream = BitstreamContext(pos, size);
                        nextrst = (nextrst + 1) & 7;
                        rstcount = rstinterval;
                        for (auto &c : comp)
                            c.dcpred = 0; // reset DC prediction to 0
                    }
                }
            }
            bitstream.refill();                      // refill buffer to make sure we have enough bits to read the last marker
            Skip(bitstream.ptr - pos);
            while ( size >= 2 &&  njDecode16(pos) == 0xffff ) Skip(1); // skip padding bytes before the last marker
            while ( size >= 2 &&  njDecode16(pos) == 0xFFD9)  Skip(2); // skip repeated EOI markers
        }
    };

    struct nj_plane {
        int width{};
        int height{};
        int stride{};
        int chroma_w_log2{};
        int chroma_h_log2{};
        std::vector<uint8_t> pixels;
    };

    struct nj_result
    {
        int width{};
        int height{};
        size_t size{};
        bool is_ycck{};
        int yuv_format{}; // native form; 444 = yuv444, 422 = yuv422, 420 = yuv420, 411 = yuv411, 400 = yuv400
        double allocation_time{};
        std::vector<nj_plane> planes{};
    };

    static nj_result decode(const uint8_t *jpeg, size_t size)
    {
        nj_context_t nj{}; // clear context
        nj.pos = jpeg;
        nj.size = size;
        nj.Decode();

        if (nj_error != __NJ_FINISHED) {
            if ( nj_error == NJ_OK )
                nj_error = NJ_NO_JPEG;
            throw nj_exception(nj_error);
        }
        nj_result res{ nj.width, nj.height, size - nj.size, nj.is_ycck && nj.ncomp == 4, nj.get_yuv_format(), nj.allocations_penalty.elapsed_total, std::vector<nj_plane>(nj.comp.size()) };
        for (int i = 0; i < nj.ncomp; ++i)
        {
            auto& src   = nj.comp[i];
            auto& dst  = res.planes[i];
            dst.width = src.width;
            dst.height = src.height;
            dst.stride = src.stride;
            dst.chroma_w_log2 = src.chroma_w_log2;
            dst.chroma_h_log2 = src.chroma_h_log2;
            std::swap( dst.pixels, src.pixels);
        }
        return res;
    }

    static void decode(const uint8_t* jpeg, size_t size,  nj_result & reuse)
    {
        nj_context_t   nj{}; // clear context
        nj.pos = jpeg;
        nj.size = size;
        nj.comp.resize( reuse.planes.size() );
        for (int i=0; i < reuse.planes.size(); ++i)
        {
           nj.comp[i].pixels = std::move(reuse.planes[i].pixels);
        }
        nj.Decode();
        if (nj_error != __NJ_FINISHED) {
            if ( nj_error == NJ_OK )
                nj_error = NJ_NO_JPEG;
            throw nj_exception(nj_error);
        }
        reuse = { nj.width, nj.height, size - nj.size, nj.is_ycck && nj.ncomp == 4, nj.get_yuv_format(), nj.allocations_penalty.elapsed_total, std::vector<nj_plane>(nj.comp.size()) };
        for (int i = 0; i < nj.ncomp; ++i)
        {
            auto& src   = nj.comp[i];
            auto& dst  = reuse.planes[i];
            dst.width = src.width;
            dst.height = src.height;
            dst.stride = src.stride;
            dst.chroma_w_log2 = src.chroma_w_log2;
            dst.chroma_h_log2 = src.chroma_h_log2;
            dst.pixels = std::move( src.pixels );
        }

    }


} // namespace nj

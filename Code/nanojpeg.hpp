// NanoJPEG++ (version 3.3) -- vovach777's JPEG Decoder based on NanoJPEG
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
//#define NJ_DISABLE_LOCKUP
#ifdef HAS_XSIMD
#include <xsimd/xsimd.hpp>
namespace xs = xsimd;
#endif
#include <simde/x86/avx.h>

namespace nanojpeg
{
#define NJ_NO_JPEG "NJ_NO_JPEG"           // not a JPEG file
#define NJ_UNSUPPORTED "NJ_UNSUPPORTED"   // unsupported format
#define NJ_OUT_OF_MEM "NJ_OUT_OF_MEM"     // out of memory
#define NJ_INTERNAL_ERR "NJ_INTERNAL_ERR" // internal error
#define NJ_SYNTAX_ERROR "NJ_SYNTAX_ERROR" // syntax error

    template <typename STR>
    [[noreturn]] inline void njThrow(STR &&e)
    {
        throw std::domain_error(std::forward<STR>(e));
    }

    inline uint8_t njClip(int x)
    {
        return std::clamp(x, 0, 0xFF); // (x < 0) ? 0 : ((x > 0xFF) ? 0xFF : (uint8_t)x);
    }

    inline int njDecode16(const uint8_t *pos)
    {
        return (pos[0] << 8) | pos[1];
    }

    struct BitstreamContext
    {

        uint64_t bits = 0; // stores bits read from the buffer
        const uint8_t *buffer_end = nullptr;
        const uint8_t *ptr = nullptr; // pointer to the position inside a buffer
        int bits_valid = 0;           // number of bits left in bits field

        void set_buffer(const uint8_t *buffer, int64_t buffer_size)
        {
            buffer_end = buffer + buffer_size;
            ptr = buffer;
            bits_valid = 0;
            bits = 0;
        }

        inline uint8_t readjbyte()
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

        void refill()
        {
            while (bits_valid <= 56)
            {
                bits |= uint64_t(readjbyte()) << (56 - bits_valid);
                bits_valid += 8;
            }
        }

        inline void skip(int n)
        {
            bits <<= n;
            bits_valid -= n;
        }

        inline uint32_t peek()
        {
            if (bits_valid<32)
                refill();
            return bits >> 32;
        }
    };

    template <int LOCKUP_SIZE>
    struct HuffCode
    {
        #ifndef NJ_DISABLE_LOCKUP
        std::vector<uint16_t> fast_lockup = std::vector<uint16_t>(  1 << LOCKUP_SIZE );
        struct {
            int huffman_code{};
            const uint8_t* counts{};
            const uint8_t* symbols{};
        } directDHT{}; //Slow huffman search
        #endif

        const uint8_t *dht{nullptr};
        int max_peek{0}; //lockup last

        void reset() {
            #ifndef NJ_DISABLE_LOCKUP
            std::fill(fast_lockup.begin(), fast_lockup.end(), 0);
            directDHT={};
            #endif
            dht = nullptr;
            max_peek = 0;
        }
        /* returns DHT item size */
        int build_lockup()
        {
            max_peek = 0;

            // Build the fast lookup table for the Huffman decoding.
            assert(dht != nullptr);
            const uint8_t *pCount = dht;
            const uint8_t *pSymbol = dht + 16;
            int huffman_code = 0;
            int remain = (1 << LOCKUP_SIZE), spread = (1 << LOCKUP_SIZE);
            int symbols_count = 0;
            for (int codelen = 1; codelen <= 16; ++codelen)
            {
                spread >>= 1;
                #ifndef NJ_DISABLE_LOCKUP
                if (codelen == LOCKUP_SIZE+1)
                    directDHT = {huffman_code,pCount,pSymbol};
                #endif
                int currcnt = *pCount++;
                if (currcnt) {
                    symbols_count += currcnt;
                    huffman_code+=currcnt;
                    #ifndef NJ_DISABLE_LOCKUP
                    if ( codelen <= LOCKUP_SIZE )
                    {
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
                    #endif
                }
                huffman_code <<= 1;
            }

            return symbols_count + 16;
        }


        uint16_t find_slow(const int peek) const noexcept
        {
            #ifdef NJ_DISABLE_LOCKUP
            const uint8_t* counts = dht;
            const uint8_t* symbols = dht+16;
            int huffman_code = 0;
            int bitlen = 1;
            #else
            const uint8_t* counts = directDHT.counts;
            const uint8_t* symbols = directDHT.symbols;
            int huffman_code = directDHT.huffman_code;
            int bitlen = LOCKUP_SIZE+1;
            #endif
            for (; bitlen <= 16; ++bitlen)
            {
                int count = *counts++;
                if (count > 0) {
                    const int peekcode = peek >> (16 - bitlen);
                    // if (peekcode < huffman_code) {
                    //     std::cout << peekcode << " < " << huffman_code << std::endl;
                    // }
                    if (peekcode >= huffman_code+count) {
                        symbols+=count;
                        huffman_code+=count;
                    }
                    else
                    {
                        for (int i = 0; i < count; ++i) {
                            int symbol = *symbols++;
                            int code   = huffman_code++;
                            if ( peekcode == code) {
                                return (symbol << 8) | bitlen; // symbol and bitlen
                            }
                        }
                    }
                }
                huffman_code <<= 1;
            }
            return 0;//error value
        }

        inline uint16_t find(const int peek16) const noexcept
        {
            #ifdef NJ_DISABLE_LOCKUP
                return find_slow(peek16);
            #else
            static_assert(LOCKUP_SIZE >= 4 && LOCKUP_SIZE <= 16, "LOCKUP_SIZE must be between 4 and 16");

            const int peek = peek16 >> (16 - LOCKUP_SIZE);
            if (peek >= max_peek)
            {
                return find_slow(peek16);
            }
            const auto symbolbits = fast_lockup[peek];
            return symbolbits;
            #endif
        }

        int njGetVLC(BitstreamContext &bs, uint8_t &code) const
        {
            uint32_t peek = bs.peek();
            const auto symbolbit = find(peek >> 16);
            uint32_t codelen = symbolbit & 0xff;
            if (codelen == 0 || codelen > 16)
                njThrow(NJ_INTERNAL_ERR);

            code = uint8_t(symbolbit >> 8);
            uint32_t bits = code & 0xfU; // number of extra bits

            if (!bits)
            {
                bs.skip(codelen); // skip the bits we just read
                return 0;
            }
            peek <<= codelen; // remove the bits we just read from the peek buffer
            int bitvalue = (peek >> (32U - bits));
            bs.skip(codelen + bits);

            if (bitvalue < (1 << (bits - 1)))
                bitvalue += ((-1) << bits) + 1;
            return bitvalue;
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
        // int qtsel{};
        // int actabsel{}, dctabsel{};
        float *qtab{};
        HuffCodeDC *dc{};
        HuffCodeAC *ac{};
        int dcpred{};
        int size{};
        std::vector<uint8_t> pixels{};
        inline void njDecodeBlock(BitstreamContext &bs, uint8_t * out)
        {

                union alignas(32) U{
                float blk[64]{};
#ifdef HAS_XSIMD
                using simd_float8 = xs::batch<float, xs::avx>;
                struct {
                    simd_float8 v0,v1,v2,v3,v4,v5,v6,v7;
                };
#endif
                struct {
                    simde__m256 r0,r1,r2,r3,r4,r5,r6,r7;
                };
                U(){}
                U(float * b) {
                std::copy(b, b+64, blk);
                }

                inline float operator[](int i) const noexcept { return blk[i]; }
                inline float& operator[](int i) noexcept { return blk[i]; }

                // https://stackoverflow.com/a/25627536/1062758
                inline void transpose() {
                    simde__m256 __t0, __t1, __t2, __t3, __t4, __t5, __t6, __t7;
                    simde__m256 __tt0, __tt1, __tt2, __tt3, __tt4, __tt5, __tt6, __tt7;
                    __t0 = simde_mm256_unpacklo_ps(r0, r1);
                    __t1 = simde_mm256_unpackhi_ps(r0, r1);
                    __t2 = simde_mm256_unpacklo_ps(r2, r3);
                    __t3 = simde_mm256_unpackhi_ps(r2, r3);
                    __t4 = simde_mm256_unpacklo_ps(r4, r5);
                    __t5 = simde_mm256_unpackhi_ps(r4, r5);
                    __t6 = simde_mm256_unpacklo_ps(r6, r7);
                    __t7 = simde_mm256_unpackhi_ps(r6, r7);
                    __tt0 = simde_mm256_shuffle_ps(__t0,__t2,_MM_SHUFFLE(1,0,1,0));
                    __tt1 = simde_mm256_shuffle_ps(__t0,__t2,_MM_SHUFFLE(3,2,3,2));
                    __tt2 = simde_mm256_shuffle_ps(__t1,__t3,_MM_SHUFFLE(1,0,1,0));
                    __tt3 = simde_mm256_shuffle_ps(__t1,__t3,_MM_SHUFFLE(3,2,3,2));
                    __tt4 = simde_mm256_shuffle_ps(__t4,__t6,_MM_SHUFFLE(1,0,1,0));
                    __tt5 = simde_mm256_shuffle_ps(__t4,__t6,_MM_SHUFFLE(3,2,3,2));
                    __tt6 = simde_mm256_shuffle_ps(__t5,__t7,_MM_SHUFFLE(1,0,1,0));
                    __tt7 = simde_mm256_shuffle_ps(__t5,__t7,_MM_SHUFFLE(3,2,3,2));
                    r0 = simde_mm256_permute2f128_ps(__tt0, __tt4, 0x20);
                    r1 = simde_mm256_permute2f128_ps(__tt1, __tt5, 0x20);
                    r2 = simde_mm256_permute2f128_ps(__tt2, __tt6, 0x20);
                    r3 = simde_mm256_permute2f128_ps(__tt3, __tt7, 0x20);
                    r4 = simde_mm256_permute2f128_ps(__tt0, __tt4, 0x31);
                    r5 = simde_mm256_permute2f128_ps(__tt1, __tt5, 0x31);
                    r6 = simde_mm256_permute2f128_ps(__tt2, __tt6, 0x31);
                    r7 = simde_mm256_permute2f128_ps(__tt3, __tt7, 0x31);
                }
#ifdef HAS_XSIMD
                inline void idct8()
                {
                    /* Even part */

                    simd_float8 tmp10 = v0 + v4; /* phase 3 */
                    simd_float8 tmp11 = v0 - v4;

                    simd_float8 tmp13 = v2 + v6;                                /* phases 5-3 */
                    simd_float8 tmp12 = (v2 - v6) * 1.414213562f - tmp13; /* 2*c4 */

                    simd_float8 tmp0 = tmp10 + tmp13; /* phase 2 */
                    simd_float8 tmp3 = tmp10 - tmp13;
                    simd_float8 tmp1 = tmp11 + tmp12;
                    simd_float8 tmp2 = tmp11 - tmp12;

                    /* Odd part */

                    simd_float8 z13 = v5 + v3; /* phase 6 */
                    simd_float8 z10 = v5 - v3;
                    simd_float8 z11 = v1 + v7;
                    simd_float8 z12 = v1 - v7;

                    simd_float8 tmp7 = z11 + z13;                         /* phase 5 */
                    tmp11 = (z11 - z13) * 1.414213562f; /* 2*c4 */

                    simd_float8 z5 = (z10 + z12) * 1.847759065f; /* 2*c2 */
                    tmp10 = z5 - z12 * 1.082392200f;       /* 2*(c2-c6) */
                    tmp12 = z5 - z10 * 2.613125930f;       /* 2*(c2+c6) */

                    simd_float8 tmp6 = tmp12 - tmp7; /* phase 2 */
                    simd_float8 tmp5 = tmp11 - tmp6;
                    simd_float8 tmp4 = tmp10 - tmp5;
                    v0 = tmp0 + tmp7;
                    v1 = tmp1 + tmp6;
                    v2 = tmp2 + tmp5;
                    v3 = tmp3 + tmp4;
                    v4 = tmp3 - tmp4;
                    v5 = tmp2 - tmp5;
                    v6 = tmp1 - tmp6;
                    v7 = tmp0 - tmp7;
                }
#else
                inline void idct8() {
                    /* Even part */

                    simde__m256 tmp10 = simde_mm256_add_ps(r0, r4); /* phase 3 */
                    simde__m256 tmp11 = simde_mm256_sub_ps(r0, r4);

                    simde__m256 tmp13 = simde_mm256_add_ps(r2, r6); /* phases 5-3 */
                    simde__m256 tmp12 = simde_mm256_sub_ps(simde_mm256_mul_ps(simde_mm256_sub_ps(r2, r6), simde_mm256_set1_ps(1.414213562f)), tmp13); /* 2*c4 */

                    simde__m256 tmp0 = simde_mm256_add_ps(tmp10, tmp13); /* phase 2 */
                    simde__m256 tmp3 = simde_mm256_sub_ps(tmp10, tmp13);
                    simde__m256 tmp1 = simde_mm256_add_ps(tmp11, tmp12);
                    simde__m256 tmp2 = simde_mm256_sub_ps(tmp11, tmp12);

                    /* Odd part */

                    simde__m256 z13 = simde_mm256_add_ps(r5, r3); /* phase 6 */
                    simde__m256 z10 = simde_mm256_sub_ps(r5, r3);
                    simde__m256 z11 = simde_mm256_add_ps(r1, r7);
                    simde__m256 z12 = simde_mm256_sub_ps(r1, r7);

                    simde__m256 tmp7 = simde_mm256_add_ps(z11, z13); /* phase 5 */
                    tmp11 = simde_mm256_mul_ps(simde_mm256_sub_ps(z11, z13), simde_mm256_set1_ps(1.414213562f)); /* 2*c4 */

                    simde__m256 z5 = simde_mm256_mul_ps(simde_mm256_add_ps(z10, z12), simde_mm256_set1_ps(1.847759065f)); /* 2*c2 */
                    tmp10 = simde_mm256_sub_ps(z5, simde_mm256_mul_ps(z12, simde_mm256_set1_ps(1.082392200f))); /* 2*(c2-c6) */
                    tmp12 = simde_mm256_sub_ps(z5, simde_mm256_mul_ps(z10, simde_mm256_set1_ps(2.613125930f))); /* 2*(c2+c6) */

                    simde__m256 tmp6 = simde_mm256_sub_ps(tmp12, tmp7); /* phase 2 */
                    simde__m256 tmp5 = simde_mm256_sub_ps(tmp11, tmp6);
                    simde__m256 tmp4 = simde_mm256_sub_ps(tmp10, tmp5);

                    r0 = simde_mm256_add_ps(tmp0, tmp7);
                    r1 = simde_mm256_add_ps(tmp1, tmp6);
                    r2 = simde_mm256_add_ps(tmp2, tmp5);
                    r3 = simde_mm256_add_ps(tmp3, tmp4);
                    r4 = simde_mm256_sub_ps(tmp3, tmp4);
                    r5 = simde_mm256_sub_ps(tmp2, tmp5);
                    r6 = simde_mm256_sub_ps(tmp1, tmp6);
                    r7 = simde_mm256_sub_ps(tmp0, tmp7);
                }
#endif

                void inv2d() {
                    //transpose(); //this transpose moved out to a ZZ order
                    idct8();
                    transpose();
                    idct8();
                }
            } block{};

            uint8_t code{};
            // DC coef
            dcpred += dc->njGetVLC(bs, code);


            alignas(32) static const auto ZZ = [] {
                        std::array<uint8_t, 64> zz_t{ 0, 1, 8, 16, 9, 2, 3, 10, 17, 24, 32, 25, 18,
                                                    11, 4, 5, 12, 19, 26, 33, 40, 48, 41, 34, 27, 20, 13, 6, 7, 14, 21, 28, 35,
                                                    42, 49, 56, 57, 50, 43, 36, 29, 22, 15, 23, 30, 37, 44, 51, 58, 59, 52, 45,
                                                    38, 31, 39, 46, 53, 60, 61, 54, 47, 55, 62, 63};
                        for (int i = 0; i < 64; i++)
                        {
                            int col = zz_t[i] % 8;
                            int row = zz_t[i] / 8;
                            std::swap(col,row); // transpose
                            zz_t[i] = row * 8 + col; // transpose back to row-major order
                        }
                        return zz_t;

                        }();

            block[0] = (dcpred)*qtab[0]; // DC component scaling and quantization

            int coef{0};

            do
            {
                int value = ac->njGetVLC(bs, code);
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

            if (coef)
            {

                block.inv2d();

                for (int i = 0, j=0; i < 64; i += 8, j += stride)
                {
                    out[j]   = njClip(block[i]   + 128.5f);
                    out[j+1] = njClip(block[i+1] + 128.5f);
                    out[j+2] = njClip(block[i+2] + 128.5f);
                    out[j+3] = njClip(block[i+3] + 128.5f);
                    out[j+4] = njClip(block[i+4] + 128.5f);
                    out[j+5] = njClip(block[i+5] + 128.5f);
                    out[j+6] = njClip(block[i+6] + 128.5f);
                    out[j+7] = njClip(block[i+7] + 128.5f);
                }
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

    };

    struct nj_context_t
    {
        const uint8_t *pos{};
        int size{};
        int length{};
        int width{}, height{};
        int mbwidth{}, mbheight{};
        int mbsizex{}, mbsizey{};
        int ncomp{};
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
            if (size < 2)
                njThrow(NJ_NO_JPEG);
            if (njDecode16(pos) != 0xFFD8)
                njThrow(NJ_NO_JPEG);
            Skip(2);
            while (1)
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
                if ( qtab.size() == 4)
                    c.qtab = qtab[ pos[2] ].data();
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
                    Skip(tab.build_lockup());
                }
                else
                {
                    auto &tab = huff_DC[idx];
                    if (tab.dht)
                        tab.reset();
                    tab.dht = dht;
                    Skip(tab.build_lockup());
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
                c.dc = &huff_DC[ (pos[1] >> 4) & 1 ];
                c.ac = &huff_AC[ (pos[1] & 1) ];
                Skip(2);
            }

            if (pos[0] || (pos[1] != 63) || pos[2])
                njThrow(NJ_UNSUPPORTED);
            Skip(length);
            BitstreamContext bitstream{};
            bitstream.set_buffer(pos, size);
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
                                c.njDecodeBlock(bitstream,  c.pixels.data() + (((mby * c.ssy + sby) * c.stride + mbx * c.ssx + sbx) << 3));
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

                        bitstream.set_buffer(pos, size);
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

    struct nj_result
    {
        int width{};
        int height{};
        size_t size{};
        bool is_ycck{};
        int yuv_format{}; // native form; 444 = yuv444, 422 = yuv422, 420 = yuv420, 411 = yuv411, 400 = yuv400
        std::vector<nj_component_t> components{};
    };

    static nj_result decode(const uint8_t *jpeg, size_t size)
    {
        nj_context_t nj{}; // clear context
        nj.pos = jpeg;
        nj.size = size;
        nj.Decode();
        return {nj.width, nj.height, size - nj.size, nj.is_ycck && nj.ncomp == 4, nj.get_yuv_format(), std::move(nj.comp)}; // return components
    }

    static void decode( const uint8_t* &jpeg, size_t &size, nj_result &reuse)
    {
        nj_context_t nj{}; // clear context
        nj.comp.resize(reuse.components.size()); // reuse components buffer
        for (int i = 0; i < reuse.components.size(); ++i)
        {
            std::swap(reuse.components[i].pixels, nj.comp[i].pixels); // reuse components buffer
        }
        nj.pos = jpeg;
        nj.size = size;
        nj.Decode();
        reuse.width = nj.width;
        reuse.height = nj.height;
        reuse.size = size - nj.size;
        reuse.is_ycck = nj.is_ycck && nj.ncomp == 4;
        reuse.yuv_format = nj.get_yuv_format();
        std::swap(nj.comp, reuse.components);
        jpeg = nj.pos;
        size = nj.size;
    }

} // namespace nj

// NanoJPEG++ (version 2.0) -- vovach777's Tiny Baseline JPEG Decoder based on
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
#include <climits>
#include <cassert>
#include <algorithm>
#include <limits>
#include <memory>
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
inline int __builtin_clz(unsigned long mask)
{
    if (mask == 0)
        return 32;
    int where = 31;
    while ((mask & (1UL << where)) == 0)
        result--;
    return 32 - where;
}
#endif
#endif
#endif

inline int leading_ones(int peek)
{
    return __builtin_clz(~(peek << 16));
}

// njDecode: Decode a JPEG image.
// Decodes a memory dump of a JPEG file into internal buffers.
// Parameters:
//   jpeg = The pointer to the memory dump.
//   size = The size of the JPEG file.
// throw exeption if error occurs

// njGetWidth: Return the width (in pixels) of the most recently decoded
// image. If njDecode() failed, the result of njGetWidth() is undefined.

// njGetHeight: Return the height (in pixels) of the most recently decoded
// image. If njDecode() failed, the result of njGetHeight() is undefined.

// njIsColor: Return 1 if the most recently decoded image is a color image
// (RGB) or 0 if it is a grayscale image. If njDecode() failed, the result
// of njGetWidth() is undefined.

// njGetComponents: Returns the decoded components. No colorspace coversion.
// Returns a reference to componets. Move object if you want to keep it before
// calling njDecode again. If njDecode() failed, the result of njGetComponents() is undefined.

inline uint8_t njClip(int x)
{
    //return (x < 0) ? 0 : ((x > 0xFF) ? 0xFF : (uint8_t)x);
    return std::clamp(x, 0, 0xFF);
}

inline void idct8(float &s0, float &s1, float &s2, float &s3, float &s4, float &s5, float &s6, float &s7)
{
    // if (s1 == 0 && s2 == 0 && s3 == 0 && s4 == 0 && s5 == 0 && s6 == 0 && s7 == 0)
    // {
    //     s1 = s2 = s3 = s4 = s5 = s6 = s7 = s0;
    //     return;
    // }
    /* Even part */

    float tmp0 = s0;
    float tmp1 = s2;
    float tmp2 = s4;
    float tmp3 = s6;

    float tmp10 = tmp0 + tmp2; /* phase 3 */
    float tmp11 = tmp0 - tmp2;

    float tmp13 = tmp1 + tmp3;                                /* phases 5-3 */
    float tmp12 = (tmp1 - tmp3) * float(1.414213562) - tmp13; /* 2*c4 */

    tmp0 = tmp10 + tmp13; /* phase 2 */
    tmp3 = tmp10 - tmp13;
    tmp1 = tmp11 + tmp12;
    tmp2 = tmp11 - tmp12;

    /* Odd part */

    float tmp4 = s1;
    float tmp5 = s3;
    float tmp6 = s5;
    float tmp7 = s7;

    float z13 = tmp6 + tmp5; /* phase 6 */
    float z10 = tmp6 - tmp5;
    float z11 = tmp4 + tmp7;
    float z12 = tmp4 - tmp7;

    tmp7 = z11 + z13;                         /* phase 5 */
    tmp11 = (z11 - z13) * float(1.414213562); /* 2*c4 */

    float z5 = (z10 + z12) * float(1.847759065); /* 2*c2 */
    tmp10 = z5 - z12 * float(1.082392200);       /* 2*(c2-c6) */
    tmp12 = z5 - z10 * float(2.613125930);       /* 2*(c2+c6) */

    tmp6 = tmp12 - tmp7; /* phase 2 */
    tmp5 = tmp11 - tmp6;
    tmp4 = tmp10 - tmp5;
    s0 = (tmp0 + tmp7);
    s1 = (tmp1 + tmp6);
    s2 = (tmp2 + tmp5);
    s3 = (tmp3 + tmp4);
    s4 = (tmp3 - tmp4);
    s5 = (tmp2 - tmp5);
    s6 = (tmp1 - tmp6);
    s7 = (tmp0 - tmp7);
}

struct NanoJpeg
{

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
                buffer_end = --ptr;
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
            refill();
            return bits >> 32;
        }
    };

    template <int LOCKUP_SIZE>
    struct HuffCode
    {
        std::array<uint16_t, 1 << LOCKUP_SIZE> fast_lockup{};

        const uint8_t *dht{nullptr};
        int max_peek{0};

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

        struct DHTItem
        {
            uint32_t mask;
            uint16_t code;
            uint16_t symbolbit;
        };
        std::array<DHTItem, 256> abc_dht{};
        std::array<uint8_t, 17> bitseek{};
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

                        if (ones > ones_pos)
                        {
                            while (ones > ones_pos)
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
                bitseek[ones_pos++] = std::min(0xff,ones_seek); // fill up the bitseek table
            }
            return dht_size;
        }

        inline uint16_t find_slow(const int peek) const noexcept
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

        inline int njGetVLC(BitstreamContext &bs, uint8_t &code) const
        {
            uint32_t peek = bs.peek();
            const auto symbolbit = find(peek >> 16);
            uint32_t codelen = symbolbit & 0xff;
            if (codelen == 0 || codelen > 16)
                njThrow(NJ_INTERNAL_ERR);

            code = uint8_t(symbolbit >> 8);
            uint32_t bits = code & 0xfU; // number of extra bits

            peek <<= codelen; // remove the bits we just read from the peek buffer

            if (!bits)
            {
                bs.skip(codelen); // skip the bits we just read
                return 0;
            }

            int bitvalue = (peek >> (32U - bits));
            bs.skip(codelen + bits);

            if (bitvalue < (1 << (bits - 1)))
                bitvalue += ((-1) << bits) + 1;
            return bitvalue;
        }
    };

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
        std::vector<uint8_t> pixels{};
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
        int qtused{}, qtavail{};
        std::vector<std::vector<int>> qtab{
            {128, 178, 178, 167, 246, 167, 151, 232,
             232, 151, 128, 209, 219, 209, 128, 101,
             178, 197, 197, 178, 101, 69, 139, 167,
             177, 167, 139, 69, 35, 96, 131, 151,
             151, 131, 96, 35, 49, 91, 118, 128,
             118, 91, 49, 46, 81, 101, 101, 81,
             46, 42, 69, 79, 69, 42, 35, 54,
             54, 35, 28, 37, 28, 19, 19, 10},
            {128, 178, 178, 167, 246, 167, 151, 232,
             232, 151, 128, 209, 219, 209, 128, 101,
             178, 197, 197, 178, 101, 69, 139, 167,
             177, 167, 139, 69, 35, 96, 131, 151,
             151, 131, 96, 35, 49, 91, 118, 128,
             118, 91, 49, 46, 81, 101, 101, 81,
             46, 42, 69, 79, 69, 42, 35, 54,
             54, 35, 28, 37, 28, 19, 19, 10},
            {128, 178, 178, 167, 246, 167, 151, 232,
             232, 151, 128, 209, 219, 209, 128, 101,
             178, 197, 197, 178, 101, 69, 139, 167,
             177, 167, 139, 69, 35, 96, 131, 151,
             151, 131, 96, 35, 49, 91, 118, 128,
             118, 91, 49, 46, 81, 101, 101, 81,
             46, 42, 69, 79, 69, 42, 35, 54,
             54, 35, 28, 37, 28, 19, 19, 10},
            {128, 178, 178, 167, 246, 167, 151, 232,
             232, 151, 128, 209, 219, 209, 128, 101,
             178, 197, 197, 178, 101, 69, 139, 167,
             177, 167, 139, 69, 35, 96, 131, 151,
             151, 131, 96, 35, 49, 91, 118, 128,
             118, 91, 49, 46, 81, 101, 101, 81,
             46, 42, 69, 79, 69, 42, 35, 54,
             54, 35, 28, 37, 28, 19, 19, 10}};

        int rstinterval{};
        std::vector<HuffCode<4>>   huff_DC{};
        std::vector<HuffCode<10>>  huff_AC{};
    };
    nj_context_t nj{};

    inline void njSkip(int count)
    {

        nj.pos += count;
        nj.size -= count;
        nj.length -= count;
        if (nj.size < 0)
            njThrow(NJ_SYNTAX_ERROR);
    }

    inline int njDecode16(const uint8_t *pos) const noexcept
    {
        return (pos[0] << 8) | pos[1];
    }

    inline void njDecodeLength(void)
    {
        if (nj.size < 2)
            njThrow(NJ_SYNTAX_ERROR);
        nj.length = njDecode16(nj.pos);
        if (nj.length > nj.size)
            njThrow(NJ_SYNTAX_ERROR);
        njSkip(2);
    }

    inline void njSkipMarker(void)
    {
        njDecodeLength();
        njSkip(nj.length);
    }

    inline void njDecodeSOF(void)
    {
        int i, ssxmax = 0, ssymax = 0;
        njDecodeLength();
        if (nj.length < 9)
            njThrow(NJ_SYNTAX_ERROR);
        if (nj.pos[0] != 8)
            njThrow(NJ_UNSUPPORTED);
        nj.height = njDecode16(nj.pos + 1);
        nj.width = njDecode16(nj.pos + 3);
        if (!nj.width || !nj.height)
            njThrow(NJ_SYNTAX_ERROR);
        nj.ncomp = nj.pos[5];
        njSkip(6);
        switch (nj.ncomp)
        {
        case 1:
        case 3:
            break;
        default:
            njThrow(NJ_UNSUPPORTED);
        }
        if (nj.length < (nj.ncomp * 3))
            njThrow(NJ_SYNTAX_ERROR);
        nj.comp.resize(nj.ncomp);
        for (auto &c : nj.comp)
        {
            c.cid = nj.pos[0];
            if (!(c.ssx = nj.pos[1] >> 4))
                njThrow(NJ_SYNTAX_ERROR);
            if (c.ssx & (c.ssx - 1))
                njThrow(NJ_UNSUPPORTED); // non-power of two
            if (!(c.ssy = nj.pos[1] & 15))
                njThrow(NJ_SYNTAX_ERROR);
            if (c.ssy & (c.ssy - 1))
                njThrow(NJ_UNSUPPORTED); // non-power of two
            if ((c.qtsel = nj.pos[2]) & 0xFC)
                njThrow(NJ_SYNTAX_ERROR);
            njSkip(3);
            nj.qtused |= 1 << c.qtsel;
            if (c.ssx > ssxmax)
                ssxmax = c.ssx;
            if (c.ssy > ssymax)
                ssymax = c.ssy;
        }
        if (nj.ncomp == 1)
        {
            nj.comp[0].ssx = nj.comp[0].ssy = ssxmax = ssymax = 1;
        }
        nj.mbsizex = ssxmax << 3;
        nj.mbsizey = ssymax << 3;
        nj.mbwidth = (nj.width + nj.mbsizex - 1) / nj.mbsizex;
        nj.mbheight = (nj.height + nj.mbsizey - 1) / nj.mbsizey;

        for (auto &c : nj.comp)
        {

            c.width = (nj.width * c.ssx + ssxmax - 1) / ssxmax;
            c.height = (nj.height * c.ssy + ssymax - 1) / ssymax;
            c.stride = nj.mbwidth * c.ssx << 3;
            if (((c.width < 3) && (c.ssx != ssxmax)) || ((c.height < 3) && (c.ssy != ssymax)))
                njThrow(NJ_UNSUPPORTED);
            c.pixels.resize(c.stride * nj.mbheight * c.ssy << 3);
        }
        int dstw = nj.comp[0].width;
        int dsth = nj.comp[0].height;

        for (int i = 1; i < nj.ncomp; ++i)
        {
            auto &c = nj.comp[i];
            int srcw = c.width;
            int srch = c.height;
            while (srcw < dstw)
            {
                c.chroma_w_log2 += 1;
                srcw <<= 1;
            }
            while (srch < dsth)
            {
                c.chroma_h_log2 += 1;
                srch <<= 1;
            }
        }

        njSkip(nj.length);
    }

    inline void njDecodeDHT(void)
    {
        njDecodeLength();
        nj.huff_DC.resize(2);
        nj.huff_AC.resize(2);

        while (nj.length >= 17)
        {
            int i = nj.pos[0];
            if (i & 0xEC)
                njThrow(NJ_SYNTAX_ERROR);
            if (i & 0x02)
                njThrow(NJ_UNSUPPORTED);
            int idx = i & 1;
            bool is_AC = i & 0x10;
            njSkip(1);
            const uint8_t *dht = nj.pos;
            if (is_AC)
            {
                auto &tab = nj.huff_AC[idx];
                if ( tab.dht )
                   njThrow(NJ_SYNTAX_ERROR);
                tab.dht = dht;
                tab.build_lockup();
                njSkip(tab.build_index());
            }
            else
            {
                auto &tab = nj.huff_DC[idx];
                if ( tab.dht )
                   njThrow(NJ_SYNTAX_ERROR);
                tab.dht = dht;
                tab.build_lockup();
                njSkip(tab.build_index());
            }
        }
        if (nj.length)
            njThrow(NJ_SYNTAX_ERROR);
    }

    inline void njDecodeDQT(void)
    {
        njDecodeLength();
        while (nj.length >= 65)
        {
            int i = nj.pos[0];
            if (i & 0xFC)
                njThrow(NJ_SYNTAX_ERROR);
            nj.qtavail |= 1 << i;
            auto p = nj.pos + 1;
            for (auto &t : nj.qtab[i])
                t *= *p++;
            njSkip(65);
        }
        if (nj.length)
            njThrow(NJ_SYNTAX_ERROR);
    }

    inline void njDecodeDRI(void)
    {
        njDecodeLength();
        if (nj.length < 2)
            njThrow(NJ_SYNTAX_ERROR);
        nj.rstinterval = njDecode16(nj.pos);
        njSkip(nj.length);
    }

    template <typename DCHuff, typename ACHuff, typename QTAB>
    static inline void njDecodeBlock(DCHuff &&dc, ACHuff &&ac, QTAB &&qtab, int &dcpred, BitstreamContext &bs, int stride, uint8_t *out)
    {
        float block[64]{};

        uint8_t code{};
        // DC coef
        dcpred += dc.njGetVLC(bs, code);

        static const uint8_t ZZ[64] = {0, 1, 8, 16, 9, 2, 3, 10, 17, 24, 32, 25, 18,
                                       11, 4, 5, 12, 19, 26, 33, 40, 48, 41, 34, 27, 20, 13, 6, 7, 14, 21, 28, 35,
                                       42, 49, 56, 57, 50, 43, 36, 29, 22, 15, 23, 30, 37, 44, 51, 58, 59, 52, 45,
                                       38, 31, 39, 46, 53, 60, 61, 54, 47, 55, 62, 63};

        block[0] = (dcpred)*qtab[0] * float(1. / 1024.); // DC component scaling and quantization

        int coef{0};

        do
        {
            int value = ac.njGetVLC(bs, code);
            if (!code)
                break; // EOB
            if (!(code & 0x0F) && (code != 0xF0))
                njThrow(NJ_SYNTAX_ERROR);
            coef += (code >> 4) + 1; // RLE jumps
            if (coef > 63)
            {
                njThrow(NJ_SYNTAX_ERROR);
            }
            block[ZZ[coef]] = value * qtab[coef] * float(1. / 1024.); // DCT coefficients scaling and quantization
        } while (coef < 63);

        if (coef)
        {

            idct8(block[0], block[1], block[2], block[3], block[4], block[5], block[6], block[7]);
            idct8(block[8], block[9], block[10], block[11], block[12], block[13], block[14], block[15]);
            idct8(block[16], block[17], block[18], block[19], block[20], block[21], block[22], block[23]);
            idct8(block[24], block[25], block[26], block[27], block[28], block[29], block[30], block[31]);
            idct8(block[32], block[33], block[34], block[35], block[36], block[37], block[38], block[39]);
            idct8(block[40], block[41], block[42], block[43], block[44], block[45], block[46], block[47]);
            idct8(block[48], block[49], block[50], block[51], block[52], block[53], block[54], block[55]);
            idct8(block[56], block[57], block[58], block[59], block[60], block[61], block[62], block[63]);

            idct8(block[0], block[8], block[16], block[24], block[32], block[40], block[48], block[56]);
            idct8(block[1], block[9], block[17], block[25], block[33], block[41], block[49], block[57]);
            idct8(block[2], block[10], block[18], block[26], block[34], block[42], block[50], block[58]);
            idct8(block[3], block[11], block[19], block[27], block[35], block[43], block[51], block[59]);
            idct8(block[4], block[12], block[20], block[28], block[36], block[44], block[52], block[60]);
            idct8(block[5], block[13], block[21], block[29], block[37], block[45], block[53], block[61]);
            idct8(block[6], block[14], block[22], block[30], block[38], block[46], block[54], block[62]);
            idct8(block[7], block[15], block[23], block[31], block[39], block[47], block[55], block[63]);

            const float *blk = block;
            for (; blk != block + 64; out += stride - 8)
            {
                *out++ = njClip(*blk++ + 128.5f);
                *out++ = njClip(*blk++ + 128.5f);
                *out++ = njClip(*blk++ + 128.5f);
                *out++ = njClip(*blk++ + 128.5f);
                *out++ = njClip(*blk++ + 128.5f);
                *out++ = njClip(*blk++ + 128.5f);
                *out++ = njClip(*blk++ + 128.5f);
                *out++ = njClip(*blk++ + 128.5f);
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

    void njDecodeScan(void)
    {
        njDecodeLength();
        if (nj.length < (4 + 2 * nj.ncomp))
            njThrow(NJ_SYNTAX_ERROR);
        if (nj.pos[0] != nj.ncomp)
            njThrow(NJ_UNSUPPORTED);
        njSkip(1);
        for (auto &c : nj.comp)
        {
            if (nj.pos[0] != c.cid)
                njThrow(NJ_SYNTAX_ERROR);
            if (nj.pos[1] & 0xEE)
                njThrow(NJ_SYNTAX_ERROR);
            c.dctabsel = (nj.pos[1] >> 4) & 1;
            c.actabsel = (nj.pos[1] & 1);
            njSkip(2);
        }

        if (nj.pos[0] || (nj.pos[1] != 63) || nj.pos[2])
            njThrow(NJ_UNSUPPORTED);
        njSkip(nj.length);
        BitstreamContext bitstream{};
        bitstream.set_buffer(nj.pos, nj.size);
        int rstcount = nj.rstinterval, nextrst = 0;

        for (int mby = 0; mby < nj.mbheight; ++mby)
            for (int mbx = 0; mbx < nj.mbwidth; ++mbx)
            {

                for (auto &c : nj.comp)
                {
                    const auto &ac = nj.huff_AC[c.actabsel]; // AC table
                    const auto &dc = nj.huff_DC[c.dctabsel]; // DC table
                    const auto &qtab = nj.qtab[c.qtsel];     // quantization table

                    for (int sby = 0; sby < c.ssy; ++sby)
                        for (int sbx = 0; sbx < c.ssx; ++sbx)
                        {
                            njDecodeBlock(dc, ac, qtab, c.dcpred, bitstream, c.stride, c.pixels.data() + (((mby * c.ssy + sby) * c.stride + mbx * c.ssx + sbx) << 3));
                        }
                }
                if (nj.rstinterval && !(--rstcount))
                {

                    bitstream.refill();

                    size_t size = bitstream.ptr + 2 - nj.pos; // 2 bytes for RST marker
                    nj.pos += size;                           // skip RST marker
                    nj.size -= size;                          // skip RST marker
                    if (nj.size < 0)
                        njThrow(NJ_SYNTAX_ERROR);
                    int code = njDecode16(bitstream.ptr);
                    if (code == 0xFFD9)
                        return; // end of image
                    bitstream.set_buffer(nj.pos, nj.size);
                    if (((code & 0xFFF8) != 0xFFD0))
                    {
                        njThrow(NJ_SYNTAX_ERROR);
                    }
                    if ((code & 7) != nextrst)
                        njThrow(NJ_SYNTAX_ERROR);
                    nextrst = (nextrst + 1) & 7;
                    rstcount = nj.rstinterval;
                    for (auto &c : nj.comp)
                        c.dcpred = 0; // reset DC prediction to 0
                }
            }
        bitstream.refill();                   // refill buffer to make sure we have enough bits to read the last marker
        size_t size = bitstream.ptr - nj.pos; // 2 bytes for RST marker
        nj.pos += size;                       // skip RST marker
        nj.size -= size;                      // skip RST marker
        if (njDecode16(nj.pos) == 0xFFD9)
            nj.pos += 2; // skip EOI marker
    }

    void njDecode(const uint8_t *jpeg, const int size)
    {
        nj = nj_context_t();

        nj.pos = jpeg;
        nj.size = size & 0x7FFFFFFF;
        if (nj.size < 2)
            njThrow(NJ_NO_JPEG);
        if ((nj.pos[0] ^ 0xFF) | (nj.pos[1] ^ 0xD8))
            njThrow(NJ_NO_JPEG);
        njSkip(2);
        while (1)
        {

            if ((nj.size < 2) || (nj.pos[0] != 0xFF))
            {
                njThrow(NJ_SYNTAX_ERROR);
            }
            njSkip(2);
            switch (nj.pos[-1])
            {
            // case 0xD9:
            //     return;
            case 0xC0:
                njDecodeSOF();
                break;
            case 0xC4:
                njDecodeDHT();
                break;
            case 0xDB:
                njDecodeDQT();
                break;
            case 0xDD:
                njDecodeDRI();
                break;
            case 0xDA:
                njDecodeScan();
                return;
            case 0xFE:
                njSkipMarker();
                break;
            default:
                if ((nj.pos[-1] & 0xF0) == 0xE0)
                    njSkipMarker();
                else
                    njThrow(NJ_UNSUPPORTED);
            }
        }
    }

    int njGetWidth(void) const { return nj.width; }
    int njGetHeight(void) const { return nj.height; }
    int njIsColor(void) const { return (nj.ncomp != 1); }
    auto &njGetComponents(void) { return nj.comp; }
};

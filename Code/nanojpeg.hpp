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

///////////////////////////////////////////////////////////////////////////////
// HEADER SECTION                                                            //
// copy and pase this into nanojpeg.h if you want                            //
///////////////////////////////////////////////////////////////////////////////

#pragma once
#include <vector>
#include <stdexcept>
#include <cstdint>
#include <climits>
#include <cassert>
#include <algorithm>
#include <limits>


#define NJ_NO_JPEG "NJ_NO_JPEG"           // not a JPEG file
#define NJ_UNSUPPORTED "NJ_UNSUPPORTED"   // unsupported format
#define NJ_OUT_OF_MEM "NJ_OUT_OF_MEM"     // out of memory
#define NJ_INTERNAL_ERR "NJ_INTERNAL_ERR" // internal error
#define NJ_SYNTAX_ERROR "NJ_SYNTAX_ERROR" // syntax error

// njInit: Initialize NanoJPEG.
// For safety reasons, this should be called at least one time before using
// using any of the other NanoJPEG functions.

// njDecode: Decode a JPEG image.
// Decodes a memory dump of a JPEG file into internal buffers.
// Parameters:
//   jpeg = The pointer to the memory dump.
//   size = The size of the JPEG file.
// Return value: The error code in case of failure, or NJ_OK (zero) on success.

// njGetWidth: Return the width (in pixels) of the most recently decoded
// image. If njDecode() failed, the result of njGetWidth() is undefined.

// njGetHeight: Return the height (in pixels) of the most recently decoded
// image. If njDecode() failed, the result of njGetHeight() is undefined.

// njIsColor: Return 1 if the most recently decoded image is a color image
// (RGB) or 0 if it is a grayscale image. If njDecode() failed, the result
// of njGetWidth() is undefined.

// njGetImage: Returns the decoded image data.
// Returns a pointer to the most recently image. The memory layout it byte-
// oriented, top-down, without any padding between lines. Pixels of color
// images will be stored as three consecutive bytes for the red, green and
// blue channels. This data format is thus compatible with the PGM or PPM
// file formats and the OpenGL texture formats GL_LUMINANCE8 or GL_RGB8.
// If njDecode() failed, the result of njGetImage() is undefined.

// njGetImageSize: Returns the size (in bytes) of the image data returned
// by njGetImage(). If njDecode() failed, the result of njGetImageSize() is
// undefined.

// njDone: Uninitialize NanoJPEG.
// Resets NanoJPEG's internal state and frees all memory that has been
// allocated at run-time by NanoJPEG. It is still possible to decode another
// image after a njDone() call.


inline uint8_t njClip(const int x)
{
    return (x < 0) ? 0 : ((x > 0xFF) ? 0xFF : (uint8_t)x);
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


// inline float b1_b3 = float(1/std::cos(4*M_PI/16));
// inline float b2 = float(1/std::cos(6*M_PI/16));
// inline float b4 = float(1/std::cos(2*M_PI/16));
// inline float b5 = float(1/( std::cos(2*M_PI/16) + std::cos(6*M_PI/16)  ));

// void idct8(float s0, float s1, float s2, float s3, float s4, float s5, float s6, float s7,
//            float& d0, float& d1, float& d2, float& d3, float& d4, float& d5, float& d6, float& d7)
// {
//         const float src4 = s5;
//         const float src7 = s3;
//         const float x4  = src4 - src7;
//         const float x7  = src4 + src7;

//         const float src5 = s1;
//         const float src6 = s7;
//         const float x5  = src5 + src6;
//         const float x6  = src5 - src6;

//         const float tmp1 = b5*(x4 - x6);
//         const float stg26 = b4*(x6) - tmp1;

//         const float x24 = tmp1 - b2*(x4);

//         const float x15 = x5 - x7;
//         const float x17 = x5 + x7;

//         const float tmp2 = stg26 - x17;
//         const float tmp3 =b1_b3*(x15) - tmp2;
//         const float x44 = tmp3 + x24;

//         const float src0 = s0;
//         const float src1 = s4;
//         const float x30 = src0 + src1;
//         const float x31 = src0 - src1;

//         const float src2 = s2;
//         const float src3 = s6;
//         const float x12 = src2 - src3;
//         const float x13 = src2 + src3;

//         const float x32 = b1_b3*(x12) - x13;

//         const float x40 = x30 + x13;
//         const float x43 = x30 - x13;
//         const float x41 = x31 + x32;
//         const float x42 = x31 - x32;

//          d0 = x40 + x17;
//          d1 = x41 + tmp2;
//          d2 = x42 + tmp3;
//          d3 = x43 - x44;
//          d4 = x43 + x44;
//          d5 = x42 - tmp3;
//          d6 = x41 - tmp2;
//          d7 = x40 - x17;
// }



struct NanoJpeg
{

struct BitstreamContextBE {
    uint64_t bits = 0; // stores bits read from the buffer
    const uint8_t *marker_ptr{};
    const uint8_t *buffer = nullptr, *buffer_end = nullptr;
    const uint8_t *ptr = nullptr; // pointer to the position inside a buffer
    int bits_valid = 0; // number of bits left in bits field

/**

 * @return

 * - 0 on successful refill

 * - a negative number when bitstream end is hit

 *

 * Always succeeds when UNCHECKED_BITSTREAM_READER is enabled.

 */
inline int bits_priv_refill_64_be()
{
    bits = 0;
    while (bits_valid < 64) {
        if (ptr >= buffer_end) {
            bits_valid = 64;
            return 0;
        }
        uint8_t value = *ptr;
        if (value == 0xff) {
            if (ptr[1] != 0) { //unsafer than ptr+1 < buffer_end
                buffer_end = ptr;
                bits_valid = 64;
                return 0;
            }
            ptr += 1;
        }
        ptr += 1;
        bits |= uint64_t(value) << (56 - bits_valid);
        bits_valid += 8;
    }
    assert(bits_valid == 64);

    return 0;
}
/**

 * @return

 * - 0 on successful refill

 * - a negative number when bitstream end is hit

 *

 * Always succeeds when UNCHECKED_BITSTREAM_READER is enabled.

 */
inline int bits_priv_refill_32_be()
{
    if ( bits_valid > 32)
       return 0;
    int end_bits_valid = bits_valid + 32;
    if ( ptr >= buffer_end ) {
       bits_valid = end_bits_valid;
       return 0;
    }
    while (bits_valid < end_bits_valid ) {
        if (ptr < buffer_end) {
            uint8_t value = *ptr;
            if (value == 0xff) {
                if (ptr[1] != 0) {
                    buffer_end = ptr;
                    bits_valid = end_bits_valid;
                    return 0;
                }
                ptr += 1;
            }
            ptr += 1;
            if (value)
                bits |= uint64_t(value) << (56 - bits_valid);
        }
        bits_valid += 8;
    }

    assert(bits_valid <= 64);
    return 0;
}
/**

 * Initialize BitstreamContext.

 * @param buffer bitstream buffer, must be AV_INPUT_BUFFER_PADDING_SIZE bytes

 *        larger than the actual read bits because some optimized bitstream

 *        readers read 32 or 64 bits at once and could read over the end

 * @param bit_size the size of the buffer in bits

 * @return 0 on success, AVERROR_INVALIDDATA if the buffer_size would overflow.

 */

void set_buffer(const uint8_t *buffer,  int64_t buffer_size)
{
    this->buffer = buffer;
    buffer_end = buffer + buffer_size;
    ptr = buffer;
    bits_valid = 0;
    bits = 0;

    if (buffer_end < ptr || buffer == nullptr) {
        throw std::length_error("Invalid data found when processing input");
    }

    bits_priv_refill_64_be();
}

inline uint64_t bits_priv_val_show_be(unsigned int n)
{
    return bits >> (64 - n);
}
inline void bits_priv_skip_remaining_be(unsigned int n)
{
    bits <<= n;
    bits_valid -= n;
}
inline uint64_t bits_priv_val_get_be(unsigned int n)
{
    uint64_t ret;
    ret = bits_priv_val_show_be(n);
    bits_priv_skip_remaining_be(n);
    return ret;
}
/**

 * Return one bit from the buffer.

 */
inline unsigned int bits_read_bit_be()
{
    if (!bits_valid && bits_priv_refill_64_be() < 0)
        return 0;
    return bits_priv_val_get_be(1);
}
/**

 * Return n bits from the buffer, n has to be in the 1-32 range.

 * May be faster than bits_read() when n is not a compile-time constant and is

 * known to be non-zero;

 */
inline uint32_t bits_read_nz_be(unsigned int n)
{
    if (n > bits_valid) {
        if (bits_valid == 0)
            bits_priv_refill_64_be();
        else
        {
            bits_priv_refill_32_be();
        }
            if ( bits_valid < n )
                throw std::runtime_error("Invalid data found when processing input");

    }
    return bits_priv_val_get_be(n);
}
/**

 * Return n bits from the buffer, n has to be in the 0-32  range.

 */
inline uint32_t bits_read_be(unsigned int n)
{
    if (!n)
        return 0;
    // if (n==1) {
    //     return bits_read_bit_be();
    // }
    return bits_read_nz_be(n);
}
/**

 * Return n bits from the buffer but do not change the buffer state.

 * n has to be in the 1-32 range. May

 */
inline uint32_t bits_peek_nz_be(unsigned int n)
{
    if (n > bits_valid)
        bits_priv_refill_32_be();
    return bits_priv_val_show_be(n);
}
/**

 * Return n bits from the buffer but do not change the buffer state.

 * n has to be in the 0-32 range.

 */
inline uint32_t bits_peek_be(unsigned int n)
{
    if (!n)
        return 0;
    return bits_peek_nz_be(n);
}
/**

 * Skip n bits in the buffer.

 */
inline void bits_skip_be(unsigned int n)
{
    if (n < bits_valid)
        bits_priv_skip_remaining_be(n);
    else {
        n -= bits_valid;
        bits = 0;
        bits_valid = 0;
        if (n >= 64) {
            unsigned int skip = n / 8;
            n -= skip * 8;
            ptr += skip;
        }
        bits_priv_refill_64_be();
        if (n)
            bits_priv_skip_remaining_be(n);
    }
}

};


struct HuffTree
{

    struct HuffNode
    {
        uint16_t left{}, right{};
        uint8_t symbol{};
    };
    int freeIndex{1};
    std::vector<HuffNode> poollist{};

    inline int parse(const uint8_t *dht)
    {
        const uint8_t *pCount = dht;
        const uint8_t *pSymbol = dht + 16;
        poollist.resize(64);
        int iCurrentCode{0};
        freeIndex = 1;
        int ret_count{};
        for (int bitLen = 1; bitLen <= 16; ++bitLen)
        {
            int count = *pCount++;
            ret_count += count;
            for (int i = 0; i < count; i++)
            {
                int iSymbol = *pSymbol++;
                int iCode = iCurrentCode++;
                insert(iCode, bitLen, iSymbol);
            }
            iCurrentCode <<= 1;
        }
        poollist.resize(freeIndex);
        assert(poollist.size() < UINT16_MAX);
        return ret_count;
    }
    inline void insert(int code, int len, int symbol)
    {
        int node = 0;
        for (int i = len - 1; i >= 0; --i)
        {
            const auto v = (code >> i) & 1;
            if (v == 0)
            {
                auto &left = poollist[node].left;
                if (left == 0)
                {
                    left = (uint16_t)freeIndex++;
                }
                node = left;
            }
            else
            {
                auto &right = poollist[node].right;
                if (right == 0)
                {
                    right = (uint16_t)freeIndex++;
                }
                node = right;
            }
            if (freeIndex >= poollist.size())
                poollist.resize(poollist.size() * 2); // double the size of pool
        }
        assert(poollist[node].left == 0);  // must be a leaf node
        assert(poollist[node].right == 0); // must be a leaf node
        poollist[node].symbol = (uint8_t)symbol;
    }
    template <typename GetBit>
    inline int find(GetBit &&getbit) const
    {
        int node = 0;
        do
        {
            node = (getbit() == 0) ? poollist[node].left : poollist[node].right;
        } while (poollist[node].right || poollist[node].left);
        return poollist[node].symbol;
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
        // nj_result_t error;
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

        // std::vector<std::vector< nj_vlc_code_t>> vlctab = std::vector(4, std::vector(65536, nj_vlc_code_t{}));
        int rstinterval{};
        std::vector<HuffTree> hufftab = std::vector(4, HuffTree{});
        BitstreamContextBE bitstream{};
    };
    nj_context_t nj{};

    template <typename STR>
    inline void njThrow(STR &&e)
    {
        throw std::domain_error(std::forward<STR>(e));
    }

    inline void njSkip(int count)
    {

        nj.pos += count;
        nj.size -= count;
        nj.length -= count;
        if (nj.size < 0)
            njThrow(NJ_SYNTAX_ERROR);
    }

    inline int njDecode16(const uint8_t *pos) const
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
        for (auto & c: nj.comp)
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

        for (int i=1; i < nj.ncomp; ++i)
        {
            auto& c=nj.comp[i];
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
        while (nj.length >= 17)
        {
            int i = nj.pos[0];
            if (i & 0xEC)
                njThrow(NJ_SYNTAX_ERROR);
            if (i & 0x02)
                njThrow(NJ_UNSUPPORTED);
            i = (i | (i >> 3)) & 3; // combined DC/AC + tableid value
            njSkip(17 + nj.hufftab[i].parse(nj.pos + 1));
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

    inline int njGetVLC(HuffTree &tree, uint8_t *code)
    {
        int value = tree.find([this]()
                              { return nj.bitstream.bits_read_bit_be(); });
        if (code)
            *code = (uint8_t)value;
        int bits = value & 15;
        if (!bits)
            return 0;
        value = nj.bitstream.bits_read_nz_be(bits);
        if (value < (1 << (bits - 1)))
            value += ((-1) << bits) + 1;
        return value;
    }

    inline void njDecodeBlock(nj_component_t &c, int offset)
    {
        uint8_t code = 0;
        uint8_t * out = c.pixels.data() + offset;
        float block[64]{};
        auto &qtab = nj.qtab[c.qtsel];
        auto &dc = nj.hufftab[c.dctabsel];
        auto &ac = nj.hufftab[c.actabsel];

        // DC coef
        int dcval = njGetVLC(dc, NULL);
        c.dcpred += dcval;

        // c->dcpred += huffman_decode_dc(dc);
        auto ZZil = {0, 1, 8, 16, 9, 2, 3, 10, 17, 24, 32, 25, 18,
                     11, 4, 5, 12, 19, 26, 33, 40, 48, 41, 34, 27, 20, 13, 6, 7, 14, 21, 28, 35,
                     42, 49, 56, 57, 50, 43, 36, 29, 22, 15, 23, 30, 37, 44, 51, 58, 59, 52, 45,
                     38, 31, 39, 46, 53, 60, 61, 54, 47, 55, 62, 63};

        auto njZZopt = std::begin(ZZil);
        bool allZiros = true;

        block[0] = (c.dcpred) * qtab[0] * float(1. / 1024.); // DC component scaling and quantization

        int coef = 0;
        do
        {
            int value = njGetVLC(ac, &code);
            if (!code)
                break; // EOB
            allZiros = false;
            if (!(code & 0x0F) && (code != 0xF0))
                njThrow(NJ_SYNTAX_ERROR);
            coef += (code >> 4) + 1; // RLE jumps
            if (coef > 63)
            {
                njThrow(NJ_SYNTAX_ERROR);
            }
            block[njZZopt[coef]] = value * qtab[coef] * float(1. / 1024.); // DCT coefficients scaling and quantization
        } while (coef < 63);
        const int stride = c.stride;

        if (!allZiros)
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

            for (int i = 0; i < 64;)
            {

                out[i & 7] = njClip(block[i] + 128.5f);
                i += 1;
                if ((i & 7) == 0)
                    out += stride; // next line
            }
        }
        else
        {
            auto value = njClip(block[0] + 128.5f);
            for (int i = 0; i < 8; ++i)
            {
                std::fill(out, out + 8, value);
                out += stride; // next line
            }
        }
    }

    inline void njDecodeScan(void)
    {
        int mbx, mby, sbx, sby;
        int rstcount = nj.rstinterval, nextrst = 0;
        njDecodeLength();
        if (nj.length < (4 + 2 * nj.ncomp))
            njThrow(NJ_SYNTAX_ERROR);
        if (nj.pos[0] != nj.ncomp)
            njThrow(NJ_UNSUPPORTED);
        njSkip(1);
        for (auto & c : nj.comp)
        {
            if (nj.pos[0] != c.cid)
                njThrow(NJ_SYNTAX_ERROR);
            if (nj.pos[1] & 0xEE)
                njThrow(NJ_SYNTAX_ERROR);
            c.dctabsel = nj.pos[1] >> 4;
            c.actabsel = (nj.pos[1] & 1) | 2;
            njSkip(2);
        }

        if (nj.pos[0] || (nj.pos[1] != 63) || nj.pos[2])
            njThrow(NJ_UNSUPPORTED);
        njSkip(nj.length);
        nj.bitstream.set_buffer(nj.pos, nj.size);
        auto prev_pos = nj.pos;

        for (mbx = mby = 0;;)
        {
            for (auto & c : nj.comp)
            {
                for (sby = 0; sby < c.ssy; ++sby)
                    for (sbx = 0; sbx < c.ssx; ++sbx)
                    {
                        njDecodeBlock(c, ((mby * c.ssy + sby) * c.stride + mbx * c.ssx + sbx) << 3);
                    }
            }

            if (++mbx >= nj.mbwidth)
            {
                mbx = 0;
                if (++mby >= nj.mbheight)
                    break;
            }
            if (nj.rstinterval && !(--rstcount))
            {
                nj.bitstream.bits_priv_refill_32_be();
                if (nj.bitstream.ptr != nj.bitstream.buffer_end)
                    njThrow(NJ_SYNTAX_ERROR);
                size_t size = nj.bitstream.ptr + 2  - nj.pos; // 2 bytes for RST marker
                nj.pos += size; // skip RST marker
                nj.size -= size; // skip RST marker
                // nj.pos = nj.bitstream.ptr + 2;
                // nj.size -= nj.pos - prev_pos;
                // prev_pos = nj.pos;
                if (nj.size < 0)
                    njThrow(NJ_SYNTAX_ERROR);
                int code = njDecode16(nj.bitstream.ptr);
                nj.bitstream.set_buffer(nj.pos, nj.size);
                if (((code & 0xFFF8) != 0xFFD0))
                {
                    njThrow(NJ_SYNTAX_ERROR);
                }
                if ((code & 7) != nextrst)
                    njThrow(NJ_SYNTAX_ERROR);
                nextrst = (nextrst + 1) & 7;
                rstcount = nj.rstinterval;
                for (auto & c : nj.comp) c.dcpred = 0; // reset DC prediction to 0
            }
            // if (bitstream.ptr == bitstream.buffer_end && (bitstream.bits_valid==0)) {
            //     if( njDecode16(bitstream.ptr) != 0xFFD9) njThrow(NJ_SYNTAX_ERROR);
            //     return;
            // }
        }
        size_t size = nj.bitstream.ptr - nj.pos; // 2 bytes for RST marker
        nj.pos += size; // skip RST marker
        nj.size -= size; // skip RST marker
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
                //if ( njDecode16(nj.pos) != 0xFFD9) njThrow(NJ_SYNTAX_ERROR);
                //std::cout << std::hex << njDecode16(nj.pos) << std::endl; // 0xFFD9
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
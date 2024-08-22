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
#define NJ_INLINE inline
#define NJ_FORCE_INLINE inline

// nj_result_t: Result codes for njDecode().
// typedef enum _nj_result {
//     NJ_OK = 0,        // no error, decoding successful
//     NJ_NO_JPEG,       // not a JPEG file
//     NJ_UNSUPPORTED,   // unsupported format
//     NJ_OUT_OF_MEM,    // out of memory
//     NJ_INTERNAL_ERR,  // internal error
//     NJ_SYNTAX_ERROR,  // syntax error
//     __NJ_FINISHED,    // used internally, will never be reported
// } nj_result_t;

#define NJ_NO_JPEG  "NJ_NO_JPEG"       // not a JPEG file
#define NJ_UNSUPPORTED "NJ_UNSUPPORTED"   // unsupported format
#define NJ_OUT_OF_MEM "NJ_OUT_OF_MEM"    // out of memory
#define NJ_INTERNAL_ERR "NJ_INTERNAL_ERR"  // internal error
#define NJ_SYNTAX_ERROR "NJ_SYNTAX_ERROR"  // syntax error

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

///////////////////////////////////////////////////////////////////////////////
// CONFIGURATION SECTION                                                     //
// adjust the default settings for the NJ_ defines here                      //
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// EXAMPLE PROGRAM                                                           //
// just define _NJ_EXAMPLE_PROGRAM to compile this (requires NJ_USE_LIBC)    //
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// IMPLEMENTATION SECTION                                                    //
// you may stop reading here                                                 //
///////////////////////////////////////////////////////////////////////////////

NJ_FORCE_INLINE uint8_t njClip(const int x) {
    return (x < 0) ? 0 : ((x > 0xFF) ? 0xFF : (uint8_t) x);
}


template <bool outclip=false,typename OT>
inline void idct8( float s0, float s1, float s2, float s3, float s4, float s5, float s6, float s7,
           OT& d0, OT& d1, OT& d2, OT& d3, OT& d4, OT& d5, OT& d6, OT& d7)
{
   /* Even part */

    float tmp0 = s0;
    float tmp1 = s2;
    float tmp2 = s4;
    float tmp3 = s6;

    float tmp10 = tmp0 + tmp2;        /* phase 3 */
    float tmp11 = tmp0 - tmp2;

    float tmp13 = tmp1 + tmp3;        /* phases 5-3 */
    float tmp12 = (tmp1 - tmp3) * float(1.414213562) - tmp13; /* 2*c4 */

    tmp0 = tmp10 + tmp13;       /* phase 2 */
    tmp3 = tmp10 - tmp13;
    tmp1 = tmp11 + tmp12;
    tmp2 = tmp11 - tmp12;

    /* Odd part */

    float tmp4 = s1;
    float tmp5 = s3;
    float tmp6 = s5;
    float tmp7 = s7;

    float z13 = tmp6 + tmp5;          /* phase 6 */
    float z10 = tmp6 - tmp5;
    float z11 = tmp4 + tmp7;
    float z12 = tmp4 - tmp7;

    tmp7 = z11 + z13;           /* phase 5 */
    tmp11 = (z11 - z13) * float(1.414213562); /* 2*c4 */

    float z5 = (z10 + z12) * float(1.847759065); /* 2*c2 */
    tmp10 = z5 - z12 * float(1.082392200); /* 2*(c2-c6) */
    tmp12 = z5 - z10 * float(2.613125930); /* 2*(c2+c6) */

    tmp6 = tmp12 - tmp7;        /* phase 2 */
    tmp5 = tmp11 - tmp6;
    tmp4 = tmp10 - tmp5;
    if constexpr (outclip) {
        d0 = njClip(tmp0 + tmp7 + 128.5f);
        d7 = njClip(tmp0 - tmp7 + 128.5f);
        d1 = njClip(tmp1 + tmp6 + 128.5f);
        d6 = njClip(tmp1 - tmp6 + 128.5f);
        d2 = njClip(tmp2 + tmp5 + 128.5f);
        d5 = njClip(tmp2 - tmp5 + 128.5f);
        d3 = njClip(tmp3 + tmp4 + 128.5f);
        d4 = njClip(tmp3 - tmp4 + 128.5f);
    } else {
        d0 = (tmp0 + tmp7);
        d7 = (tmp0 - tmp7);
        d1 = (tmp1 + tmp6);
        d6 = (tmp1 - tmp6);
        d2 = (tmp2 + tmp5);
        d5 = (tmp2 - tmp5);
        d3 = (tmp3 + tmp4);
        d4 = (tmp3 - tmp4);
    }
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



struct NanoJpeg {

struct nj_vlc_code_t {
    uint8_t bits, code;
};

struct nj_component_t {
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

struct nj_context_t {
    //nj_result_t error;
    const uint8_t *pos{};
    int size{};
    int length{};
    int width{}, height{};
    int mbwidth{}, mbheight{};
    int mbsizex{}, mbsizey{};
    int ncomp{};
    std::vector< nj_component_t> comp{};
    int qtused{}, qtavail{};
    std::vector<std::vector<int>> qtab{
    {128,  178,  178,  167,  246,  167,  151,  232,
   232,  151,  128,  209,  219,  209,  128,  101,
   178,  197,  197,  178,  101,   69,  139,  167,
   177,  167,  139,   69,   35,   96,  131,  151,
   151,  131,   96,   35,   49,   91,  118,  128,
   118,   91,   49,   46,   81,  101,  101,   81,
   46,   42,   69,   79,   69,   42,   35,   54,
   54,   35,   28,   37,   28,   19,   19,   10},
    {128,  178,  178,  167,  246,  167,  151,  232,
   232,  151,  128,  209,  219,  209,  128,  101,
   178,  197,  197,  178,  101,   69,  139,  167,
   177,  167,  139,   69,   35,   96,  131,  151,
   151,  131,   96,   35,   49,   91,  118,  128,
   118,   91,   49,   46,   81,  101,  101,   81,
   46,   42,   69,   79,   69,   42,   35,   54,
   54,   35,   28,   37,   28,   19,   19,   10},
    {128,  178,  178,  167,  246,  167,  151,  232,
   232,  151,  128,  209,  219,  209,  128,  101,
   178,  197,  197,  178,  101,   69,  139,  167,
   177,  167,  139,   69,   35,   96,  131,  151,
   151,  131,   96,   35,   49,   91,  118,  128,
   118,   91,   49,   46,   81,  101,  101,   81,
   46,   42,   69,   79,   69,   42,   35,   54,
   54,   35,   28,   37,   28,   19,   19,   10},
    {128,  178,  178,  167,  246,  167,  151,  232,
   232,  151,  128,  209,  219,  209,  128,  101,
   178,  197,  197,  178,  101,   69,  139,  167,
   177,  167,  139,   69,   35,   96,  131,  151,
   151,  131,   96,   35,   49,   91,  118,  128,
   118,   91,   49,   46,   81,  101,  101,   81,
   46,   42,   69,   79,   69,   42,   35,   54,
   54,   35,   28,   37,   28,   19,   19,   10}};

    std::vector<std::vector< nj_vlc_code_t>> vlctab = std::vector(4, std::vector(65536, nj_vlc_code_t{}));
    int buf{}, bufbits{};
    int rstinterval{};
//    nj_context_t() : vlctab(4, std::vector< nj_vlc_code_t>(65536)), qtab(4,std::vector<char>(64)) {}
};
    nj_context_t nj{};





template <typename STR>
inline void njThrow(STR && e)  {
    throw std::domain_error( std::forward<STR>(e));
}
#define njCheckError()

NJ_INLINE int njShowBits(int bits) {
    uint8_t newbyte;
    if (!bits) return 0;
    while (nj.bufbits < bits) {
        if (nj.size <= 0) {
            nj.buf = (nj.buf << 8) | 0xFF;
            nj.bufbits += 8;
            continue;
        }
        newbyte = *nj.pos++;
        nj.size--;
        nj.bufbits += 8;
        nj.buf = (nj.buf << 8) | newbyte;
        if (newbyte == 0xFF) {
            if (nj.size) {
                uint8_t marker = *nj.pos++;
                nj.size--;
                switch (marker) {
                    case 0x00:
                    case 0xFF:
                        break;
                    case 0xD9: nj.size = 0; break;
                    default:
                        if ((marker & 0xF8) != 0xD0)
                            njThrow(NJ_SYNTAX_ERROR);
                        else {
                            nj.buf = (nj.buf << 8) | marker;
                            nj.bufbits += 8;
                        }
                }
            } else
               njThrow( NJ_SYNTAX_ERROR );
        }
    }
    return (nj.buf >> (nj.bufbits - bits)) & ((1 << bits) - 1);
}

NJ_INLINE void njSkipBits(int bits) {
    if (nj.bufbits < bits)
        (void) njShowBits(bits);
    nj.bufbits -= bits;
}

NJ_INLINE int njGetBits(int bits) {
    int res = njShowBits(bits);
    njSkipBits(bits);
    return res;
}

NJ_INLINE void njByteAlign(void) {
    nj.bufbits &= 0xF8;
}

NJ_INLINE void njSkip(int count) {
    nj.pos += count;
    nj.size -= count;
    nj.length -= count;
    if (nj.size < 0) njThrow(NJ_SYNTAX_ERROR);
}

NJ_INLINE unsigned short njDecode16(const uint8_t *pos) {
    return (pos[0] << 8) | pos[1];
}

NJ_INLINE void njDecodeLength(void) {
    if (nj.size < 2) njThrow(NJ_SYNTAX_ERROR);
    nj.length = njDecode16(nj.pos);
    if (nj.length > nj.size) njThrow(NJ_SYNTAX_ERROR);
    njSkip(2);
}

NJ_INLINE void njSkipMarker(void) {
    njDecodeLength();
    njSkip(nj.length);
}

NJ_INLINE void njDecodeSOF(void) {
    int i, ssxmax = 0, ssymax = 0;
    njDecodeLength();
    njCheckError();
    if (nj.length < 9) njThrow(NJ_SYNTAX_ERROR);
    if (nj.pos[0] != 8) njThrow(NJ_UNSUPPORTED);
    nj.height = njDecode16(nj.pos+1);
    nj.width = njDecode16(nj.pos+3);
    if (!nj.width || !nj.height) njThrow(NJ_SYNTAX_ERROR);
    nj.ncomp = nj.pos[5];
    nj.comp.resize( nj.ncomp );
    njSkip(6);
    switch (nj.ncomp) {
        case 1:
        case 3:
            break;
        default:
            njThrow(NJ_UNSUPPORTED);
    }
    if (nj.length < (nj.ncomp * 3)) njThrow(NJ_SYNTAX_ERROR);
    for (i = 0;  i < nj.ncomp;  ++i) {
        auto c = &nj.comp[i];
        c->cid = nj.pos[0];
        if (!(c->ssx = nj.pos[1] >> 4)) njThrow(NJ_SYNTAX_ERROR);
        if (c->ssx & (c->ssx - 1)) njThrow(NJ_UNSUPPORTED);  // non-power of two
        if (!(c->ssy = nj.pos[1] & 15)) njThrow(NJ_SYNTAX_ERROR);
        if (c->ssy & (c->ssy - 1)) njThrow(NJ_UNSUPPORTED);  // non-power of two
        if ((c->qtsel = nj.pos[2]) & 0xFC) njThrow(NJ_SYNTAX_ERROR);
        njSkip(3);
        nj.qtused |= 1 << c->qtsel;
        if (c->ssx > ssxmax) ssxmax = c->ssx;
        if (c->ssy > ssymax) ssymax = c->ssy;
    }
    if (nj.ncomp == 1) {
        auto c = &nj.comp[0];
        c->ssx = c->ssy = ssxmax = ssymax = 1;
    }
    nj.mbsizex = ssxmax << 3;
    nj.mbsizey = ssymax << 3;
    nj.mbwidth = (nj.width + nj.mbsizex - 1) / nj.mbsizex;
    nj.mbheight = (nj.height + nj.mbsizey - 1) / nj.mbsizey;
    for (i = 0;  i < nj.ncomp;  ++i) {
        auto c = &nj.comp[i];
        c->width = (nj.width * c->ssx + ssxmax - 1) / ssxmax;
        c->height = (nj.height * c->ssy + ssymax - 1) / ssymax;
        c->stride = nj.mbwidth * c->ssx << 3;
        if (((c->width < 3) && (c->ssx != ssxmax)) || ((c->height < 3) && (c->ssy != ssymax))) njThrow(NJ_UNSUPPORTED);
        //if (!(c->pixels = (uint8_t*) njAllocMem(c->stride * nj.mbheight * c->ssy << 3))) njThrow(NJ_OUT_OF_MEM);
        c->pixels.resize( c->stride * nj.mbheight * c->ssy << 3);
        if (i > 0) {
            int srcw=c->width,dstw=nj.comp[0].width,srch=c->height,dsth=nj.comp[0].height;
            while ( srcw < dstw ) {
                c->chroma_w_log2 += 1;
                srcw <<= 1;
            }
            while ( srch < dsth ) {
                c->chroma_h_log2 += 1;
                srch <<= 1;
            }

        }
    }
    njSkip(nj.length);
}

NJ_INLINE void njDecodeDHT(void) {
    int codelen, currcnt, remain, spread, i, j;
    nj_vlc_code_t *vlc;
    uint8_t counts[16];
    njDecodeLength();
    njCheckError();
    while (nj.length >= 17) {
        i = nj.pos[0];
        if (i & 0xEC) njThrow(NJ_SYNTAX_ERROR);
        if (i & 0x02) njThrow(NJ_UNSUPPORTED);
        i = (i | (i >> 3)) & 3;  // combined DC/AC + tableid value
        for (codelen = 1;  codelen <= 16;  ++codelen)
            counts[codelen - 1] = nj.pos[codelen];
        njSkip(17);
        vlc = &nj.vlctab[i][0];
        remain = spread = 65536;
        for (codelen = 1;  codelen <= 16;  ++codelen) {
            spread >>= 1;
            currcnt = counts[codelen - 1];
            if (!currcnt) continue;
            if (nj.length < currcnt) njThrow(NJ_SYNTAX_ERROR);
            remain -= currcnt << (16 - codelen);
            if (remain < 0) njThrow(NJ_SYNTAX_ERROR);
            for (i = 0;  i < currcnt;  ++i) {
                uint8_t code = nj.pos[i];
                for (j = spread;  j;  --j) {
                    vlc->bits = (uint8_t) codelen;
                    vlc->code = code;
                    ++vlc;
                }
            }
            njSkip(currcnt);
        }
        while (remain--) {
            vlc->bits = 0;
            ++vlc;
        }
    }
    if (nj.length) njThrow(NJ_SYNTAX_ERROR);
}

NJ_INLINE void njDecodeDQT(void) {
    int i;
    uint8_t *t;
    njDecodeLength();
    njCheckError();
    while (nj.length >= 65) {
        i = nj.pos[0];
        if (i & 0xFC) njThrow(NJ_SYNTAX_ERROR);
        nj.qtavail |= 1 << i;
        auto& t = nj.qtab[i];
        for (i = 0;  i < 64;  ++i)
            t[i] *= nj.pos[i + 1];
        njSkip(65);
    }
    if (nj.length) njThrow(NJ_SYNTAX_ERROR);
}

NJ_INLINE void njDecodeDRI(void) {
    njDecodeLength();
    njCheckError();
    if (nj.length < 2) njThrow(NJ_SYNTAX_ERROR);
    nj.rstinterval = njDecode16(nj.pos);
    njSkip(nj.length);
}

NJ_INLINE int njGetVLC(nj_vlc_code_t* vlc, uint8_t* code) {
    int value = njShowBits(16);
    int bits = vlc[value].bits;
    if (!bits) {
        njThrow( NJ_SYNTAX_ERROR);
         }
    njSkipBits(bits);
    value = vlc[value].code;
    if (code) *code = (uint8_t) value;
    bits = value & 15;
    if (!bits) return 0;
    value = njGetBits(bits);
    if (value < (1 << (bits - 1)))
        value += ((-1) << bits) + 1;
    return value;
}

NJ_INLINE void njDecodeBlock(nj_component_t* c, uint8_t* out) {
    uint8_t code = 0;
    int value, coef = 0;
    float block[64]{};
    auto& qtab = nj.qtab[c->qtsel];

    c->dcpred += njGetVLC(&nj.vlctab[c->dctabsel][0], NULL);
auto ZZil = { 0, 1, 8, 16, 9, 2, 3, 10, 17, 24, 32, 25, 18,
11, 4, 5, 12, 19, 26, 33, 40, 48, 41, 34, 27, 20, 13, 6, 7, 14, 21, 28, 35,
42, 49, 56, 57, 50, 43, 36, 29, 22, 15, 23, 30, 37, 44, 51, 58, 59, 52, 45,
38, 31, 39, 46, 53, 60, 61, 54, 47, 55, 62, 63 };

    auto njZZopt = std::begin(ZZil);
    bool allZiros = true;

    block[0] = (c->dcpred) * qtab[0] * float(1./1024.); // DC component scaling and quantization

    do {
        value = njGetVLC(&nj.vlctab[c->actabsel][0], &code);
        if (!code) break;  // EOB
        allZiros = false;
        if (!(code & 0x0F) && (code != 0xF0)) njThrow(NJ_SYNTAX_ERROR);
        coef += (code >> 4) + 1; //RLE jumps
        if (coef > 63) {
            njThrow(NJ_SYNTAX_ERROR);
        }
        block[ njZZopt[coef] ] = value * qtab[coef] * float(1./1024.); // DCT coefficients scaling and quantization
    } while (coef < 63);
    const int stride = c->stride;

    if (!allZiros) {



        // static const double aanscalefactor[8] = {
        //     1.0, 1.387039845, 1.306562965, 1.175875602,
        //     1.0, 0.785694958, 0.541196100, 0.275899379
        // };

        // int i = 0;
        // for (int row = 0; row < 8; row++) {
        //     for (int col = 0; col < 8; col++) {
        //         block[i] *= float(aanscalefactor[row] * aanscalefactor[col] * 0.125);
        //         i++;
        //     }
        // }

    idct8(block[0] , block[1], block[2], block[3], block[4], block[5], block[6], block[7], block[0], block[1], block[2], block[3], block[4], block[5], block[6], block[7]);
    idct8(block[8] , block[9], block[10], block[11], block[12], block[13], block[14], block[15], block[8], block[9], block[10], block[11], block[12], block[13], block[14], block[15]);
    idct8(block[16], block[17], block[18], block[19], block[20], block[21], block[22], block[23], block[16], block[17], block[18], block[19], block[20], block[21], block[22], block[23]);
    idct8(block[24], block[25], block[26], block[27], block[28], block[29], block[30], block[31], block[24], block[25], block[26], block[27], block[28], block[29], block[30], block[31]);
    idct8(block[32], block[33], block[34], block[35], block[36], block[37], block[38], block[39], block[32], block[33], block[34], block[35], block[36], block[37], block[38], block[39]);
    idct8(block[40], block[41], block[42], block[43], block[44], block[45], block[46], block[47], block[40], block[41], block[42], block[43], block[44], block[45], block[46], block[47]);
    idct8(block[48], block[49], block[50], block[51], block[52], block[53], block[54], block[55], block[48], block[49], block[50], block[51], block[52], block[53], block[54], block[55]);
    idct8(block[56], block[57], block[58], block[59], block[60], block[61], block[62], block[63], block[56], block[57], block[58], block[59], block[60], block[61], block[62], block[63]);

    idct8(block[0], block[8], block[16], block[24], block[32], block[40], block[48], block[56], block[0], block[8], block[16], block[24], block[32], block[40], block[48], block[56]);
    idct8(block[1], block[9], block[17], block[25], block[33], block[41], block[49], block[57], block[1], block[9], block[17], block[25], block[33], block[41], block[49], block[57]);
    idct8(block[2], block[10], block[18], block[26], block[34], block[42], block[50], block[58], block[2], block[10], block[18], block[26], block[34], block[42], block[50], block[58]);
    idct8(block[3], block[11], block[19], block[27], block[35], block[43], block[51], block[59], block[3], block[11], block[19], block[27], block[35], block[43], block[51], block[59]);
    idct8(block[4], block[12], block[20], block[28], block[36], block[44], block[52], block[60], block[4], block[12], block[20], block[28], block[36], block[44], block[52], block[60]);
    idct8(block[5], block[13], block[21], block[29], block[37], block[45], block[53], block[61], block[5], block[13], block[21], block[29], block[37], block[45], block[53], block[61]);
    idct8(block[6], block[14], block[22], block[30], block[38], block[46], block[54], block[62], block[6], block[14], block[22], block[30], block[38], block[46], block[54], block[62]);
    idct8(block[7], block[15], block[23], block[31], block[39], block[47], block[55], block[63], block[7], block[15], block[23], block[31], block[39], block[47], block[55], block[63]);






        for (int i = 0; i < 64;) {

            out[ i & 7 ] = njClip( block[i] + 128.5f );
            i += 1;
            if ((i & 7) == 0) out += stride; // next line

        }

    // idct8<true>(block[0], block[8], block[16], block[24], block[32], block[40], block[48], block[56],  out[0], out[0+stride], out[0+2*stride], out[0+3*stride], out[0+4*stride], out[0+5*stride], out[0+6*stride], out[0+7*stride]);
    // idct8<true>(block[1], block[9], block[17], block[25], block[33], block[41], block[49], block[57],  out[1], out[1+stride], out[1+2*stride], out[1+3*stride], out[1+4*stride], out[1+5*stride], out[1+6*stride], out[1+7*stride]);
    // idct8<true>(block[2], block[10], block[18], block[26], block[34], block[42], block[50], block[58], out[2], out[2+stride], out[2+2*stride], out[2+3*stride], out[2+4*stride], out[2+5*stride], out[2+6*stride], out[2+7*stride]);
    // idct8<true>(block[3], block[11], block[19], block[27], block[35], block[43], block[51], block[59], out[3], out[3+stride], out[3+2*stride], out[3+3*stride], out[3+4*stride], out[3+5*stride], out[3+6*stride], out[3+7*stride]);
    // idct8<true>(block[4], block[12], block[20], block[28], block[36], block[44], block[52], block[60], out[4], out[4+stride], out[4+2*stride], out[4+3*stride], out[4+4*stride], out[4+5*stride], out[4+6*stride], out[4+7*stride]);
    // idct8<true>(block[5], block[13], block[21], block[29], block[37], block[45], block[53], block[61], out[5], out[5+stride], out[5+2*stride], out[5+3*stride], out[5+4*stride], out[5+5*stride], out[5+6*stride], out[5+7*stride]);
    // idct8<true>(block[6], block[14], block[22], block[30], block[38], block[46], block[54], block[62], out[6], out[6+stride], out[6+2*stride], out[6+3*stride], out[6+4*stride], out[6+5*stride], out[6+6*stride], out[6+7*stride]);
    // idct8<true>(block[7], block[15], block[23], block[31], block[39], block[47], block[55], block[63], out[7], out[7+stride], out[7+2*stride], out[7+3*stride], out[7+4*stride], out[7+5*stride], out[7+6*stride], out[7+7*stride]);



    } else {
        auto value = njClip(block[0] + 128.5f);
        for (int i=0; i < 8;++i) {
            std::fill(out, out+8, value );
            out += stride; // next line
        }
    }



    // idct8(block[0], block[1], block[2], block[3], block[4], block[5], block[6], block[7], block[0], block[1], block[2], block[3], block[4], block[5], block[6], block[7]);
    // idct8(block[8], block[9], block[10], block[11], block[12], block[13], block[14], block[15], block[8], block[9], block[10], block[11], block[12], block[13], block[14], block[15]);
    // idct8(block[16], block[17], block[18], block[19], block[20], block[21], block[22], block[23], block[16], block[17], block[18], block[19], block[20], block[21], block[22], block[23]);
    // idct8(block[24], block[25], block[26], block[27], block[28], block[29], block[30], block[31], block[24], block[25], block[26], block[27], block[28], block[29], block[30], block[31]);
    // idct8(block[32], block[33], block[34], block[35], block[36], block[37], block[38], block[39], block[32], block[33], block[34], block[35], block[36], block[37], block[38], block[39]);
    // idct8(block[40], block[41], block[42], block[43], block[44], block[45], block[46], block[47], block[40], block[41], block[42], block[43], block[44], block[45], block[46], block[47]);
    // idct8(block[48], block[49], block[50], block[51], block[52], block[53], block[54], block[55], block[48], block[49], block[50], block[51], block[52], block[53], block[54], block[55]);
    // idct8(block[56], block[57], block[58], block[59], block[60], block[61], block[62], block[63], block[56], block[57], block[58], block[59], block[60], block[61], block[62], block[63]);

    // idct8(block[0], block[8], block[16], block[24], block[32], block[40], block[48], block[56], block[0], block[8], block[16], block[24], block[32], block[40], block[48], block[56]);
    // idct8(block[1], block[9], block[17], block[25], block[33], block[41], block[49], block[57], block[1], block[9], block[17], block[25], block[33], block[41], block[49], block[57]);
    // idct8(block[2], block[10], block[18], block[26], block[34], block[42], block[50], block[58], block[2], block[10], block[18], block[26], block[34], block[42], block[50], block[58]);
    // idct8(block[3], block[11], block[19], block[27], block[35], block[43], block[51], block[59], block[3], block[11], block[19], block[27], block[35], block[43], block[51], block[59]);
    // idct8(block[4], block[12], block[20], block[28], block[36], block[44], block[52], block[60], block[4], block[12], block[20], block[28], block[36], block[44], block[52], block[60]);
    // idct8(block[5], block[13], block[21], block[29], block[37], block[45], block[53], block[61], block[5], block[13], block[21], block[29], block[37], block[45], block[53], block[61]);
    // idct8(block[6], block[14], block[22], block[30], block[38], block[46], block[54], block[62], block[6], block[14], block[22], block[30], block[38], block[46], block[54], block[62]);
    // idct8(block[7], block[15], block[23], block[31], block[39], block[47], block[55], block[63], block[7], block[15], block[23], block[31], block[39], block[47], block[55], block[63]);



    // for (int src = 0, dst=0; src < 64;  dst += stride, src += 8) {
    //     out[dst+0] = njClip( block[src+0] / 4 + 128 );
    //     out[dst+1] = njClip( block[src+1] / 4 + 128 );
    //     out[dst+2] = njClip( block[src+2] / 4 + 128 );
    //     out[dst+3] = njClip( block[src+3] / 4 + 128 );
    //     out[dst+4] = njClip( block[src+4] / 4 + 128 );
    //     out[dst+5] = njClip( block[src+5] / 4 + 128 );
    //     out[dst+6] = njClip( block[src+6] / 4 + 128 );
    //     out[dst+7] = njClip( block[src+7] / 4 + 128 );
    // }


    // for (coef = 0;  coef < 64;  coef += 8)
    //     njRowIDCT(&block[coef]);
    // for (coef = 0;  coef < 8;  ++coef)
    //     njColIDCT(&block[coef], &out[coef], c->stride);
}

NJ_INLINE void njDecodeScan(void) {
    int i, mbx, mby, sbx, sby;
    int rstcount = nj.rstinterval, nextrst = 0;
    njDecodeLength();
    njCheckError();
    if (nj.length < (4 + 2 * nj.ncomp)) njThrow(NJ_SYNTAX_ERROR);
    if (nj.pos[0] != nj.ncomp) njThrow(NJ_UNSUPPORTED);
    njSkip(1);
    for (i = 0;  i < nj.ncomp;  ++i) {
        auto c = &nj.comp[i];
        if (nj.pos[0] != c->cid) njThrow(NJ_SYNTAX_ERROR);
        if (nj.pos[1] & 0xEE) njThrow(NJ_SYNTAX_ERROR);
        c->dctabsel = nj.pos[1] >> 4;
        c->actabsel = (nj.pos[1] & 1) | 2;
        njSkip(2);
    }
    if (nj.pos[0] || (nj.pos[1] != 63) || nj.pos[2]) njThrow(NJ_UNSUPPORTED);
    njSkip(nj.length);
    for (mbx = mby = 0;;) {
        for (i = 0;  i < nj.ncomp;  ++i) {
            auto c = &nj.comp[i];
            for (sby = 0;  sby < c->ssy;  ++sby)
                for (sbx = 0;  sbx < c->ssx;  ++sbx) {
                    njDecodeBlock(c, &c->pixels[((mby * c->ssy + sby) * c->stride + mbx * c->ssx + sbx) << 3]);
                    njCheckError();
                }
        }
        if (++mbx >= nj.mbwidth) {
            mbx = 0;
            if (++mby >= nj.mbheight) break;
        }
        if (nj.rstinterval && !(--rstcount)) {
            njByteAlign();
            i = njGetBits(16);
            if (((i & 0xFFF8) != 0xFFD0) || ((i & 7) != nextrst)) njThrow(NJ_SYNTAX_ERROR);
            nextrst = (nextrst + 1) & 7;
            rstcount = nj.rstinterval;
            for (i = 0;  i < 3;  ++i)
                nj.comp[i].dcpred = 0;
        }
    }
    //nj.error = __NJ_FINISHED;
}


void njDecode(const uint8_t* jpeg, const int size) {
    nj = nj_context_t();
    nj.pos = jpeg;
    nj.size = size & 0x7FFFFFFF;
    if (nj.size < 2) njThrow( NJ_NO_JPEG );
    if ((nj.pos[0] ^ 0xFF) | (nj.pos[1] ^ 0xD8)) njThrow(NJ_NO_JPEG);
    njSkip(2);
    while (1) {

        if ((nj.size < 2) || (nj.pos[0] != 0xFF)) {
            njThrow(NJ_SYNTAX_ERROR);
            //return;
        }
        njSkip(2);
        switch (nj.pos[-1]) {
            // case 0xD9:
            //     return;
            case 0xC0: njDecodeSOF();  break;
            case 0xC4: njDecodeDHT();  break;
            case 0xDB: njDecodeDQT();  break;
            case 0xDD: njDecodeDRI();  break;
            case 0xDA: njDecodeScan();
                       return;
            break;
            case 0xFE: njSkipMarker(); break;
            default:
                if ((nj.pos[-1] & 0xF0) == 0xE0)
                    njSkipMarker();
                else
                    njThrow( NJ_UNSUPPORTED );
        }
    }
}

int njGetWidth(void)   const         { return nj.width; }
int njGetHeight(void)  const         { return nj.height; }
int njIsColor(void)    const         { return (nj.ncomp != 1); }
auto& njGetComponents(void) { return nj.comp; }


};
#include <benchmark/benchmark.h>
#include <immintrin.h>
#include <cstdlib>
#include <cstdint>
#include <immintrin.h>
#include <algorithm>

constexpr inline uint8_t njClip(int x)
{
    return std::clamp(x, 0, 0xFF);
}

void idct8x8(const float *v, uint8_t* out, int stride)
{
    __m256 row0,row1,row2,row3,row4,row5,row6,row7;
    __m256 tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7;
    __m256 tmp10,tmp11,tmp12,tmp13;
    __m256 z5,z10,z11,z12,z13;

    row0 = _mm256_load_ps(v);
    row1 = _mm256_load_ps(v+8);
    row2 = _mm256_load_ps(v+16);
    row3 = _mm256_load_ps(v+24);
    row4 = _mm256_load_ps(v+32);
    row5 = _mm256_load_ps(v+40);
    row6 = _mm256_load_ps(v+48);
    row7 = _mm256_load_ps(v+56);
    /* Even part */

    tmp10 = _mm256_add_ps(row0, row4); /* phase 3 */
    tmp11 = _mm256_sub_ps(row0, row4);

    tmp13 = _mm256_add_ps(row2, row6); /* phases 5-3 */
    tmp12 = _mm256_sub_ps( _mm256_mul_ps( _mm256_sub_ps(row2, row6), _mm256_set1_ps( 1.414213562f )), tmp13); /* 2*c4 */

    tmp0 = _mm256_add_ps(tmp10, tmp13); /* phase 2 */
    tmp3 = _mm256_sub_ps(tmp10, tmp13);
    tmp1 = _mm256_add_ps(tmp11, tmp12);
    tmp2 = _mm256_sub_ps(tmp11, tmp12);

    /* Odd part */

    z13 = _mm256_add_ps(row5, row3); /* phase 6 */
    z10 = _mm256_sub_ps(row5, row3);
    z11 = _mm256_add_ps(row1, row7);
    z12 = _mm256_sub_ps(row1, row7);

    tmp7 = _mm256_add_ps(z11, z13);                         /* phase 5 */
    tmp11 = _mm256_mul_ps( _mm256_sub_ps(z11,z13), _mm256_set1_ps(1.414213562f)); /* 2*c4 */

    z5 = _mm256_mul_ps( _mm256_add_ps(z10, z12), _mm256_set1_ps(1.847759065f)); /* 2*c2 */
    tmp10 = _mm256_sub_ps(z5, _mm256_mul_ps(z12, _mm256_set1_ps(1.082392200f)));       /* 2*(c2-c6) */
    tmp12 = _mm256_sub_ps(z5, _mm256_mul_ps(z10, _mm256_set1_ps(2.613125930f)));       /* 2*(c2+c6) */

    tmp6 = _mm256_sub_ps(tmp12, tmp7); /* phase 2 */
    tmp5 = _mm256_sub_ps(tmp11, tmp6);
    tmp4 = _mm256_sub_ps(tmp10, tmp5);

    row0 = _mm256_add_ps(tmp0,tmp7);
    row1 = _mm256_add_ps(tmp1,tmp6);
    row2 = _mm256_add_ps(tmp2,tmp5);
    row3 = _mm256_add_ps(tmp3,tmp4);
    row4 = _mm256_sub_ps(tmp3,tmp4);
    row5 = _mm256_sub_ps(tmp2,tmp5);
    row6 = _mm256_sub_ps(tmp1,tmp6);
    row7 = _mm256_sub_ps(tmp0,tmp7);

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

    /* Even part */

    tmp10 = _mm256_add_ps(row0, row4); /* phase 3 */
    tmp11 = _mm256_sub_ps(row0, row4);

    tmp13 = _mm256_add_ps(row2, row6); /* phases 5-3 */
    tmp12 = _mm256_sub_ps( _mm256_mul_ps( _mm256_sub_ps(row2, row6), _mm256_set1_ps( 1.414213562f )), tmp13); /* 2*c4 */

    tmp0 = _mm256_add_ps(tmp10, tmp13); /* phase 2 */
    tmp3 = _mm256_sub_ps(tmp10, tmp13);
    tmp1 = _mm256_add_ps(tmp11, tmp12);
    tmp2 = _mm256_sub_ps(tmp11, tmp12);

    /* Odd part */

    z13 = _mm256_add_ps(row5, row3); /* phase 6 */
    z10 = _mm256_sub_ps(row5, row3);
    z11 = _mm256_add_ps(row1, row7);
    z12 = _mm256_sub_ps(row1, row7);

    tmp7 = _mm256_add_ps(z11, z13);                         /* phase 5 */
    tmp11 = _mm256_mul_ps( _mm256_sub_ps(z11,z13), _mm256_set1_ps(1.414213562f)); /* 2*c4 */

    z5 = _mm256_mul_ps( _mm256_add_ps(z10, z12), _mm256_set1_ps(1.847759065f)); /* 2*c2 */
    tmp10 = _mm256_sub_ps(z5, _mm256_mul_ps(z12, _mm256_set1_ps(1.082392200f)));       /* 2*(c2-c6) */
    tmp12 = _mm256_sub_ps(z5, _mm256_mul_ps(z10, _mm256_set1_ps(2.613125930f)));       /* 2*(c2+c6) */

    tmp6 = _mm256_sub_ps(tmp12, tmp7); /* phase 2 */
    tmp5 = _mm256_sub_ps(tmp11, tmp6);
    tmp4 = _mm256_sub_ps(tmp10, tmp5);

    __m256i i0 =_mm256_cvtps_epi32( _mm256_add_ps(_mm256_add_ps(tmp0,tmp7),_mm256_set1_ps(128.5f)));
    __m256i i1 =_mm256_cvtps_epi32( _mm256_add_ps(_mm256_add_ps(tmp1,tmp6),_mm256_set1_ps(128.5f)));
    __m256i i2 =_mm256_cvtps_epi32( _mm256_add_ps(_mm256_add_ps(tmp2,tmp5),_mm256_set1_ps(128.5f)));
    __m256i i3 =_mm256_cvtps_epi32( _mm256_add_ps(_mm256_add_ps(tmp3,tmp4),_mm256_set1_ps(128.5f)));
    __m256i i4 =_mm256_cvtps_epi32( _mm256_add_ps(_mm256_sub_ps(tmp3,tmp4),_mm256_set1_ps(128.5f)));
    __m256i i5 =_mm256_cvtps_epi32( _mm256_add_ps(_mm256_sub_ps(tmp2,tmp5),_mm256_set1_ps(128.5f)));
    __m256i i6 =_mm256_cvtps_epi32( _mm256_add_ps(_mm256_sub_ps(tmp1,tmp6),_mm256_set1_ps(128.5f)));
    __m256i i7 =_mm256_cvtps_epi32( _mm256_add_ps(_mm256_sub_ps(tmp0,tmp7),_mm256_set1_ps(128.5f)));
    alignas(32) int out32[64];

    _mm256_stream_si256((__m256i*)out32, i0);
    _mm256_stream_si256((__m256i*)(out32+8), i1);
    _mm256_stream_si256((__m256i*)(out32+16), i2);
    _mm256_stream_si256((__m256i*)(out32+24), i3);
    _mm256_stream_si256((__m256i*)(out32+32), i4);
    _mm256_stream_si256((__m256i*)(out32+40), i5);
    _mm256_stream_si256((__m256i*)(out32+48), i6);
    _mm256_stream_si256((__m256i*)(out32+56), i7);

      for (int i = 0, j=0; i < 64; i += 8, j += stride)
      {
          out[j]   = njClip(out32[i]   );
          out[j+1] = njClip(out32[i+1] );
          out[j+2] = njClip(out32[i+2] );
          out[j+3] = njClip(out32[i+3] );
          out[j+4] = njClip(out32[i+4] );
          out[j+5] = njClip(out32[i+5] );
          out[j+6] = njClip(out32[i+6] );
          out[j+7] = njClip(out32[i+7] );
      }
  }



template <unsigned i=0, int fract>
inline __m128i _mm_mul_epi16_const(__m128i t0)  {
    __m128i t1 = _mm_mulhrs_epi16(t0, _mm_set1_epi16(fract));
    if constexpr (i != 0) {

        if constexpr (i == 2) {
            t0 = _mm_adds_epi16( t0, t0 );
        }
        if constexpr (  fract < 0) {
            t1 = _mm_subs_epi16( t1, t0 );

        } else {
            t1 = _mm_adds_epi16( t1, t0 );
        }
    }
    return t1;
}
constexpr int fix_pass = 4;
    constexpr int fix_qtab = (10-fix_pass);


void idct8x8_epi16(const int16_t *v, uint8_t * out, int stride)
{
    __m128i row0,row1,row2,row3,row4,row5,row6,row7;
    __m128i tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7;
    __m128i tmp10,tmp11,tmp12,tmp13;
    __m128i z5,z10,z11,z12,z13;

    row0 = _mm_load_si128( (const __m128i*)(v) );
    row1 = _mm_load_si128( (const __m128i*)(v+8));
    row2 = _mm_load_si128( (const __m128i*)(v+16));
    row3 = _mm_load_si128( (const __m128i*)(v+24));
    row4 = _mm_load_si128( (const __m128i*)(v+32));
    row5 = _mm_load_si128( (const __m128i*)(v+40));
    row6 = _mm_load_si128( (const __m128i*)(v+48));
    row7 = _mm_load_si128( (const __m128i*)(v+56));
    /* Even part */
    tmp10 = _mm_adds_epi16(row0, row4); /* phase 3 */
    tmp11 = _mm_subs_epi16(row0, row4);

    tmp13 = _mm_adds_epi16(row2, row6); /* phases 5-3 */
    tmp12 = _mm_mul_epi16_const<1,int(0.414213562*0x8000)>(_mm_subs_epi16(row2, row6));
    tmp12 = _mm_subs_epi16(tmp12,tmp13);

    tmp0 = _mm_adds_epi16(tmp10, tmp13); /* phase 2 */
    tmp3 = _mm_subs_epi16(tmp10, tmp13);
    tmp1 = _mm_adds_epi16(tmp11, tmp12);
    tmp2 = _mm_subs_epi16(tmp11, tmp12);

    /* Odd part */

    z13 = _mm_adds_epi16(row5, row3); /* phase 6 */
    z10 = _mm_subs_epi16(row5, row3);
    z11 = _mm_adds_epi16(row1, row7);
    z12 = _mm_subs_epi16(row1, row7);

    tmp7 = _mm_adds_epi16(z11, z13);  /* phase 5 */

    tmp11 = _mm_mul_epi16_const<1,int(0.414213562*0x8000)>( _mm_subs_epi16(z11,z13) );

    z5 = _mm_mul_epi16_const<1,int(0.847759065*0x8000)>( _mm_adds_epi16(z10,z12) );

    tmp10 = _mm_subs_epi16(z5, _mm_mul_epi16_const<1,int(0.082392200*0x8000)>(z12) );
    tmp12 = _mm_subs_epi16(z5, _mm_mul_epi16_const<2,int(0.613125930*0x8000)>(z10) );

    tmp6 = _mm_subs_epi16(tmp12, tmp7); /* phase 2 */
    tmp5 = _mm_subs_epi16(tmp11, tmp6);
    tmp4 = _mm_subs_epi16(tmp10, tmp5);

    row0 = _mm_adds_epi16(tmp0,tmp7);
    row1 = _mm_adds_epi16(tmp1,tmp6);
    row2 = _mm_adds_epi16(tmp2,tmp5);
    row3 = _mm_adds_epi16(tmp3,tmp4);
    row4 = _mm_subs_epi16(tmp3,tmp4);
    row5 = _mm_subs_epi16(tmp2,tmp5);
    row6 = _mm_subs_epi16(tmp1,tmp6);
    row7 = _mm_subs_epi16(tmp0,tmp7);
    /* transpose here row[0-7] */

    // Step 1: Unpack and interleave adjacent rows
    __m128i t0 = _mm_unpacklo_epi16(row0, row1); // t0: [00, 10, 01, 11]
    __m128i t1 = _mm_unpackhi_epi16(row0, row1); // t1: [02, 12, 03, 13]
    __m128i t2 = _mm_unpacklo_epi16(row2, row3); // t2: [20, 30, 21, 31]
    __m128i t3 = _mm_unpackhi_epi16(row2, row3); // t3: [22, 32, 23, 33]
    __m128i t4 = _mm_unpacklo_epi16(row4, row5); // t4: [40, 50, 41, 51]
    __m128i t5 = _mm_unpackhi_epi16(row4, row5); // t5: [42, 52, 43, 53]
    __m128i t6 = _mm_unpacklo_epi16(row6, row7); // t6: [60, 70, 61, 71]
    __m128i t7 = _mm_unpackhi_epi16(row6, row7); // t7: [62, 72, 63, 73]

    // Step 2: Continue unpacking into 32-bit segments
    __m128i u0 = _mm_unpacklo_epi32(t0, t2); // u0: [00, 10, 20, 30]
    __m128i u1 = _mm_unpackhi_epi32(t0, t2); // u1: [01, 11, 21, 31]
    __m128i u2 = _mm_unpacklo_epi32(t1, t3); // u2: [02, 12, 22, 32]
    __m128i u3 = _mm_unpackhi_epi32(t1, t3); // u3: [03, 13, 23, 33]
    __m128i u4 = _mm_unpacklo_epi32(t4, t6); // u4: [40, 50, 60, 70]
    __m128i u5 = _mm_unpackhi_epi32(t4, t6); // u5: [41, 51, 61, 71]
    __m128i u6 = _mm_unpacklo_epi32(t5, t7); // u6: [42, 52, 62, 72]
    __m128i u7 = _mm_unpackhi_epi32(t5, t7); // u7: [43, 53, 63, 73]

    // Step 3: Unpack into 64-bit segments
    row0 = _mm_unpacklo_epi64(u0, u4); // row[0]: [00, 10, 20, 30, 40, 50, 60, 70]
    row1 = _mm_unpackhi_epi64(u0, u4); // row[1]: [01, 11, 21, 31, 41, 51, 61, 71]
    row2 = _mm_unpacklo_epi64(u1, u5); // row[2]: [02, 12, 22, 32, 42, 52, 62, 72]
    row3 = _mm_unpackhi_epi64(u1, u5); // row[3]: [03, 13, 23, 33, 43, 53, 63, 73]
    row4 = _mm_unpacklo_epi64(u2, u6); // row[4]: [04, 14, 24, 34, 44, 54, 64, 74]
    row5 = _mm_unpackhi_epi64(u2, u6); // row[5]: [05, 15, 25, 35, 45, 55, 65, 75]
    row6 = _mm_unpacklo_epi64(u3, u7); // row[6]: [06, 16, 26, 36, 46, 56, 66, 76]
    row7 = _mm_unpackhi_epi64(u3, u7); // row[7]: [07, 17, 27, 37, 47, 57, 67, 77]

    /* Even part */
    tmp10 = _mm_adds_epi16(row0, row4); /* phase 3 */
    tmp11 = _mm_subs_epi16(row0, row4);

    tmp13 = _mm_adds_epi16(row2, row6); /* phases 5-3 */
    tmp12 = _mm_mul_epi16_const<1,int(0.414213562*0x8000)>(_mm_subs_epi16(row2, row6));
    tmp12 = _mm_subs_epi16(tmp12,tmp13);

    tmp0 = _mm_adds_epi16(tmp10, tmp13); /* phase 2 */
    tmp3 = _mm_subs_epi16(tmp10, tmp13);
    tmp1 = _mm_adds_epi16(tmp11, tmp12);
    tmp2 = _mm_subs_epi16(tmp11, tmp12);

    /* Odd part */

    z13 = _mm_adds_epi16(row5, row3); /* phase 6 */
    z10 = _mm_subs_epi16(row5, row3);
    z11 = _mm_adds_epi16(row1, row7);
    z12 = _mm_subs_epi16(row1, row7);

    tmp7 = _mm_adds_epi16(z11, z13);  /* phase 5 */

    tmp11 = _mm_mul_epi16_const<1,int(0.414213562*0x8000)>( _mm_subs_epi16(z11,z13) );

    z5 = _mm_mul_epi16_const<1,int(0.847759065*0x8000)>( _mm_adds_epi16(z10,z12) );

    tmp10 = _mm_subs_epi16(z5, _mm_mul_epi16_const<1,int(0.082392200*0x8000)>(z12) );
    tmp12 = _mm_subs_epi16(z5, _mm_mul_epi16_const<2,int(0.613125930*0x8000)>(z10) );

    tmp6 = _mm_subs_epi16(tmp12, tmp7); /* phase 2 */
    tmp5 = _mm_subs_epi16(tmp11, tmp6);
    tmp4 = _mm_subs_epi16(tmp10, tmp5);

    row0 = _mm_adds_epi16(tmp0,tmp7);
    row1 = _mm_adds_epi16(tmp1,tmp6);
    row2 = _mm_adds_epi16(tmp2,tmp5);
    row3 = _mm_adds_epi16(tmp3,tmp4);
    row4 = _mm_subs_epi16(tmp3,tmp4);
    row5 = _mm_subs_epi16(tmp2,tmp5);
    row6 = _mm_subs_epi16(tmp1,tmp6);
    row7 = _mm_subs_epi16(tmp0,tmp7);

    row0 = _mm_srai_epi16( row0, fix_pass);
    row1 = _mm_srai_epi16( row1, fix_pass);
    row2 = _mm_srai_epi16( row2, fix_pass);
    row3 = _mm_srai_epi16( row3, fix_pass);
    row4 = _mm_srai_epi16( row4, fix_pass);
    row5 = _mm_srai_epi16( row5, fix_pass);
    row6 = _mm_srai_epi16( row6, fix_pass);
    row7 = _mm_srai_epi16( row7, fix_pass);
    row0 = _mm_packs_epi16( row0, row1 );
    row1 = _mm_set1_epi8(-128);
    row2 = _mm_packs_epi16( row2, row3 );
    row4 = _mm_packs_epi16( row4, row5 );
    row6 = _mm_packs_epi16( row6, row7 );
    row0 = _mm_xor_si128( row0, row1);
    row2 = _mm_xor_si128( row2, row1);
    row4 = _mm_xor_si128( row4, row1);
    row6 = _mm_xor_si128( row6, row1);

    alignas(16) uint8_t m16[8*2];

    _mm_stream_si128( (__m128i*) m16, row0);
    std::copy(m16,m16+8, out); out+=stride; std::copy(m16+8,m16+16, out); out+=stride;
    _mm_stream_si128( (__m128i*) m16, row2);
    std::copy(m16,m16+8, out); out+=stride; std::copy(m16+8,m16+16, out); out+=stride;
    _mm_stream_si128( (__m128i*) m16, row4);
    std::copy(m16,m16+8, out); out+=stride; std::copy(m16+8,m16+16, out); out+=stride;
    _mm_stream_si128( (__m128i*) m16, row6);
    std::copy(m16,m16+8, out); out+=stride; std::copy(m16+8,m16+16, out); out+=stride;

}

alignas(32) float block_ps[64]{1024.0f};
alignas(16) int16_t block_i[64]{1024};
uint8_t out[320*200];

void idct_0(benchmark::State& state) {
    for (auto _ : state) {
        // Приостановить замер времени для инициализации
        state.PauseTiming();
        static float value{1.0f};

        // Инициализация данных перед каждой итерацией
        std::fill(std::begin(block_ps), std::end(block_ps), value);
        value += 0.001;

        state.ResumeTiming();  // Возобновить замер времени

        // Замер основного кода
        idct8x8(block_ps, out, 320);

        benchmark::DoNotOptimize(out);
        benchmark::ClobberMemory();
    }
}

void idct_1(benchmark::State& state) {
    for (auto _ : state) {
        state.PauseTiming();
        static int16_t value{1};

        // Инициализация данных перед каждой итерацией
        std::fill(std::begin(block_i), std::end(block_i), value);
        value += 1;

        state.ResumeTiming();

        // Замер основного кода
        idct8x8_epi16(block_i, out, 320);

        benchmark::DoNotOptimize(out);
        benchmark::ClobberMemory();
    }
}


BENCHMARK(idct_0);
BENCHMARK(idct_1);
BENCHMARK_MAIN();

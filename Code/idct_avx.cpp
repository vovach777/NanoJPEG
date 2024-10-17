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

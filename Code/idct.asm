%ifdef WIN64
    %define ARG1 rcx
    %define ARG2 rdx
    %define ARG3 r8
    %define ARG4 r9
%else
    %define ARG1 rdi
    %define ARG2 rsi
    %define ARG3 rdx
    %define ARG4 rcx
    %define ARG5 r8
    %define ARG6 r9
%endif
default rel
section .text
bits 64
global idct8x8_epi16
;SSSE3 + SSE4.1
;idct8x8_epi16(short const* block_aligned16_transposed_shifted_left_by_4_AAN, unsigned char* out, int stride):
idct8x8_epi16:
        movdqa  xmm0, [ARG1 + 64]
        movdqa  xmm1, [ARG1 + 96]
        movdqa  xmm5, [ARG1]
        movdqa  xmm3, [ARG1 + 16]
        movdqa  xmm2, [ARG1 + 32]
        movdqa  xmm9, [ARG1 + 48]
        movdqa  xmm4, xmm5
        paddsw  xmm4, xmm0
        psubsw  xmm5, xmm0
        movdqa  xmm6, xmm2
        paddsw  xmm6, xmm1
        psubsw  xmm2, xmm1
        movdqa  xmm0, [LCPI0_0]
        movdqa  xmm1, xmm2
        pmulhrsw xmm1, xmm0
        paddsw  xmm1, xmm2
        psubsw  xmm1, xmm6
        movdqa  xmm8, xmm4
        paddsw  xmm8, xmm6
        psubsw  xmm4, xmm6
        movdqa  xmm6, xmm5
        paddsw  xmm6, xmm1
        psubsw  xmm5, xmm1
        movdqa  xmm7, [ARG1 + 80]
        movdqa  xmm1, [ARG1 + 112]
        movdqa  xmm2, xmm7
        paddsw  xmm2, xmm9
        psubsw  xmm7, xmm9
        movdqa  xmm11, xmm3
        paddsw  xmm11, xmm1
        psubsw  xmm3, xmm1
        movdqa  xmm9, xmm11
        paddsw  xmm9, xmm2
        psubsw  xmm11, xmm2
        movdqa  xmm10, xmm11
        pmulhrsw xmm10, xmm0
        paddsw  xmm10, xmm11
        movdqa  xmm2, xmm7
        paddsw  xmm2, xmm3
        movdqa  xmm1, [LCPI0_1]
        movdqa  xmm11, xmm2
        pmulhrsw xmm11, xmm1
        paddsw  xmm11, xmm2
        movdqa  xmm2, [LCPI0_2]
        movdqa  xmm12, xmm3
        pmulhrsw xmm12, xmm2
        paddsw  xmm12, xmm3
        movdqa  xmm13, xmm11
        psubsw  xmm13, xmm12
        movdqa  xmm3, [LCPI0_3]
        movdqa  xmm12, xmm7
        pmulhrsw xmm12, xmm3
        paddsw  xmm7, xmm7
        paddsw  xmm7, xmm12
        psubsw  xmm11, xmm7
        psubsw  xmm11, xmm9
        psubsw  xmm10, xmm11
        psubsw  xmm13, xmm10
        movdqa  xmm7, xmm8
        paddsw  xmm7, xmm9
        movdqa  xmm14, xmm6
        paddsw  xmm14, xmm11
        movdqa  xmm12, xmm5
        paddsw  xmm12, xmm10
        movdqa  xmm15, xmm4
        paddsw  xmm15, xmm13
        psubsw  xmm4, xmm13
        psubsw  xmm5, xmm10
        psubsw  xmm6, xmm11
        psubsw  xmm8, xmm9
        movdqa  xmm11, xmm7
        punpcklwd   xmm11, xmm14
        punpckhwd   xmm7, xmm14
        movdqa      xmm9, xmm12
        punpcklwd   xmm9, xmm15
        punpckhwd   xmm12, xmm15
        movdqa      xmm13, xmm4
        punpcklwd   xmm13, xmm5
        punpckhwd   xmm4, xmm5
        movdqa      xmm5, xmm6
        punpcklwd   xmm5, xmm8
        punpckhwd   xmm6, xmm8
        movdqa      xmm10, xmm11
        punpckldq   xmm10, xmm9
        punpckhdq   xmm11, xmm9
        movdqa      xmm9, xmm7
        punpckldq   xmm9, xmm12
        punpckhdq   xmm7, xmm12
        movdqa      xmm8, xmm13
        punpckldq   xmm8, xmm5
        punpckhdq   xmm13, xmm5
        movdqa      xmm5, xmm4
        punpckldq   xmm5, xmm6
        punpckhdq   xmm4, xmm6
        movdqa      xmm6, xmm10
        punpcklqdq  xmm6, xmm8
        punpckhqdq  xmm10, xmm8
        movdqa      xmm8, xmm11
        punpcklqdq  xmm8, xmm13
        punpckhqdq  xmm11, xmm13
        movdqa      xmm12, xmm9
        punpcklqdq  xmm12, xmm5
        punpckhqdq  xmm9, xmm5
        movdqa      xmm5, xmm7
        punpcklqdq  xmm5, xmm4
        punpckhqdq  xmm7, xmm4
        movdqa  xmm4, xmm6
        paddsw  xmm4, xmm12
        psubsw  xmm6, xmm12
        movdqa  xmm12, xmm8
        paddsw  xmm12, xmm5
        psubsw  xmm8, xmm5
        movdqa  xmm13, xmm8
        pmulhrsw        xmm13, xmm0
        paddsw  xmm13, xmm8
        psubsw  xmm13, xmm12
        movdqa  xmm8, xmm4
        paddsw  xmm8, xmm12
        psubsw  xmm4, xmm12
        movdqa  xmm5, xmm6
        paddsw  xmm5, xmm13
        psubsw  xmm6, xmm13
        movdqa  xmm12, xmm9
        paddsw  xmm12, xmm11
        psubsw  xmm9, xmm11
        movdqa  xmm11, xmm10
        paddsw  xmm11, xmm7
        psubsw  xmm10, xmm7
        movdqa  xmm7, xmm11
        paddsw  xmm7, xmm12
        psubsw  xmm11, xmm12
        pmulhrsw        xmm0, xmm11
        paddsw  xmm0, xmm11
        movdqa  xmm11, xmm9
        paddsw  xmm11, xmm10
        pmulhrsw        xmm1, xmm11
        paddsw  xmm1, xmm11
        pmulhrsw        xmm2, xmm10
        paddsw  xmm2, xmm10
        movdqa  xmm10, xmm1
        psubsw  xmm10, xmm2
        pmulhrsw        xmm3, xmm9
        paddsw  xmm9, xmm9
        paddsw  xmm9, xmm3
        psubsw  xmm1, xmm9
        psubsw  xmm1, xmm7
        psubsw  xmm0, xmm1
        psubsw  xmm10, xmm0
        movdqa  xmm3, xmm8
        paddsw  xmm3, xmm7
        movdqa  xmm9, xmm5
        paddsw  xmm9, xmm1
        movdqa  xmm2, xmm6
        paddsw  xmm2, xmm0
        movdqa  xmm11, xmm4
        paddsw  xmm11, xmm10
        psubsw  xmm4, xmm10
        psubsw  xmm6, xmm0
        psubsw  xmm5, xmm1
        psubsw  xmm8, xmm7
        psraw   xmm3, 4
        psraw   xmm9, 4
        packsswb        xmm3, xmm9
        psraw   xmm2, 4
        psraw   xmm11, 4
        packsswb        xmm2, xmm11
        psraw   xmm4, 4
        psraw   xmm6, 4
        packsswb        xmm4, xmm6
        psraw   xmm5, 4
        psraw   xmm8, 4
        packsswb        xmm5, xmm8
        movdqa  xmm0, [LCPI0_4]
        pxor    xmm3, xmm0
        movq    [ARG2], xmm3
        add     ARG2,ARG3
        pextrq  [ARG2], xmm3, 1
        add     ARG2,ARG3
        pxor    xmm2, xmm0
        movq    [ARG2], xmm2
        add     ARG2,ARG3
        pextrq  [ARG2], xmm2, 1
        add     ARG2,ARG3
        pxor    xmm4, xmm0
        movq    [ARG2], xmm4
        add     ARG2, ARG3
        pextrq  [ARG2], xmm4, 1
        add     ARG2, ARG3
        pxor    xmm5, xmm0
        movq    [ARG2], xmm5
        add     ARG2, ARG3
        pextrq  [ARG2], xmm5, 1
        ret

section .data
    LCPI0_0:
        dw  13572, 13572, 13572, 13572, 13572, 13572, 13572, 13572
    LCPI0_1:
        dw  27779, 27779, 27779, 27779, 27779, 27779, 27779, 27779
    LCPI0_2:
        dw  2699, 2699, 2699, 2699, 2699, 2699, 2699, 2699
    LCPI0_3:
        dw  20090, 20090, 20090, 20090, 20090, 20090, 20090, 20090
    LCPI0_4:
        times 16 db 0x80

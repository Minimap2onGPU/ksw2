#include <stdio.h>
#include <stdlib.h>
#include <smmintrin.h>
#include <math.h>
#include <assert.h>

#ifdef __SSE2__
#ifdef USE_SIMDE
#include <simde/x86/sse2.h>
#else
#include <emmintrin.h>
#endif

#ifdef KSW_SSE2_ONLY
#undef __SSE4_1__
#endif

#ifdef __SSE4_1__
#ifdef USE_SIMDE
#include <simde/x86/sse4.1.h>
#else
#include <smmintrin.h>
#endif
#endif

__inline__ __m128i copy(const vec in) {
    __m128i out = _mm_load_si128(&in);
    return out;
}

__inline__ void printm128(const __m128i a) {
    u_int8_t *x = (u_int8_t*)&a;
    for (int i = 0; i < size; ++i) {
        printf("%u ", x[i]);
    }
    printf("\n");

}

__inline__ void print(const vec *a) {
    for (int i = 0; i < size; ++i) {
        printf("%u ", a->x[i]);
    }
    printf("\n");

}

#include "vector.h"

#define TEST() printf("%32s ", __func__);

__inline__ void compare(const __m128i *a, const vec *b) {
    u_int8_t *x = (u_int8_t*)a;
    for (int i = 0; i < size; ++ i) {
        if (x[i] != b->x[i]) {
            printf("a[%d] = %d b[%d] = %d \n", i, x[i], i, b->x[i]);
        }
        assert(x[i] == b->x[i]);
    }
    printf("OK\n");

}

void init_random(int *a, int n) {
    for (int i = 0; i < n; ++i) {
        a[i] = rand();
    }
}

int test_load(const int *a) {
    TEST()
    __m128i x = _mm_load_si128(a);
    vec y = load(a);
    compare(&x, &y);

    return 0;
}

int test_loadu(const int *a) {
    TEST()
    __m128i x = _mm_loadu_si128(a);
    vec y = load(a);
    compare(&x, &y);

    return 0;
}

// 1100111
// 1000101

int test_shift_right(const int *a) {
    TEST()
    __m128i x = _mm_load_si128(a);
    x = _mm_srli_si128(x, 15);

    vec y = load(a);
    y = shift_right(*(vec*)a);
    compare(&x, &y);

    return 0;
}

int test_shift_left(const int *a) {
    TEST()
    __m128i x = _mm_load_si128(a);
    x = _mm_slli_si128(x, 1);

    vec y = load(a);
    y = shift_left(*(vec*)a);
    compare(&x, &y);

    return 0;
}

int test_bit_or(const int *a, const int *b) {
    TEST()
    __m128i x1 = _mm_load_si128(a);
    __m128i x2 = _mm_load_si128(b);
    __m128i x = _mm_or_si128(x1, x2);

    vec y = bit_or(*(vec*)a, *(vec*)b);
    compare(&x, &y);

    return 0;
}

int test_bit_and(const int *a, const int *b) {
    TEST()
    __m128i x1 = _mm_load_si128(a);
    __m128i x2 = _mm_load_si128(b);
    __m128i x = _mm_and_si128(x1, x2);

    vec y = bit_and(*(vec*)a, *(vec*)b);
    compare(&x, &y);

    return 0;
}

int test_add(const int *a, const int *b) {
    TEST()
    __m128i x1 = _mm_load_si128(a);
    __m128i x2 = _mm_load_si128(b);
    __m128i x = _mm_add_epi8(x1, x2);

    vec y = add(*(vec*)a, *(vec*)b);
    compare(&x, &y);

    return 0;
}

int test_add32(const int *a, const int *b) {
    TEST()
    __m128i x1 = _mm_load_si128(a);
    __m128i x2 = _mm_load_si128(b);
    __m128i x = _mm_add_epi32(x1, x2);

    vec y = add32(*(vec*)a, *(vec*)b);
    compare(&x, &y);

    return 0;
}

int test_sub(const int *a, const int *b) {
    TEST()
    __m128i x1 = _mm_load_si128(a);
    __m128i x2 = _mm_load_si128(b);
    __m128i x = _mm_sub_epi8(x1, x2);

    vec y = sub(*(vec*)a, *(vec*)b);
    compare(&x, &y);

    return 0;
}

int test_sub32(const int *a, const int *b) {
    TEST()
    __m128i x1 = _mm_load_si128(a);
    __m128i x2 = _mm_load_si128(b);
    __m128i x = _mm_sub_epi32(x1, x2);

    vec y = sub32(*(vec*)a, *(vec*)b);
    compare(&x, &y);

    return 0;
}

int test_andnot(const int *a, const int *b) {
    TEST()
    __m128i x1 = _mm_load_si128(a);
    __m128i x2 = _mm_load_si128(b);
    __m128i x = _mm_andnot_si128(x1, x2);

    vec y = andnot(*(vec*)a, *(vec*)b);
    compare(&x, &y);

    return 0;
}

int test_cmpeq(const int *a, const int *b) {
    TEST()
    __m128i x1 = _mm_load_si128(a);
    __m128i x2 = _mm_load_si128(b);
    __m128i x = _mm_cmpeq_epi8(x1, x2);

    vec y = cmpeq(*(vec*)a, *(vec*)b);
    compare(&x, &y);

    return 0;
}

int test_cmpgt(const int *a, const int *b) {
    TEST()
    __m128i x1 = _mm_load_si128(a);
    __m128i x2 = _mm_load_si128(b);
    __m128i x = _mm_cmpgt_epi8(x1, x2);

    vec y = cmpgt(*(vec*)a, *(vec*)b);
    compare(&x, &y);

    return 0;
}

int test_cmpgt32(const int *a, const int *b) {
    TEST()
    __m128i x1 = _mm_load_si128(a);
    __m128i x2 = _mm_load_si128(b);
    __m128i x = _mm_cmpgt_epi32(x1, x2);

    vec y = cmpgt32(*(vec*)a, *(vec*)b);
    compare(&x, &y);

    return 0;
}

int test_cmplt(const int *a, const int *b) {
    TEST()
    __m128i x1 = _mm_load_si128(a);
    __m128i x2 = _mm_load_si128(b);
    __m128i x = _mm_cmplt_epi8(x1, x2);

    vec y = cmplt(*(vec*)a, *(vec*)b);
    compare(&x, &y);

    return 0;
}

int test_cvtsi32(const int *a) {
    TEST()
    __m128i x = _mm_cvtsi32_si128(a[0]);
    vec y = cvtsi32(a[0]);
    compare(&x, &y);

    return 0;
}

 int test_set1(const int *a) {
    TEST()
    u_int8_t val = (u_int8_t)a[0];
    __m128i x = _mm_set1_epi8(val);
    vec y = set1(val);
    compare(&x, &y);

    return 0;
}

int test_set32(const int *a) {
    TEST()
    u_int8_t val = (u_int8_t)a[0];
    __m128i x = _mm_set1_epi8(val);
    vec y = set1(val);
    compare(&x, &y);

    return 0;
}

int test_setr32(const int *a) {
    TEST()
    __m128i x = _mm_setr_epi32(a[0], a[1], a[2], a[3]);
    vec y = setr32(a[0], a[1], a[2], a[3]);
    compare(&x, &y);

    return 0;
}

int test_store(const int *a) {
    TEST()
    __m128i x;
    _mm_store_si128(&x, ((__m128i*)a)[0]);
    vec y;
    store(&y, ((vec*)a)[0]);
    compare(&x, &y);

    return 0;
}

int test_storeu(const int *a) {
    TEST()
    __m128i x;
    _mm_storeu_si128(&x, ((__m128i*)a)[0]);
    vec y;
    store(&y, ((vec*)a)[0]);
    compare(&x, &y);

    return 0;
}

int test_copy(const int *a) {
    TEST()
    vec x = load(a);
    __m128i y = copy(x);

    compare(&x, &y);

    return 0;
}




int main() {

    int a[size];
    int b[size];
    for (int i = 0; i < 100; ++i) {
        printf("Test case: %d \n", i);
        init_random(a, size);
        init_random(b, size);

        test_load(a);
        test_shift_right(a);
        test_shift_left(a);
        test_bit_or(a, b);
        test_bit_and(a, b);
        test_add(a, b);
        test_add32(a, b);
        test_sub(a, b);
        test_sub32(a, b);
        test_andnot(a, b);
        test_cmpeq(a, b);
        test_cmpgt(a, b);
        test_cmpgt32(a, b);
        test_cmplt(a, b);
        test_cvtsi32(a);
        test_set1(a);
        test_set32(a);
        test_loadu(a);
        test_setr32(a);
        test_store(a);
        test_storeu(a);
        test_copy(a);
    }

}

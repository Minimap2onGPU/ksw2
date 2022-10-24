#pragma once
#include <limits.h>

#define size 16
#define size32 4
typedef union {
    u_int8_t x[size];
    u_int32_t y[size32];
    int8_t z[size];
    int32_t w[size32];
} vec;


// _mm_load_si128
__inline__ vec load(const vec *v) {
    vec out;
    for (int i = 0; i < size; ++i)
        out.x[i] = v->x[i];
    return out;
}

// _mm_srli_si128 with imm8 = size - 1
__inline__ vec shift_right(const vec v) {
    vec out = { .x = {0}};
    out.x[0] = v.x[size - 1];
    return out;
}


// _mm_slli_si128 with imm8 = 1
__inline__ vec shift_left(const vec v) {
    vec out = { .x = {0}};
    out.x[0] = 0;
    for (int i = 1; i < size; ++i)
        out.x[i] = v.x[i-1];
    return out;
}


//  _mm_or_si128
__inline__ vec bit_or(const vec a, const vec b) {
    vec out;
    for (int i = 0; i < size; ++ i) {
        out.x[i] = a.x[i] | b.x[i];
    }
    return out;
}

// _mm_add_epi8
__inline__ vec add(const vec a, const vec b) {
    vec out;
    for (int i = 0; i < size; ++ i) {
        out.x[i] = a.x[i] + b.x[i];
    }
    return out;
}

// _mm_add_epi32
__inline__ vec add32(const vec a, const vec b) {
    vec out;
    for (int i = 0; i < size32; ++ i) {
        out.y[i] = a.y[i] + b.y[i];
    }
    return out;
}
 
// _mm_sub_epi8
__inline__ vec sub(const vec a, const vec b) {
    vec out;
    for (int i = 0; i < size; ++ i) {
        out.x[i] = a.x[i] - b.x[i];
    }
    return out;
}

// _mm_sub_epi32
__inline__ vec sub32(const vec a, const vec b) {
    vec out;
    for (int i = 0; i < size32; ++ i) {
        out.y[i] = a.y[i] - b.y[i];
    }
    return out;
}

// _mm_andnot_si128
__inline__ vec andnot(const vec a, const vec b) {
    vec out;
    for (int i = 0; i < size; ++ i) {
        out.x[i] = ~a.x[i] & b.x[i];
    }
    return out;
}

// _mm_and_si128
__inline__ vec bit_and(const vec a, const vec b) {
    vec out;
    for (int i = 0; i < size; ++ i) {
        out.x[i] = a.x[i] & b.x[i];
    }
    return out;
}

// _mm_cmpeq_epi8
__inline__ vec cmpeq(const vec a, const vec b) {
    vec out;
    for (int i = 0; i < size; ++ i) {
        out.x[i] = (a.x[i] == b.x[i]) ? 0xFF : 0;
    }
    return out;
}

// _mm_cmpgt_epi8
__inline__ vec cmpgt(const vec a, const vec b) {
    vec out;
    for (int i = 0; i < size; ++ i) {
        out.z[i] = (a.z[i] > b.z[i]) ? 0xFF : 0;
    }
    return out;
}

// _mm_cmpgt_epi32
__inline__ vec cmpgt32(const vec a, const vec b) {
    vec out;
    for (int i = 0; i < size32; ++ i) {
        out.w[i] = (a.w[i] > b.w[i]) ? 0xFFFFFFFF : 0;
    }
    return out;
}


// _mm_cmplt_epi8
__inline__ vec cmplt(const vec a, const vec b) {
    vec out;
    for (int i = 0; i < size; ++ i) {
        out.z[i] = (a.z[i] < b.z[i]) ? 0xFF : 0;
    }
    return out;
}

// _mm_cvtsi32_si128
__inline__ vec cvtsi32(const int a) {
    vec out = {.y = {0}};
    out.y[0] = a;
    return out;
}


// _mm_set1_epi8
__inline__ vec set1(const u_int8_t a) {
    vec out;
    for (int i = 0; i < size; ++i)
        out.x[i] = a;
    return out;
}

// _mm_set1_epi32
__inline__ vec set32(const int a) {
    vec out;
    for (int i = 0; i < size32; ++i)
        out.y[i] = a;
    return out;
}

// _mm_setr_epi32 (int e3, int e2, int e1, int e0)
__inline__ vec setr32(int e3, int e2, int e1, int e0) {
    vec out;
    out.y[0] = e3;
    out.y[1] = e2;
    out.y[2] = e1;
    out.y[3] = e0;
    return out;
}

// _mm_store_si128
__inline__ void store(vec *out, const vec in) {
    for (int i = 0; i < size32; ++i)
        out->y[i] = in.y[i];
}

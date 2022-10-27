#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "ksw2.h"

#include <vector.h>

void ksw_extd2_no_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
                   int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
{
#define __dp_code_block1 \
    vz = load(&s[t]); \
    vxt1 = load(&x[t]);                     /* xt1 <- x[r-1][t..t+15] */ \
    vvt1 = load(&v[t]); \
    vut = load(&u[t]); \
    vx2t1 = load(&x2[t]); \
    vtmp = shift_right(vxt1);                   /* tmp <- x[r-1][t+15] */ \
    vxt1 = bit_or(shift_left(vxt1), vx1_); \
    vx1_ = vtmp; \
    vtmp = shift_right(vvt1); \
    vvt1 = bit_or(shift_left(vvt1), vv1_); \
    vv1_ = vtmp; \
    va = add(vxt1, vvt1); \
    vb = add(load(&y[t]), vut); \
    vtmp = shift_right(vx2t1); \
    vx2t1 = bit_or(shift_left(vx2t1), vx21_); \
    vx21_= vtmp; \
    va2 = add(vx2t1, vvt1); \
    vb2 = add(load(&y2[t]), vut); 

#define __dp_code_block2 \
    store(&u[t], sub(vz, vvt1));    /* u[r][t..t+15] <- z - v[r-1][t-1..t+14] */ \
    store(&v[t], sub(vz, vut));     /* v[r][t..t+15] <- z - u[r-1][t..t+15] */ \
    vtmp = sub(vz, vq_); \
    va = sub(va, vtmp); \
    vb = sub(vb, vtmp); \
    vtmp = sub(vz, vq2_); \
    va2= sub(va2, vtmp); \
    vb2= sub(vb2, vtmp); \

    int r, t, qe = q + e, n_col_, *off = 0, *off_end = 0, tlen_, qlen_, last_st, last_en, wl, wr, max_sc, min_sc, long_thres, long_diff;
    int with_cigar = !(flag&KSW_EZ_SCORE_ONLY), approx_max = !!(flag&KSW_EZ_APPROX_MAX);
    int32_t *H = 0, H0 = 0, last_H0_t = 0;
    uint8_t *qr, *sf, *mem, *mem2 = 0;
    vec vq_, vq2_, vqe_, vqe2_, vzero_, vsc_mch_, vsc_mis_, vm1_, vsc_N_;
    vec *u, *v, *x, *y, *x2, *y2, *s, *p = 0;

    ksw_reset_extz(ez);
    if (m <= 1 || qlen <= 0 || tlen <= 0) return;

    if (q2 + e2 < q + e) t = q, q = q2, q2 = t, t = e, e = e2, e2 = t; // make sure q+e no larger than q2+e2

    vzero_   = set1(0);
    vq_      = set1(q);
    vq2_     = set1(q2);
    vqe_     = set1(q + e);
    vqe2_    = set1(q2 + e2);
    vsc_mch_ = set1(mat[0]);
    vsc_mis_ = set1(mat[1]);
    vsc_N_   = mat[m*m-1] == 0? set1(-e2) : set1(mat[m*m-1]);
    vm1_     = set1(m - 1); // wildcard


    if (w < 0) w = tlen > qlen? tlen : qlen;
    wl = wr = w;
    tlen_ = (tlen + 15) / 16;
    n_col_ = qlen < tlen? qlen : tlen;
    n_col_ = ((n_col_ < w + 1? n_col_ : w + 1) + 15) / 16 + 1;
    qlen_ = (qlen + 15) / 16;
    for (t = 1, max_sc = mat[0], min_sc = mat[1]; t < m * m; ++t) {
        max_sc = max_sc > mat[t]? max_sc : mat[t];
        min_sc = min_sc < mat[t]? min_sc : mat[t];
    }
    if (-min_sc > 2 * (q + e)) return; // otherwise, we won't see any mismatches

    long_thres = e != e2? (q2 - q) / (e - e2) - 1 : 0;
    if (q2 + e2 + long_thres * e2 > q + e + long_thres * e)
        ++long_thres;
    long_diff = long_thres * (e - e2) - (q2 - q) - e2;

    mem = (uint8_t*)kcalloc(km, tlen_ * 8 + qlen_ + 1, 16);
    u = (vec*)(((size_t)mem + 15) >> 4 << 4); // 16-byte aligned
    v = u + tlen_, x = v + tlen_, y = x + tlen_, x2 = y + tlen_, y2 = x2 + tlen_;
    s = y2 + tlen_, sf = (uint8_t*)(s + tlen_), qr = sf + tlen_ * 16;
    memset(u,  -q  - e,  tlen_ * 16);
    memset(v,  -q  - e,  tlen_ * 16);
    memset(x,  -q  - e,  tlen_ * 16);
    memset(y,  -q  - e,  tlen_ * 16);
    memset(x2, -q2 - e2, tlen_ * 16);
    memset(y2, -q2 - e2, tlen_ * 16);

    if (!approx_max) {
        H = (int32_t*)kmalloc(km, tlen_ * 16 * 4);
        for (t = 0; t < tlen_ * 16; ++t) H[t] = KSW_NEG_INF;
    }
    if (with_cigar) {
        mem2 = (uint8_t*)kmalloc(km, ((size_t)(qlen + tlen - 1) * n_col_ + 1) * 16);
        p = (vec*)(((size_t)mem2 + 15) >> 4 << 4);
        off = (int*)kmalloc(km, (qlen + tlen - 1) * sizeof(int) * 2);
        off_end = off + qlen + tlen - 1;
    }

    for (t = 0; t < qlen; ++t) qr[t] = query[qlen - 1 - t];
    memcpy(sf, target, tlen);

    for (r = 0, last_st = last_en = -1; r < qlen + tlen - 1; ++r) {
        int st = 0, en = tlen - 1, st0, en0, st_, en_;
        int8_t x1, x21, v1;
        uint8_t *qrr = qr + (qlen - 1 - r);
        int8_t *u8 = (int8_t*)u, *v8 = (int8_t*)v, *x8 = (int8_t*)x, *x28 = (int8_t*)x2;
        vec vx1_, vx21_, vv1_;
        // find the boundaries
        if (st < r - qlen + 1) st = r - qlen + 1;
        if (en > r) en = r;
        if (st < (r-wr+1)>>1) st = (r-wr+1)>>1; // take the ceil
        if (en > (r+wl)>>1) en = (r+wl)>>1; // take the floor
        if (st > en) {
            ez->zdropped = 1;
            break;
        }
        st0 = st, en0 = en;
        st = st / 16 * 16, en = (en + 16) / 16 * 16 - 1;
        // set boundary conditions
        if (st > 0) {
            if (st - 1 >= last_st && st - 1 <= last_en) {
                x1 = x8[st - 1], x21 = x28[st - 1], v1 = v8[st - 1]; // (r-1,s-1) calculated in the last round
            } else {
                x1 = -q - e, x21 = -q2 - e2;
                v1 = -q - e;
            }
        } else {
            x1 = -q - e, x21 = -q2 - e2;
            v1 = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : -e2;
        }
        if (en >= r) {
            ((int8_t*)y)[r] = -q - e, ((int8_t*)y2)[r] = -q2 - e2;
            u8[r] = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : -e2;
        }
        // loop fission: set scores first
        if (!(flag & KSW_EZ_GENERIC_SC)) {
            for (t = st0; t <= en0; t += 16) {
                vec vsq, vst, vtmp, vmask;
                vsq = load(&sf[t]);
                vst = load(&qrr[t]);
                vmask = bit_or(cmpeq(vsq, vm1_), cmpeq(vst, vm1_));
                vtmp = cmpeq(vsq, vst);
                vtmp = bit_or(andnot(vtmp,  vsc_mis_), bit_and(vtmp,  vsc_mch_));
                vtmp = bit_or(andnot(vmask, vtmp),     bit_and(vmask, vsc_N_));
                store((vec*)((int8_t*)s + t), vtmp);
            }
        } else {
            for (t = st0; t <= en0; ++t)
                ((uint8_t*)s)[t] = mat[sf[t] * m + qrr[t]];
        }
        // core loop
        vx1_  = cvtsi32((uint8_t)x1);
        vx21_ = cvtsi32((uint8_t)x21);
        vv1_  = cvtsi32((uint8_t)v1);

        st_ = st / 16, en_ = en / 16;
        assert(en_ - st_ + 1 <= n_col_);
        if (!with_cigar) { // score only
            for (t = st_; t <= en_; ++t) {
                vec vz, va, vb, va2, vb2, vxt1, vx2t1, vvt1, vut, vtmp;
                __dp_code_block1;

                vtmp = cmpgt(va, vz); 
                vz = bit_or(andnot(vtmp, vz), bit_and(vtmp, va));

                vtmp = cmpgt(vb,  vz);
                vz = bit_or(andnot(vtmp, vz), bit_and(vtmp, vb));

                vtmp = cmpgt(va2,  vz);
                vz = bit_or(andnot(vtmp, vz), bit_and(vtmp, va2));

                vtmp = cmpgt(vb2, vz);
                vz = bit_or(andnot(vtmp, vz), bit_and(vtmp, vb2));

                vtmp = cmplt(vsc_mch_, vz);
                vz = bit_or(bit_and(vtmp, vsc_mch_), andnot(vtmp, vz));

                __dp_code_block2;

                vtmp = cmpgt(va, vzero_);
                store(&x[t], sub(bit_and(vtmp, va),  vqe_));

                vtmp = cmpgt(vb, vzero_);
                store(&y[t],  sub(bit_and(vtmp, vb),  vqe_));

                vtmp = cmpgt(va2, vzero_);
                store(&x2[t], sub(bit_and(vtmp, va2), vqe2_));

                vtmp = cmpgt(vb2, vzero_);
                store(&y2[t], sub(bit_and(vtmp, vb2), vqe2_));

            }
        } else if (!(flag&KSW_EZ_RIGHT)) { // gap left-alignment
            vec *pr = p + (size_t)r * n_col_ - st_;
            off[r] = st, off_end[r] = en;
            for (t = st_; t <= en_; ++t) {
                vec vd, vz, va, vb, va2, vb2, vxt1, vx2t1, vvt1, vut, vtmp;
                __dp_code_block1;

                vtmp = cmpgt(va,  vz);

                vd = bit_and(vtmp, set1(1));
                vz = bit_or(andnot(vtmp, vz), bit_and(vtmp, va));
                vtmp = cmpgt(vb,  vz);
                
                vd = bit_or(andnot(vtmp, vd), bit_and(vtmp, set1(2)));
                vz = bit_or(andnot(vtmp, vz), bit_and(vtmp, vb));
                vtmp = cmpgt(va2, vz);

                vd = bit_or(andnot(vtmp, vd), bit_and(vtmp, set1(3)));
                vz = bit_or(andnot(vtmp, vz), bit_and(vtmp, va2));
                vtmp = cmpgt(vb2, vz);
                
                vd = bit_or(andnot(vtmp, vd), bit_and(vtmp, set1(4)));
                vz = bit_or(andnot(vtmp, vz), bit_and(vtmp, vb2));
                vtmp = cmplt(vsc_mch_, vz);
                
                vz = bit_or(bit_and(vtmp, vsc_mch_), andnot(vtmp, vz));

                __dp_code_block2;
                
                vtmp = cmpgt(va, vzero_);
                store(&x[t],  sub(bit_and(vtmp, va),  vqe_));
                
                vd = bit_or(vd, bit_and(vtmp, set1(0x08))); // d = a > 0? 1<<3 : 0
                vtmp = cmpgt(vb, vzero_);
                store(&y[t],  sub(bit_and(vtmp, vb),  vqe_));
                
                vd = bit_or(vd, bit_and(vtmp, set1(0x10))); // d = b > 0? 1<<4 : 0
                vtmp = cmpgt(va2, vzero_);
                store(&x2[t], sub(bit_and(vtmp, va2), vqe2_));
                
                vd = bit_or(vd, bit_and(vtmp, set1(0x20))); // d = a > 0? 1<<5 : 0
                vtmp = cmpgt(vb2, vzero_);
                store(&y2[t], sub(bit_and(vtmp, vb2), vqe2_));

                vd = bit_or(vd, bit_and(vtmp, set1(0x40))); // d = b > 0? 1<<6 : 0
                store(&pr[t], vd);
            }
        } else { // gap right-alignment
            vec *pr = p + (size_t)r * n_col_ - st_;
            off[r] = st, off_end[r] = en;
            for (t = st_; t <= en_; ++t) {
                vec vd, vz, va, vb, va2, vb2, vxt1, vx2t1, vvt1, vut, vtmp;
                __dp_code_block1;

                vtmp = cmpgt(vz, va);
                vd = andnot(vtmp, set1(1));
                vz = bit_or(bit_and(vtmp, vz), andnot(vtmp, va));
                vtmp = cmpgt(vz, vb);
                
                vd = bit_or(bit_and(vtmp, vd), andnot(vtmp, set1(2)));
                vz = bit_or(bit_and(vtmp, vz), andnot(vtmp, vb));
                vtmp = cmpgt(vz, va2);
                
                vd = bit_or(bit_and(vtmp, vd), andnot(vtmp, set1(3)));
                vz = bit_or(bit_and(vtmp, vz), andnot(vtmp, va2));
                vtmp = cmpgt(vz, vb2);
                
                vd = bit_or(bit_and(vtmp, vd), andnot(vtmp, set1(4)));
                vz = bit_or(bit_and(vtmp, vz), andnot(vtmp, vb2));
                vtmp = cmplt(vsc_mch_, vz);
                vz = bit_or(bit_and(vtmp, vsc_mch_), andnot(vtmp, vz));

                __dp_code_block2;

                vtmp = cmpgt(vzero_, va);
                store(&x[t],  sub(andnot(vtmp, va),  vqe_));

                vd = bit_or(vd, andnot(vtmp, set1(0x08))); // d = a > 0? 1<<3 : 0
                vtmp = cmpgt(vzero_, vb);
                store(&y[t],  sub(andnot(vtmp, vb),  vqe_));
                
                vd = bit_or(vd, andnot(vtmp, set1(0x10))); // d = b > 0? 1<<4 : 0
                vtmp = cmpgt(vzero_, va2);
                store(&x2[t], sub(andnot(vtmp, va2), vqe2_));
                
                vd = bit_or(vd, andnot(vtmp, set1(0x20))); // d = a > 0? 1<<5 : 0
                vtmp = cmpgt(vzero_, vb2);
                store(&y2[t], sub(andnot(vtmp, vb2), vqe2_));
                
                vd = bit_or(vd, andnot(vtmp, set1(0x40))); // d = b > 0? 1<<6 : 0
                store(&pr[t], vd);
            }
        }
        if (!approx_max) { // find the exact max with a 32-bit score array
            int32_t max_H, max_t;
            // compute H[], max_H and max_t
            if (r > 0) {
                int32_t HH[4], tt[4], en1 = st0 + (en0 - st0) / 4 * 4, i;
                vec vmax_H_, vmax_t_;
                max_H = H[en0] = en0 > 0? H[en0-1] + u8[en0] : H[en0] + v8[en0]; // special casing the last element
                max_t = en0;

                vmax_H_ = set32(max_H);
                vmax_t_ = set32(max_t);

                for (t = st0; t < en1; t += 4) { // this implements: H[t]+=v8[t]-qe; if(H[t]>max_H) max_H=H[t],max_t=t;
                    vec vH1, vtmp, vt_;
                    vH1 = load(&H[t]);

                    vt_ = setr32(v8[t], v8[t+1], v8[t+2], v8[t+3]);
                    vH1 = add32(vH1, vt_);
                    store(&H[t], vH1);

                    vt_ = set32(t);
                    vtmp = cmpgt32(vH1, vmax_H_);

                    vmax_H_ = bit_or(bit_and(vtmp, vH1), andnot(vtmp, vmax_H_));
                    vmax_t_ = bit_or(bit_and(vtmp, vt_), andnot(vtmp, vmax_t_));

                }

                store(HH, vmax_H_);
                store(tt, vmax_t_);

                for (i = 0; i < 4; ++i)
                    if (max_H < HH[i]) max_H = HH[i], max_t = tt[i] + i;
                for (; t < en0; ++t) { // for the rest of values that haven't been computed with SSE
                    H[t] += (int32_t)v8[t];
                    if (H[t] > max_H)
                        max_H = H[t], max_t = t;
                }
            } else H[0] = v8[0] - qe, max_H = H[0], max_t = 0; // special casing r==0
            // update ez
            if (en0 == tlen - 1 && H[en0] > ez->mte)
                ez->mte = H[en0], ez->mte_q = r - en;
            if (r - st0 == qlen - 1 && H[st0] > ez->mqe)
                ez->mqe = H[st0], ez->mqe_t = st0;
            if (ksw_apply_zdrop(ez, 1, max_H, r, max_t, zdrop, e2)) break;
            if (r == qlen + tlen - 2 && en0 == tlen - 1)
                ez->score = H[tlen - 1];
        } else { // find approximate max; Z-drop might be inaccurate, too.
            if (r > 0) {
                if (last_H0_t >= st0 && last_H0_t <= en0 && last_H0_t + 1 >= st0 && last_H0_t + 1 <= en0) {
                    int32_t d0 = v8[last_H0_t];
                    int32_t d1 = u8[last_H0_t + 1];
                    if (d0 > d1) H0 += d0;
                    else H0 += d1, ++last_H0_t;
                } else if (last_H0_t >= st0 && last_H0_t <= en0) {
                    H0 += v8[last_H0_t];
                } else {
                    ++last_H0_t, H0 += u8[last_H0_t];
                }
            } else H0 = v8[0] - qe, last_H0_t = 0;
            if ((flag & KSW_EZ_APPROX_DROP) && ksw_apply_zdrop(ez, 1, H0, r, last_H0_t, zdrop, e2)) break;
            if (r == qlen + tlen - 2 && en0 == tlen - 1)
                ez->score = H0;
        }
        last_st = st, last_en = en;
        //for (t = st0; t <= en0; ++t) printf("(%d,%d)\t(%d,%d,%d,%d)\t%d\n", r, t, ((int8_t*)u)[t], ((int8_t*)v)[t], ((int8_t*)x)[t], ((int8_t*)y)[t], H[t]); // for debugging
    }
    kfree(km, mem);
    if (!approx_max) kfree(km, H);
    if (with_cigar) { // backtrack
        int rev_cigar = !!(flag & KSW_EZ_REV_CIGAR);
        if (!ez->zdropped && !(flag&KSW_EZ_EXTZ_ONLY)) {
            ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*16, tlen-1, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
        } else if (!ez->zdropped && (flag&KSW_EZ_EXTZ_ONLY) && ez->mqe + end_bonus > (int)ez->max) {
            ez->reach_end = 1;
            ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*16, ez->mqe_t, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
        } else if (ez->max_t >= 0 && ez->max_q >= 0) {
            ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*16, ez->max_t, ez->max_q, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
        }
        if (flag & KSW_EZ_EQX) {
            int32_t nc0 = ez->n_cigar;
            uint32_t *ci0;
            ci0 = (uint32_t*)kmalloc(km, nc0 * sizeof(uint32_t));
            memcpy(ci0, ez->cigar, nc0 * sizeof(uint32_t));
            ksw_cigar2eqx(km, query, target, nc0, ci0, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
            kfree(km, ci0);
        }
        kfree(km, mem2); kfree(km, off);
    }
}

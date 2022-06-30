#ifndef _SCHEME_H
#define _SCHEME_H 1

#include <math.h>
#include "util.h"
#include "basic.h"

class Scheme {
public:

static real_t ffvc_muscl(real_t uc0, real_t ue1, real_t ue2, real_t un1, real_t un2, real_t ut1, real_t ut2, real_t uw1, real_t uw2, real_t us1, real_t us2, real_t ub1, real_t ub2, real_t ufe, real_t vfn, real_t wft, real_t ufw, real_t vfs, real_t wfb, real_t det, int epe, int epn, int ept, int epw, int eps, int epb, real_t kappa, real_t beta, real_t m1, real_t m2) {
    real_t d1, d2, d3, d4;
    real_t s1, s2, s3, s4;
    real_t g1, g2, g3, g4, g5, g6;
    real_t u11, u10, u01, u00;
    real_t f1, f0;
    real_t adv;

    d4 = ue2 - ue1;
    d3 = ue1 - uc0;
    d2 = uc0 - uw1;
    d1 = uw1 - uw2;
    s4 = copysign(1.0, d4);
    s3 = copysign(1.0, d3);
    s2 = copysign(1.0, d2);
    s1 = copysign(1.0, d1);
    g6  = s4 * Util::max(0.0, Util::min(fabs(d4), s4 * beta * d3));
    g5  = s3 * Util::max(0.0, Util::min(fabs(d3), s3 * beta * d4));
    g4  = s3 * Util::max(0.0, Util::min(fabs(d3), s3 * beta * d2));
    g3  = s2 * Util::max(0.0, Util::min(fabs(d2), s2 * beta * d3));
    g2  = s2 * Util::max(0.0, Util::min(fabs(d2), s2 * beta * d1));
    g1  = s1 * Util::max(0.0, Util::min(fabs(d1), s1 * beta * d2));
    u11 = ue1 - 0.25 * (m1 * g6 + m2 * g5) * epe;
    u10 = uc0 + 0.25 * (m1 * g3 + m2 * g4);
    u01 = uc0 - 0.25 * (m1 * g4 + m2 * g3);
    u00 = uw1 + 0.25 * (m1 * g1 + m2 * g2) * epw;
    f1  = 0.5 * (ufe * (u11 + u10) - fabs(ufe) * (u11 - u10));
    f0  = 0.5 * (ufw * (u01 + u00) - fabs(ufw) * (u01 - u00));
    adv = (f1 - f0) / det;

    d4  = un2 - un1;
    d3  = un1 - uc0;
    d2  = uc0 - us1;
    d1  = us1 - us2;
    s4  = copysign(1.0, d4);
    s3  = copysign(1.0, d3);
    s2  = copysign(1.0, d2);
    s1  = copysign(1.0, d1);
    g6  = s4 * Util::max(0.0, Util::min(fabs(d4), s4 * beta * d3));
    g5  = s3 * Util::max(0.0, Util::min(fabs(d3), s3 * beta * d4));
    g4  = s3 * Util::max(0.0, Util::min(fabs(d3), s3 * beta * d2));
    g3  = s2 * Util::max(0.0, Util::min(fabs(d2), s2 * beta * d3));
    g2  = s2 * Util::max(0.0, Util::min(fabs(d2), s2 * beta * d1));
    g1  = s1 * Util::max(0.0, Util::min(fabs(d1), s1 * beta * d2));
    u11 = un1 - 0.25 * (m1 * g6 + m2 * g5) * epn;
    u10 = uc0 + 0.25 * (m1 * g3 + m2 * g4);
    u01 = uc0 - 0.25 * (m1 * g4 + m2 * g3);
    u00 = us1 + 0.25 * (m1 * g1 + m2 * g2) * eps;
    f1  = 0.5 * (vfn * (u11 + u10) - fabs(vfn) * (u11 - u10));
    f0  = 0.5 * (vfs * (u01 + u00) - fabs(vfs) * (u01 - u00));
    adv+= (f1 - f0) / det;

    d4  = ut2 - ut1;
    d3  = ut1 - uc0;
    d2  = uc0 - ub1;
    d1  = ub1 - ub2;
    s4  = copysign(1.0, d4);
    s3  = copysign(1.0, d3);
    s2  = copysign(1.0, d2);
    s1  = copysign(1.0, d1);
    g6  = s4 * Util::max(0.0, Util::min(fabs(d4), s4 * beta * d3));
    g5  = s3 * Util::max(0.0, Util::min(fabs(d3), s3 * beta * d4));
    g4  = s3 * Util::max(0.0, Util::min(fabs(d3), s3 * beta * d2));
    g3  = s2 * Util::max(0.0, Util::min(fabs(d2), s2 * beta * d3));
    g2  = s2 * Util::max(0.0, Util::min(fabs(d2), s2 * beta * d1));
    g1  = s1 * Util::max(0.0, Util::min(fabs(d1), s1 * beta * d2));
    u11 = ut1 - 0.25 * (m1 * g6 + m2 * g5) * ept;
    u10 = uc0 + 0.25 * (m1 * g3 + m2 * g4);
    u01 = uc0 - 0.25 * (m1 * g4 + m2 * g3);
    u00 = ub1 + 0.25 * (m1 * g1 + m2 * g2) * epb;
    f1  = 0.5 * (wft * (u11 + u10) - fabs(wft) * (u11 - u10));
    f0  = 0.5 * (wfb * (u01 + u00) - fabs(wfb) * (u01 - u00));
    adv+= (f1 - f0) / det;

    return adv;
}

static real_t viscosity(unsigned w13, unsigned w23, unsigned w33, unsigned w12, unsigned w22, unsigned w32, real_t mag, real_t uc0, real_t ue1, real_t un1, real_t ut1, real_t uw1, real_t us1, real_t ub1, real_t de1, real_t dn1, real_t dt1, real_t dw1, real_t ds1, real_t db1, real_t nc0, real_t ne1, real_t nn1, real_t nt1, real_t nw1, real_t ns1, real_t nb1, real_t det, real_t g1c, real_t g2c, real_t g3c, real_t g1e, real_t g2n, real_t g3t, real_t g1w, real_t g2s, real_t g3b, real_t ri) {
    real_t ffe, ffn, fft, ffw, ffs, ffb;
    ffe = (ri + 0.5 * (nc0 + ne1)) * 0.5 * (g1e + g1c) * (ue1 - uc0);
    ffn = (ri + 0.5 * (nc0 + nn1)) * 0.5 * (g2n + g2c) * (un1 - uc0);
    fft = (ri + 0.5 * (nc0 + nt1)) * 0.5 * (g3t + g3c) * (ut1 - uc0);
    ffw = (ri + 0.5 * (nc0 + nw1)) * 0.5 * (g1c + g1w) * (uc0 - uw1);
    ffs = (ri + 0.5 * (nc0 + ns1)) * 0.5 * (g2c + g2s) * (uc0 - us1);
    ffb = (ri + 0.5 * (nc0 + nb1)) * 0.5 * (g3c + g3b) * (uc0 - ub1);

    return (ffe + ffn + fft - ffw - ffs - ffb) / det;
}

};

#endif

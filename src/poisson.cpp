#include <math.h>
#include "../include/poisson.h"
#include "../include/flag.h"
#include "../include/boundary.h"

void Poisson::sor(scalar_field<real_t> &p, scalar_field<real_t> &div, Dom &dom) {
    real_t                  err = 0;
    int                     cnt = 0;
    Ctrl                   &c   = dom.c;
    scalar_field<unsigned> &f   = dom.f;
    vector_field<real_t>   &g   = dom.g;
    scalar_field<real_t>   &ja  = dom.ja;

    c.poisson.it = 0;
    do {
        #pragma acc kernels loop independent collapse(3) reduction(+:err, cnt) present(p, div, f, g, ja, c, dom) copy(err, cnt)
        for (int i = GUIDE; i < dom.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < dom.size[1] - GUIDE; j ++) {
                for (int k = GUIDE; k < dom.size[2] - GUIDE; k++) {
                    if (Util::ibsee(f.m[id3(i,j,k,f.size)], Flag::Active, Util::Mask1) && (i+j+k) % 2 == 0) {
                        unsigned f12, f13;
                        unsigned f22, f23;
                        unsigned f32, f33;
                        unsigned m12, m13;
                        unsigned m22, m23;
                        unsigned m32, m33;
                        unsigned b12, b13;
                        unsigned b22, b23;
                        unsigned b32, b33;
                        real_t   ac0;
                        real_t   ae1, an1, at1;
                        real_t   aw1, as1, ab1;
                        real_t   pc0;
                        real_t   pe1, pn1, pt1;
                        real_t   pw1, ps1, pb1;
                        real_t   g1c, g2c, g3c;
                        real_t   g1e, g2n, g3t;
                        real_t   g1w, g2s, g3b;
                        real_t   det;
                        real_t   ref, dis;
                        real_t   rhs, rc0;

                        rhs = div.m[id3(i  ,j  ,k  ,  div.size)] / c.time.dt;
                        det =  ja.m[id3(i  ,j  ,k  ,   ja.size)];
                        g1c =   g.m[id4(i  ,j  ,k  ,0,  g.size)];
                        g2c =   g.m[id4(i  ,j  ,k  ,1,  g.size)];
                        g3c =   g.m[id4(i  ,j  ,k  ,2,  g.size)];
                        g1e =   g.m[id4(i+1,j  ,k  ,0,  g.size)];
                        g2n =   g.m[id4(i  ,j+1,k  ,1,  g.size)];
                        g3t =   g.m[id4(i  ,j  ,k+1,2,  g.size)];
                        g1w =   g.m[id4(i-1,j  ,k  ,0,  g.size)];
                        g2s =   g.m[id4(i  ,j-1,k  ,1,  g.size)];
                        g3b =   g.m[id4(i  ,j  ,k-1,2,  g.size)];
                        pc0 =   p.m[id3(i  ,j  ,k  ,    p.size)];
                        pe1 =   p.m[id3(i+1,j  ,k  ,    p.size)];
                        pn1 =   p.m[id3(i  ,j+1,k  ,    p.size)];
                        pt1 =   p.m[id3(i  ,j  ,k+1,    p.size)];
                        pw1 =   p.m[id3(i-1,j  ,k  ,    p.size)];
                        ps1 =   p.m[id3(i  ,j-1,k  ,    p.size)];
                        pb1 =   p.m[id3(i  ,j  ,k-1,    p.size)];
                        f13 =   f.m[id3(i  ,j  ,k  ,    f.size)];
                        f23 =   f.m[id3(i  ,j  ,k  ,    f.size)];
                        f33 =   f.m[id3(i  ,j  ,k  ,    f.size)];
                        f12 =   f.m[id3(i-1,j  ,k  ,    f.size)];
                        f22 =   f.m[id3(i  ,j-1,k  ,    f.size)];
                        f32 =   f.m[id3(i  ,j  ,k-1,    f.size)];
                        m12 = Util::ibsee(f12, Flag::Me, Util::Mask1);
                        m13 = Util::ibsee(f13, Flag::Me, Util::Mask1);
                        m22 = Util::ibsee(f22, Flag::Mn, Util::Mask1);
                        m23 = Util::ibsee(f23, Flag::Mn, Util::Mask1);
                        m32 = Util::ibsee(f32, Flag::Mt, Util::Mask1);
                        m33 = Util::ibsee(f33, Flag::Mt, Util::Mask1);
                        f12 = Util::ibsee(f12, Flag::Fe, Util::Mask8);
                        f13 = Util::ibsee(f13, Flag::Fe, Util::Mask8);
                        f22 = Util::ibsee(f22, Flag::Fn, Util::Mask8);
                        f23 = Util::ibsee(f23, Flag::Fn, Util::Mask8);
                        f32 = Util::ibsee(f32, Flag::Ft, Util::Mask8);
                        f33 = Util::ibsee(f33, Flag::Ft, Util::Mask8);
                        b12 = Util::ibsee(p.bflag[f12], 0, Util::Mask8);
                        b13 = Util::ibsee(p.bflag[f13], 0, Util::Mask8);
                        b22 = Util::ibsee(p.bflag[f22], 0, Util::Mask8);
                        b23 = Util::ibsee(p.bflag[f23], 0, Util::Mask8);
                        b32 = Util::ibsee(p.bflag[f32], 0, Util::Mask8);
                        b33 = Util::ibsee(p.bflag[f33], 0, Util::Mask8);
                        if (b13) {
                            BB::pre(m13, pc0, pe1, ref, dis);
                            pe1 = 2 * BB::eva(b13, ref, dis, p.b[f13]) - pc0;
                        }
                        if (b23) {
                            BB::pre(m23, pc0, pn1, ref, dis);
                            pn1 = 2 * BB::eva(b23, ref, dis, p.b[f23]) - pc0;
                        }
                        if (b33) {
                            BB::pre(m33, pc0, pt1, ref, dis);
                            pt1 = 2 * BB::eva(b33, ref, dis, p.b[f33]) - pc0;
                        }
                        if (b12) {
                            BB::pre(m12, pw1, pc0, ref, dis);
                            pw1 = 2 * BB::eva(b12, ref, dis, p.b[f12]) - pc0;
                        }
                        if (b22) {
                            BB::pre(m22, ps1, pc0, ref, dis);
                            ps1 = 2 * BB::eva(b22, ref, dis, p.b[f22]) - pc0;
                        }
                        if (b32) {
                            BB::pre(m32, pb1, pc0, ref, dis);
                            pb1 = 2 * BB::eva(b32, ref, dis, p.b[f32]) - pc0;
                        }

                        real_t ddet = 2 * det;
                        ae1 = (g1e + g1c) / ddet;
                        an1 = (g2n + g2c) / ddet;
                        at1 = (g3t + g3c) / ddet;
                        aw1 = (g1w + g1c) / ddet;
                        as1 = (g2s + g2c) / ddet;
                        ab1 = (g3b + g3c) / ddet;
                        ac0 = - (ae1 + an1 + at1 + aw1 + as1 + ab1);
                        rc0 = (rhs - ae1 * pe1 - an1 * pn1 - at1 * pt1 - aw1 * pw1 - as1 * ps1 - ab1 * pb1) / ac0 - pc0;
                        err = err + rc0 * rc0;

                        p.m[id3(i,j,k,p.size)] = pc0 + c.poisson.omega * rc0;
                        cnt += 1;
                    }
                }
            }
        }
        #pragma acc kernels loop independent collapse(3) reduction(+:err, cnt) present(p, div, f, g, ja, c, dom) copy(err, cnt)
        for (int i = GUIDE; i < dom.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < dom.size[1] - GUIDE; j ++) {
                for (int k = GUIDE; k < dom.size[2] - GUIDE; k++) {
                    if (Util::ibsee(f.m[id3(i,j,k,f.size)], Flag::Active, Util::Mask1) && (i+j+k) % 2 == 1) {
                        unsigned f12, f13;
                        unsigned f22, f23;
                        unsigned f32, f33;
                        unsigned m12, m13;
                        unsigned m22, m23;
                        unsigned m32, m33;
                        unsigned b12, b13;
                        unsigned b22, b23;
                        unsigned b32, b33;
                        real_t   ac0;
                        real_t   ae1, an1, at1;
                        real_t   aw1, as1, ab1;
                        real_t   pc0;
                        real_t   pe1, pn1, pt1;
                        real_t   pw1, ps1, pb1;
                        real_t   g1c, g2c, g3c;
                        real_t   g1e, g2n, g3t;
                        real_t   g1w, g2s, g3b;
                        real_t   det;
                        real_t   ref, dis;
                        real_t   rhs, rc0;

                        rhs = div.m[id3(i  ,j  ,k  ,  div.size)] / c.time.dt;
                        det =  ja.m[id3(i  ,j  ,k  ,   ja.size)];
                        g1c =   g.m[id4(i  ,j  ,k  ,0,  g.size)];
                        g2c =   g.m[id4(i  ,j  ,k  ,1,  g.size)];
                        g3c =   g.m[id4(i  ,j  ,k  ,2,  g.size)];
                        g1e =   g.m[id4(i+1,j  ,k  ,0,  g.size)];
                        g2n =   g.m[id4(i  ,j+1,k  ,1,  g.size)];
                        g3t =   g.m[id4(i  ,j  ,k+1,2,  g.size)];
                        g1w =   g.m[id4(i-1,j  ,k  ,0,  g.size)];
                        g2s =   g.m[id4(i  ,j-1,k  ,1,  g.size)];
                        g3b =   g.m[id4(i  ,j  ,k-1,2,  g.size)];
                        pc0 =   p.m[id3(i  ,j  ,k  ,    p.size)];
                        pe1 =   p.m[id3(i+1,j  ,k  ,    p.size)];
                        pn1 =   p.m[id3(i  ,j+1,k  ,    p.size)];
                        pt1 =   p.m[id3(i  ,j  ,k+1,    p.size)];
                        pw1 =   p.m[id3(i-1,j  ,k  ,    p.size)];
                        ps1 =   p.m[id3(i  ,j-1,k  ,    p.size)];
                        pb1 =   p.m[id3(i  ,j  ,k-1,    p.size)];
                        f13 =   f.m[id3(i  ,j  ,k  ,    f.size)];
                        f23 =   f.m[id3(i  ,j  ,k  ,    f.size)];
                        f33 =   f.m[id3(i  ,j  ,k  ,    f.size)];
                        f12 =   f.m[id3(i-1,j  ,k  ,    f.size)];
                        f22 =   f.m[id3(i  ,j-1,k  ,    f.size)];
                        f32 =   f.m[id3(i  ,j  ,k-1,    f.size)];
                        m12 = Util::ibsee(f12, Flag::Me, Util::Mask1);
                        m13 = Util::ibsee(f13, Flag::Me, Util::Mask1);
                        m22 = Util::ibsee(f22, Flag::Mn, Util::Mask1);
                        m23 = Util::ibsee(f23, Flag::Mn, Util::Mask1);
                        m32 = Util::ibsee(f32, Flag::Mt, Util::Mask1);
                        m33 = Util::ibsee(f33, Flag::Mt, Util::Mask1);
                        f12 = Util::ibsee(f12, Flag::Fe, Util::Mask8);
                        f13 = Util::ibsee(f13, Flag::Fe, Util::Mask8);
                        f22 = Util::ibsee(f22, Flag::Fn, Util::Mask8);
                        f23 = Util::ibsee(f23, Flag::Fn, Util::Mask8);
                        f32 = Util::ibsee(f32, Flag::Ft, Util::Mask8);
                        f33 = Util::ibsee(f33, Flag::Ft, Util::Mask8);
                        b12 = Util::ibsee(p.bflag[f12], 0, Util::Mask8);
                        b13 = Util::ibsee(p.bflag[f13], 0, Util::Mask8);
                        b22 = Util::ibsee(p.bflag[f22], 0, Util::Mask8);
                        b23 = Util::ibsee(p.bflag[f23], 0, Util::Mask8);
                        b32 = Util::ibsee(p.bflag[f32], 0, Util::Mask8);
                        b33 = Util::ibsee(p.bflag[f33], 0, Util::Mask8);
                        if (b13) {
                            BB::pre(m13, pc0, pe1, ref, dis);
                            pe1 = 2 * BB::eva(b13, ref, dis, p.b[f13]) - pc0;
                        }
                        if (b23) {
                            BB::pre(m23, pc0, pn1, ref, dis);
                            pn1 = 2 * BB::eva(b23, ref, dis, p.b[f23]) - pc0;
                        }
                        if (b33) {
                            BB::pre(m33, pc0, pt1, ref, dis);
                            pt1 = 2 * BB::eva(b33, ref, dis, p.b[f33]) - pc0;
                        }
                        if (b12) {
                            BB::pre(m12, pw1, pc0, ref, dis);
                            pw1 = 2 * BB::eva(b12, ref, dis, p.b[f12]) - pc0;
                        }
                        if (b22) {
                            BB::pre(m22, ps1, pc0, ref, dis);
                            ps1 = 2 * BB::eva(b22, ref, dis, p.b[f22]) - pc0;
                        }
                        if (b32) {
                            BB::pre(m32, pb1, pc0, ref, dis);
                            pb1 = 2 * BB::eva(b32, ref, dis, p.b[f32]) - pc0;
                        }

                        real_t ddet = 2 * det;
                        ae1 = (g1e + g1c) / ddet;
                        an1 = (g2n + g2c) / ddet;
                        at1 = (g3t + g3c) / ddet;
                        aw1 = (g1w + g1c) / ddet;
                        as1 = (g2s + g2c) / ddet;
                        ab1 = (g3b + g3c) / ddet;
                        ac0 = - (ae1 + an1 + at1 + aw1 + as1 + ab1);
                        rc0 = (rhs - ae1 * pe1 - an1 * pn1 - at1 * pt1 - aw1 * pw1 - as1 * ps1 - ab1 * pb1) / ac0 - pc0;
                        err = err + rc0 * rc0;

                        p.m[id3(i,j,k,p.size)] = pc0 + c.poisson.omega * rc0;
                        cnt += 1;
                    }
                }
            }
        }
        c.poisson.res = sqrt(err / cnt);

        BB::scalar_outer(p, dom);
    } while (++c.poisson.it < c.poisson.maxit && c.poisson.res > c.poisson.epsilon);
}

void Poisson::jacobi(scalar_field<real_t> &p, scalar_field<real_t> &div, Dom &dom) {
    real_t                  err = 0;
    int                     cnt = 0;
    Ctrl                   &c   = dom.c;
    scalar_field<unsigned> &f   = dom.f;
    vector_field<real_t>   &g   = dom.g;
    scalar_field<real_t>   &ja  = dom.ja;

    c.poisson.it = 0;
    do {
        #pragma acc kernels loop independent collapse(3) present(this[0:1], p, pd)
        for (int i = 0; i < p.size[0]; i ++) {
            for (int j = 0; j < p.size[1]; j ++) {
                for (int k = 0; k < p.size[2]; k++) {
                    pd.m[id3(i,j,k,pd.size)] = p.m[id3(i,j,k,p.size)];
                }
            }
        }
        #pragma acc kernels loop independent collapse(3) reduction(+:err, cnt) present(this[0:1], p, pd, div, f, g, ja, c, dom) copy(err, cnt)
        for (int i = GUIDE; i < dom.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < dom.size[1] - GUIDE; j ++) {
                for (int k = GUIDE; k < dom.size[2] - GUIDE; k++) {
                    if (Util::ibsee(f.m[id3(i,j,k,f.size)], Flag::Active, Util::Mask1)) {
                        unsigned f12, f13;
                        unsigned f22, f23;
                        unsigned f32, f33;
                        unsigned m12, m13;
                        unsigned m22, m23;
                        unsigned m32, m33;
                        unsigned b12, b13;
                        unsigned b22, b23;
                        unsigned b32, b33;
                        real_t   ac0;
                        real_t   ae1, an1, at1;
                        real_t   aw1, as1, ab1;
                        real_t   pc0;
                        real_t   pe1, pn1, pt1;
                        real_t   pw1, ps1, pb1;
                        real_t   g1c, g2c, g3c;
                        real_t   g1e, g2n, g3t;
                        real_t   g1w, g2s, g3b;
                        real_t   det;
                        real_t   ref, dis;
                        real_t   rhs, rc0;
    
                        rhs = div.m[id3(i  ,j  ,k  ,  div.size)] / c.time.dt;
                        det =  ja.m[id3(i  ,j  ,k  ,   ja.size)];
                        g1c =   g.m[id4(i  ,j  ,k  ,0,  g.size)];
                        g2c =   g.m[id4(i  ,j  ,k  ,1,  g.size)];
                        g3c =   g.m[id4(i  ,j  ,k  ,2,  g.size)];
                        g1e =   g.m[id4(i+1,j  ,k  ,0,  g.size)];
                        g2n =   g.m[id4(i  ,j+1,k  ,1,  g.size)];
                        g3t =   g.m[id4(i  ,j  ,k+1,2,  g.size)];
                        g1w =   g.m[id4(i-1,j  ,k  ,0,  g.size)];
                        g2s =   g.m[id4(i  ,j-1,k  ,1,  g.size)];
                        g3b =   g.m[id4(i  ,j  ,k-1,2,  g.size)];
                        pc0 =  pd.m[id3(i  ,j  ,k  ,   pd.size)];
                        pe1 =  pd.m[id3(i+1,j  ,k  ,   pd.size)];
                        pn1 =  pd.m[id3(i  ,j+1,k  ,   pd.size)];
                        pt1 =  pd.m[id3(i  ,j  ,k+1,   pd.size)];
                        pw1 =  pd.m[id3(i-1,j  ,k  ,   pd.size)];
                        ps1 =  pd.m[id3(i  ,j-1,k  ,   pd.size)];
                        pb1 =  pd.m[id3(i  ,j  ,k-1,   pd.size)];
                        f13 =   f.m[id3(i  ,j  ,k  ,    f.size)];
                        f23 =   f.m[id3(i  ,j  ,k  ,    f.size)];
                        f33 =   f.m[id3(i  ,j  ,k  ,    f.size)];
                        f12 =   f.m[id3(i-1,j  ,k  ,    f.size)];
                        f22 =   f.m[id3(i  ,j-1,k  ,    f.size)];
                        f32 =   f.m[id3(i  ,j  ,k-1,    f.size)];
                        m12 = Util::ibsee(f12, Flag::Me, Util::Mask1);
                        m13 = Util::ibsee(f13, Flag::Me, Util::Mask1);
                        m22 = Util::ibsee(f22, Flag::Mn, Util::Mask1);
                        m23 = Util::ibsee(f23, Flag::Mn, Util::Mask1);
                        m32 = Util::ibsee(f32, Flag::Mt, Util::Mask1);
                        m33 = Util::ibsee(f33, Flag::Mt, Util::Mask1);
                        f12 = Util::ibsee(f12, Flag::Fe, Util::Mask8);
                        f13 = Util::ibsee(f13, Flag::Fe, Util::Mask8);
                        f22 = Util::ibsee(f22, Flag::Fn, Util::Mask8);
                        f23 = Util::ibsee(f23, Flag::Fn, Util::Mask8);
                        f32 = Util::ibsee(f32, Flag::Ft, Util::Mask8);
                        f33 = Util::ibsee(f33, Flag::Ft, Util::Mask8);
                        b12 = Util::ibsee(pd.bflag[f12], 0, Util::Mask8);
                        b13 = Util::ibsee(pd.bflag[f13], 0, Util::Mask8);
                        b22 = Util::ibsee(pd.bflag[f22], 0, Util::Mask8);
                        b23 = Util::ibsee(pd.bflag[f23], 0, Util::Mask8);
                        b32 = Util::ibsee(pd.bflag[f32], 0, Util::Mask8);
                        b33 = Util::ibsee(pd.bflag[f33], 0, Util::Mask8);
                        if (b13) {
                            BB::pre(m13, pc0, pe1, ref, dis);
                            pe1 = 2 * BB::eva(b13, ref, dis, pd.b[f13]) - pc0;
                        }
                        if (b23) {
                            BB::pre(m23, pc0, pn1, ref, dis);
                            pn1 = 2 * BB::eva(b23, ref, dis, pd.b[f23]) - pc0;
                        }
                        if (b33) {
                            BB::pre(m33, pc0, pt1, ref, dis);
                            pt1 = 2 * BB::eva(b33, ref, dis, pd.b[f33]) - pc0;
                        }
                        if (b12) {
                            BB::pre(m12, pw1, pc0, ref, dis);
                            pw1 = 2 * BB::eva(b12, ref, dis, pd.b[f12]) - pc0;
                        }
                        if (b22) {
                            BB::pre(m22, ps1, pc0, ref, dis);
                            ps1 = 2 * BB::eva(b22, ref, dis, pd.b[f22]) - pc0;
                        }
                        if (b32) {
                            BB::pre(m32, pb1, pc0, ref, dis);
                            pb1 = 2 * BB::eva(b32, ref, dis, pd.b[f32]) - pc0;
                        }
    
                        real_t ddet = 2 * det;
                        ae1 = (g1e + g1c) / ddet;
                        an1 = (g2n + g2c) / ddet;
                        at1 = (g3t + g3c) / ddet;
                        aw1 = (g1w + g1c) / ddet;
                        as1 = (g2s + g2c) / ddet;
                        ab1 = (g3b + g3c) / ddet;
                        ac0 = - (ae1 + an1 + at1 + aw1 + as1 + ab1);
                        rc0 = (rhs - ae1 * pe1 - an1 * pn1 - at1 * pt1 - aw1 * pw1 - as1 * ps1 - ab1 * pb1) / ac0 - pc0;
                        err = err + rc0 * rc0;
    
                        p.m[id3(i,j,k,p.size)] = pc0 + rc0;
                        cnt += 1;
                    }
                }
            }
        }
        c.poisson.res = sqrt(err / cnt);

        BB::scalar_outer(p, dom);
    } while (++c.poisson.it < c.poisson.maxit && c.poisson.res > c.poisson.epsilon);
}

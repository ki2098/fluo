#include <math.h>
#include <stdio.h>
#include <float.h>
#include "../include/mmac.h"
#include "../include/flag.h"
#include "../include/boundary.h"

void MMAC::Poisson::jacobi_core(scalar_field<real_t> &p, scalar_field<real_t> &psi, Dom &dom, int maxit, real_t torelence, bool ignore, int &it, real_t &err) {
    scalar_field<unsigned> &f   = dom.f;
    vector_field<real_t>   &g   = dom.g;
    scalar_field<real_t>   &ja  = dom.ja;

    it = 0;
    do {
        real_t _err = 0;
        int    _cnt = 0;
        pd.copy(p);
        #pragma acc kernels loop independent collapse(3) reduction(+:_err, _cnt) present(this[0:1], p, pd, psi, f, g, ja, dom) copy(_err, _cnt)
        for (int i = GUIDE; i < dom.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < dom.size[1] - GUIDE; j ++) {
                for (int k = GUIDE; k < dom.size[2] - GUIDE; k++) {
                    if (Util::ibsee(f.m[id3(i,j,k,f.size)], Cell::Active, Util::Mask1)) {
                        unsigned f12, f13;
                        unsigned f22, f23;
                        unsigned f32, f33;
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
                        real_t   rhs, rc0;
    
                        rhs = psi.m[id3(i  ,j  ,k  ,  psi.size)];
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
                        f12 = Util::ibsee(f12, Cell::Fe, Util::Mask8);
                        f13 = Util::ibsee(f13, Cell::Fe, Util::Mask8);
                        f22 = Util::ibsee(f22, Cell::Fn, Util::Mask8);
                        f23 = Util::ibsee(f23, Cell::Fn, Util::Mask8);
                        f32 = Util::ibsee(f32, Cell::Ft, Util::Mask8);
                        f33 = Util::ibsee(f33, Cell::Ft, Util::Mask8);
                        b12 = Util::ibsee(dom.p.bflag[f12], 0, Util::Mask8);
                        b13 = Util::ibsee(dom.p.bflag[f13], 0, Util::Mask8);
                        b22 = Util::ibsee(dom.p.bflag[f22], 0, Util::Mask8);
                        b23 = Util::ibsee(dom.p.bflag[f23], 0, Util::Mask8);
                        b32 = Util::ibsee(dom.p.bflag[f32], 0, Util::Mask8);
                        b33 = Util::ibsee(dom.p.bflag[f33], 0, Util::Mask8);
                        unsigned feN, feD, fnN, fnD, ftN, ftD;
                        unsigned fwN, fwD, fsN, fsD, fbN, fbD;
                        feN = Util::ibsee(b13, BB::neumann,   Util::Mask1);
                        feD = Util::ibsee(b13, BB::dirichlet, Util::Mask1);
                        fnN = Util::ibsee(b23, BB::neumann,   Util::Mask1);
                        fnD = Util::ibsee(b23, BB::dirichlet, Util::Mask1);
                        ftN = Util::ibsee(b33, BB::neumann,   Util::Mask1);
                        ftD = Util::ibsee(b33, BB::dirichlet, Util::Mask1);
                        fwN = Util::ibsee(b12, BB::neumann,   Util::Mask1);
                        fwD = Util::ibsee(b12, BB::dirichlet, Util::Mask1);
                        fsN = Util::ibsee(b22, BB::neumann,   Util::Mask1);
                        fsD = Util::ibsee(b22, BB::dirichlet, Util::Mask1);
                        fbN = Util::ibsee(b32, BB::neumann,   Util::Mask1);
                        fbD = Util::ibsee(b32, BB::dirichlet, Util::Mask1);
                        
                        real_t jge, jgw, jgn, jgs, jgt, jgb;
                        jge = 0.5 * (g1e + g1c);
                        jgn = 0.5 * (g2n + g2c);
                        jgt = 0.5 * (g3t + g3c);
                        jgw = 0.5 * (g1w + g1c);
                        jgs = 0.5 * (g2s + g2c);
                        jgb = 0.5 * (g3b + g3c);
                        
                        ae1 = (1 - feN) * (1 - feD) * (jge) / det;
                        an1 = (1 - fnN) * (1 - fnD) * (jgn) / det;
                        at1 = (1 - ftN) * (1 - ftD) * (jgt) / det;
                        aw1 = (1 - fwN) * (1 - fwD) * (jgw) / det;
                        as1 = (1 - fsN) * (1 - fsD) * (jgs) / det;
                        ab1 = (1 - fbN) * (1 - fbD) * (jgb) / det;
                        ac0 = - (ae1 + an1 + at1 + aw1 + as1 + ab1) - 2 * (jge * feD + jgw * fwD + jgn * fnD + jgs * fsD + jgt * ftD + jgb * fbD);
                        rc0 = (rhs - ae1 * pe1 - an1 * pn1 - at1 * pt1 - aw1 * pw1 - as1 * ps1 - ab1 * pb1) / ac0 - pc0;
                        
                        p.m[id3(i,j,k,p.size)] = pc0 + rc0;
                        _err += rc0 * rc0;
                        _cnt += 1;
                    }
                }
            }
        }
        err = sqrt(_err / _cnt);
        BB::scalar_outer(p, dom);
    } while (++it < maxit && (ignore || (!ignore && err > torelence)));
}

void MMAC::Poisson::sor_core(scalar_field<real_t> &p, scalar_field<real_t> &psi, Dom &dom, real_t omega, int maxit, real_t torelence, bool ignore, int &it, real_t &err) {
    scalar_field<unsigned> &f   = dom.f;
    vector_field<real_t>   &g   = dom.g;
    scalar_field<real_t>   &ja  = dom.ja;

    it = 0;
    do {
        real_t _err = 0;
        int    _cnt = 0;
        #pragma acc kernels loop independent collapse(3) reduction(+:_err, _cnt) present(p, psi, f, g, ja, dom) copy(_err, _cnt)
        for (int i = GUIDE; i < dom.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < dom.size[1] - GUIDE; j ++) {
                for (int k = GUIDE; k < dom.size[2] - GUIDE; k++) {
                    if (Util::ibsee(f.m[id3(i,j,k,f.size)], Cell::Active, Util::Mask1) && (i+j+k) % 2 == 0) {
                        unsigned f12, f13;
                        unsigned f22, f23;
                        unsigned f32, f33;
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
                        real_t   rhs, rc0;
    
                        rhs = psi.m[id3(i  ,j  ,k  ,  psi.size)];
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
                        f12 = Util::ibsee(f12, Cell::Fe, Util::Mask8);
                        f13 = Util::ibsee(f13, Cell::Fe, Util::Mask8);
                        f22 = Util::ibsee(f22, Cell::Fn, Util::Mask8);
                        f23 = Util::ibsee(f23, Cell::Fn, Util::Mask8);
                        f32 = Util::ibsee(f32, Cell::Ft, Util::Mask8);
                        f33 = Util::ibsee(f33, Cell::Ft, Util::Mask8);
                        b12 = Util::ibsee(dom.p.bflag[f12], 0, Util::Mask8);
                        b13 = Util::ibsee(dom.p.bflag[f13], 0, Util::Mask8);
                        b22 = Util::ibsee(dom.p.bflag[f22], 0, Util::Mask8);
                        b23 = Util::ibsee(dom.p.bflag[f23], 0, Util::Mask8);
                        b32 = Util::ibsee(dom.p.bflag[f32], 0, Util::Mask8);
                        b33 = Util::ibsee(dom.p.bflag[f33], 0, Util::Mask8);
                        unsigned feN, feD, fnN, fnD, ftN, ftD;
                        unsigned fwN, fwD, fsN, fsD, fbN, fbD;
                        feN = Util::ibsee(b13, BB::neumann,   Util::Mask1);
                        feD = Util::ibsee(b13, BB::dirichlet, Util::Mask1);
                        fnN = Util::ibsee(b23, BB::neumann,   Util::Mask1);
                        fnD = Util::ibsee(b23, BB::dirichlet, Util::Mask1);
                        ftN = Util::ibsee(b33, BB::neumann,   Util::Mask1);
                        ftD = Util::ibsee(b33, BB::dirichlet, Util::Mask1);
                        fwN = Util::ibsee(b12, BB::neumann,   Util::Mask1);
                        fwD = Util::ibsee(b12, BB::dirichlet, Util::Mask1);
                        fsN = Util::ibsee(b22, BB::neumann,   Util::Mask1);
                        fsD = Util::ibsee(b22, BB::dirichlet, Util::Mask1);
                        fbN = Util::ibsee(b32, BB::neumann,   Util::Mask1);
                        fbD = Util::ibsee(b32, BB::dirichlet, Util::Mask1);
                        
                        real_t jge, jgw, jgn, jgs, jgt, jgb;
                        jge = 0.5 * (g1e + g1c);
                        jgn = 0.5 * (g2n + g2c);
                        jgt = 0.5 * (g3t + g3c);
                        jgw = 0.5 * (g1w + g1c);
                        jgs = 0.5 * (g2s + g2c);
                        jgb = 0.5 * (g3b + g3c);
                        
                        ae1 = (1 - feN) * (1 - feD) * (jge) / det;
                        an1 = (1 - fnN) * (1 - fnD) * (jgn) / det;
                        at1 = (1 - ftN) * (1 - ftD) * (jgt) / det;
                        aw1 = (1 - fwN) * (1 - fwD) * (jgw) / det;
                        as1 = (1 - fsN) * (1 - fsD) * (jgs) / det;
                        ab1 = (1 - fbN) * (1 - fbD) * (jgb) / det;
                        ac0 = - (ae1 + an1 + at1 + aw1 + as1 + ab1) - 2 * (jge * feD + jgw * fwD + jgn * fnD + jgs * fsD + jgt * ftD + jgb * fbD);
                        rc0 = (rhs - ae1 * pe1 - an1 * pn1 - at1 * pt1 - aw1 * pw1 - as1 * ps1 - ab1 * pb1) / ac0 - pc0;
                        
                        p.m[id3(i,j,k,p.size)] = pc0 + omega * rc0;
                        _err += rc0 * rc0;
                        _cnt += 1;
                    }
                }
            }
        }
        #pragma acc kernels loop independent collapse(3) reduction(+:_err, _cnt) present(p, psi, f, g, ja, dom) copy(_err, _cnt)
        for (int i = GUIDE; i < dom.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < dom.size[1] - GUIDE; j ++) {
                for (int k = GUIDE; k < dom.size[2] - GUIDE; k++) {
                    if (Util::ibsee(f.m[id3(i,j,k,f.size)], Cell::Active, Util::Mask1) && (i+j+k) % 2 == 1) {
                        unsigned f12, f13;
                        unsigned f22, f23;
                        unsigned f32, f33;
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
                        real_t   rhs, rc0;
    
                        rhs = psi.m[id3(i  ,j  ,k  ,  psi.size)];
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
                        f12 = Util::ibsee(f12, Cell::Fe, Util::Mask8);
                        f13 = Util::ibsee(f13, Cell::Fe, Util::Mask8);
                        f22 = Util::ibsee(f22, Cell::Fn, Util::Mask8);
                        f23 = Util::ibsee(f23, Cell::Fn, Util::Mask8);
                        f32 = Util::ibsee(f32, Cell::Ft, Util::Mask8);
                        f33 = Util::ibsee(f33, Cell::Ft, Util::Mask8);
                        b12 = Util::ibsee(dom.p.bflag[f12], 0, Util::Mask8);
                        b13 = Util::ibsee(dom.p.bflag[f13], 0, Util::Mask8);
                        b22 = Util::ibsee(dom.p.bflag[f22], 0, Util::Mask8);
                        b23 = Util::ibsee(dom.p.bflag[f23], 0, Util::Mask8);
                        b32 = Util::ibsee(dom.p.bflag[f32], 0, Util::Mask8);
                        b33 = Util::ibsee(dom.p.bflag[f33], 0, Util::Mask8);
                        unsigned feN, feD, fnN, fnD, ftN, ftD;
                        unsigned fwN, fwD, fsN, fsD, fbN, fbD;
                        feN = Util::ibsee(b13, BB::neumann,   Util::Mask1);
                        feD = Util::ibsee(b13, BB::dirichlet, Util::Mask1);
                        fnN = Util::ibsee(b23, BB::neumann,   Util::Mask1);
                        fnD = Util::ibsee(b23, BB::dirichlet, Util::Mask1);
                        ftN = Util::ibsee(b33, BB::neumann,   Util::Mask1);
                        ftD = Util::ibsee(b33, BB::dirichlet, Util::Mask1);
                        fwN = Util::ibsee(b12, BB::neumann,   Util::Mask1);
                        fwD = Util::ibsee(b12, BB::dirichlet, Util::Mask1);
                        fsN = Util::ibsee(b22, BB::neumann,   Util::Mask1);
                        fsD = Util::ibsee(b22, BB::dirichlet, Util::Mask1);
                        fbN = Util::ibsee(b32, BB::neumann,   Util::Mask1);
                        fbD = Util::ibsee(b32, BB::dirichlet, Util::Mask1);
                        
                        real_t jge, jgw, jgn, jgs, jgt, jgb;
                        jge = 0.5 * (g1e + g1c);
                        jgn = 0.5 * (g2n + g2c);
                        jgt = 0.5 * (g3t + g3c);
                        jgw = 0.5 * (g1w + g1c);
                        jgs = 0.5 * (g2s + g2c);
                        jgb = 0.5 * (g3b + g3c);
                        
                        ae1 = (1 - feN) * (1 - feD) * (jge) / det;
                        an1 = (1 - fnN) * (1 - fnD) * (jgn) / det;
                        at1 = (1 - ftN) * (1 - ftD) * (jgt) / det;
                        aw1 = (1 - fwN) * (1 - fwD) * (jgw) / det;
                        as1 = (1 - fsN) * (1 - fsD) * (jgs) / det;
                        ab1 = (1 - fbN) * (1 - fbD) * (jgb) / det;
                        ac0 = - (ae1 + an1 + at1 + aw1 + as1 + ab1) - 2 * (jge * feD + jgw * fwD + jgn * fnD + jgs * fsD + jgt * ftD + jgb * fbD);
                        rc0 = (rhs - ae1 * pe1 - an1 * pn1 - at1 * pt1 - aw1 * pw1 - as1 * ps1 - ab1 * pb1) / ac0 - pc0;
                        
                        p.m[id3(i,j,k,p.size)] = pc0 + omega * rc0;
                        _err += rc0 * rc0;
                        _cnt += 1;
                    }
                }
            }
        }
        err = sqrt(_err / _cnt);
        BB::scalar_outer(p, dom);
    } while (++it < maxit && (ignore || (!ignore && err > torelence)));
}

void MMAC::Poisson::sor(scalar_field<real_t> &p, scalar_field<real_t> &psi, Dom &dom) {
    Ctrl &c = dom.c;
    sor_core(p, psi, dom, c.poisson.omega, c.poisson.maxit, c.poisson.epsilon, false, c.poisson.it, c.poisson.err);
}

void MMAC::Poisson::jacobi(scalar_field<real_t> &p, scalar_field<real_t> &psi, Dom &dom) {
    Ctrl &c = dom.c;
    jacobi_core(p, psi, dom, c.poisson.maxit, c.poisson.epsilon, false, c.poisson.it, c.poisson.err);
}

real_t MMAC::Poisson::dot(scalar_field<real_t> &v1, scalar_field<real_t> &v2, Dom &dom) {
    scalar_field<unsigned> &f = dom.f;
    real_t                sum = 0;
    #pragma acc kernels loop independent collapse(3) reduction(+:sum) present(v1, v2, v3, f, dom) copy(sum)
    for (int i = GUIDE; i < dom.size[0] - GUIDE; i ++) {
        for (int j = GUIDE; j < dom.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < dom.size[2] - GUIDE; k ++) {
                if (Util::ibsee(f.m[id3(i,j,k,f.size)], Cell::Active, Util::Mask1)) {
                    sum += v1.m[id3(i,j,k,v1.size)] * v2.m[id3(i,j,k,v2.size)];
                }
            }
        }
    }
    return sum;
}

void MMAC::Poisson::calc_ax(scalar_field<real_t> &v1, scalar_field<real_t> &v2, Dom &dom) {
    scalar_field<unsigned> &f = dom.f;
    scalar_field<real_t>  &ja = dom.ja;
    vector_field<real_t>   &g = dom.g;
    #pragma acc kernels loop independent collapse(3) present(v1, v2, f, ja, g, dom)
    for (int i = GUIDE; i < dom.size[0] - GUIDE; i ++) {
        for (int j = GUIDE; j < dom.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < dom.size[2] - GUIDE; k ++) {
                if (Util::ibsee(f.m[id3(i,j,k,f.size)], Cell::Active, Util::Mask1)) {
                    unsigned f12, f13;
                    unsigned f22, f23;
                    unsigned f32, f33;
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
                    pc0 =  v1.m[id3(i  ,j  ,k  ,   v1.size)];
                    pe1 =  v1.m[id3(i+1,j  ,k  ,   v1.size)];
                    pn1 =  v1.m[id3(i  ,j+1,k  ,   v1.size)];
                    pt1 =  v1.m[id3(i  ,j  ,k+1,   v1.size)];
                    pw1 =  v1.m[id3(i-1,j  ,k  ,   v1.size)];
                    ps1 =  v1.m[id3(i  ,j-1,k  ,   v1.size)];
                    pb1 =  v1.m[id3(i  ,j  ,k-1,   v1.size)];
                    f13 =   f.m[id3(i  ,j  ,k  ,    f.size)];
                    f23 =   f.m[id3(i  ,j  ,k  ,    f.size)];
                    f33 =   f.m[id3(i  ,j  ,k  ,    f.size)];
                    f12 =   f.m[id3(i-1,j  ,k  ,    f.size)];
                    f22 =   f.m[id3(i  ,j-1,k  ,    f.size)];
                    f32 =   f.m[id3(i  ,j  ,k-1,    f.size)];
                    f12 = Util::ibsee(f12, Cell::Fe, Util::Mask8);
                    f13 = Util::ibsee(f13, Cell::Fe, Util::Mask8);
                    f22 = Util::ibsee(f22, Cell::Fn, Util::Mask8);
                    f23 = Util::ibsee(f23, Cell::Fn, Util::Mask8);
                    f32 = Util::ibsee(f32, Cell::Ft, Util::Mask8);
                    f33 = Util::ibsee(f33, Cell::Ft, Util::Mask8);
                    b12 = Util::ibsee(dom.p.bflag[f12], 0, Util::Mask8);
                    b13 = Util::ibsee(dom.p.bflag[f13], 0, Util::Mask8);
                    b22 = Util::ibsee(dom.p.bflag[f22], 0, Util::Mask8);
                    b23 = Util::ibsee(dom.p.bflag[f23], 0, Util::Mask8);
                    b32 = Util::ibsee(dom.p.bflag[f32], 0, Util::Mask8);
                    b33 = Util::ibsee(dom.p.bflag[f33], 0, Util::Mask8);
                    unsigned feN, feD, fnN, fnD, ftN, ftD;
                    unsigned fwN, fwD, fsN, fsD, fbN, fbD;
                    feN = Util::ibsee(b13, BB::neumann,   Util::Mask1);
                    feD = Util::ibsee(b13, BB::dirichlet, Util::Mask1);
                    fnN = Util::ibsee(b23, BB::neumann,   Util::Mask1);
                    fnD = Util::ibsee(b23, BB::dirichlet, Util::Mask1);
                    ftN = Util::ibsee(b33, BB::neumann,   Util::Mask1);
                    ftD = Util::ibsee(b33, BB::dirichlet, Util::Mask1);
                    fwN = Util::ibsee(b12, BB::neumann,   Util::Mask1);
                    fwD = Util::ibsee(b12, BB::dirichlet, Util::Mask1);
                    fsN = Util::ibsee(b22, BB::neumann,   Util::Mask1);
                    fsD = Util::ibsee(b22, BB::dirichlet, Util::Mask1);
                    fbN = Util::ibsee(b32, BB::neumann,   Util::Mask1);
                    fbD = Util::ibsee(b32, BB::dirichlet, Util::Mask1);
                        
                    real_t jge, jgw, jgn, jgs, jgt, jgb;
                    jge = 0.5 * (g1e + g1c);
                    jgn = 0.5 * (g2n + g2c);
                    jgt = 0.5 * (g3t + g3c);
                    jgw = 0.5 * (g1w + g1c);
                    jgs = 0.5 * (g2s + g2c);
                    jgb = 0.5 * (g3b + g3c);
                        
                    ae1 = (1 - feN) * (1 - feD) * (jge) / det;
                    an1 = (1 - fnN) * (1 - fnD) * (jgn) / det;
                    at1 = (1 - ftN) * (1 - ftD) * (jgt) / det;
                    aw1 = (1 - fwN) * (1 - fwD) * (jgw) / det;
                    as1 = (1 - fsN) * (1 - fsD) * (jgs) / det;
                    ab1 = (1 - fbN) * (1 - fbD) * (jgb) / det;
                    ac0 = - (ae1 + an1 + at1 + aw1 + as1 + ab1) - 2 * (jge * feD + jgw * fwD + jgn * fnD + jgs * fsD + jgt * ftD + jgb * fbD);

                    v2.m[id3(i,j,k,v2.size)] = ac0 * pc0 + ae1 * pe1 + an1 * pn1 + at1 * pt1 + aw1 * pw1 + as1 * ps1 + ab1 * pb1;
                }
            }
        }
    }
}

void MMAC::Poisson::triad(scalar_field<real_t> &v1, scalar_field<real_t> &v2, scalar_field<real_t> &v3, real_t a, Dom &dom) {
    #pragma acc kernels loop independent collapse(3) present(v1, v2, v3, dom) copyin(a)
    for (int i = 0; i < dom.size[0]; i ++) {
        for (int j = 0; j < dom.size[1]; j ++) {
            for (int k = 0; k < dom.size[2]; k ++) {
                v3.m[id3(i,j,k,v3.size)] = v1.m[id3(i,j,k,v1.size)] + a * v2.m[id3(i,j,k,v2.size)];
            }
        }
    }
}

void MMAC::Poisson::calc_res(scalar_field<real_t> &p, scalar_field<real_t> &b, scalar_field<real_t> &r, Dom &dom) {
    scalar_field<unsigned> &f = dom.f;
    scalar_field<real_t>  &ja = dom.ja;
    vector_field<real_t>   &g = dom.g;
    #pragma acc kernels loop independent collapse(3) present(p, b, r, f, ja, g, dom)
    for (int i = GUIDE; i < dom.size[0] - GUIDE; i ++) {
        for (int j = GUIDE; j < dom.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < dom.size[2] - GUIDE; k ++) {
                if (Util::ibsee(f.m[id3(i,j,k,f.size)], Cell::Active, Util::Mask1)) {
                    unsigned f12, f13;
                    unsigned f22, f23;
                    unsigned f32, f33;
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
                    f12 = Util::ibsee(f12, Cell::Fe, Util::Mask8);
                    f13 = Util::ibsee(f13, Cell::Fe, Util::Mask8);
                    f22 = Util::ibsee(f22, Cell::Fn, Util::Mask8);
                    f23 = Util::ibsee(f23, Cell::Fn, Util::Mask8);
                    f32 = Util::ibsee(f32, Cell::Ft, Util::Mask8);
                    f33 = Util::ibsee(f33, Cell::Ft, Util::Mask8);
                    b12 = Util::ibsee(dom.p.bflag[f12], 0, Util::Mask8);
                    b13 = Util::ibsee(dom.p.bflag[f13], 0, Util::Mask8);
                    b22 = Util::ibsee(dom.p.bflag[f22], 0, Util::Mask8);
                    b23 = Util::ibsee(dom.p.bflag[f23], 0, Util::Mask8);
                    b32 = Util::ibsee(dom.p.bflag[f32], 0, Util::Mask8);
                    b33 = Util::ibsee(dom.p.bflag[f33], 0, Util::Mask8);
                    unsigned feN, feD, fnN, fnD, ftN, ftD;
                    unsigned fwN, fwD, fsN, fsD, fbN, fbD;
                    feN = Util::ibsee(b13, BB::neumann,   Util::Mask1);
                    feD = Util::ibsee(b13, BB::dirichlet, Util::Mask1);
                    fnN = Util::ibsee(b23, BB::neumann,   Util::Mask1);
                    fnD = Util::ibsee(b23, BB::dirichlet, Util::Mask1);
                    ftN = Util::ibsee(b33, BB::neumann,   Util::Mask1);
                    ftD = Util::ibsee(b33, BB::dirichlet, Util::Mask1);
                    fwN = Util::ibsee(b12, BB::neumann,   Util::Mask1);
                    fwD = Util::ibsee(b12, BB::dirichlet, Util::Mask1);
                    fsN = Util::ibsee(b22, BB::neumann,   Util::Mask1);
                    fsD = Util::ibsee(b22, BB::dirichlet, Util::Mask1);
                    fbN = Util::ibsee(b32, BB::neumann,   Util::Mask1);
                    fbD = Util::ibsee(b32, BB::dirichlet, Util::Mask1);
                        
                    real_t jge, jgw, jgn, jgs, jgt, jgb;
                    jge = 0.5 * (g1e + g1c);
                    jgn = 0.5 * (g2n + g2c);
                    jgt = 0.5 * (g3t + g3c);
                    jgw = 0.5 * (g1w + g1c);
                    jgs = 0.5 * (g2s + g2c);
                    jgb = 0.5 * (g3b + g3c);
                        
                    ae1 = (1 - feN) * (1 - feD) * (jge) / det;
                    an1 = (1 - fnN) * (1 - fnD) * (jgn) / det;
                    at1 = (1 - ftN) * (1 - ftD) * (jgt) / det;
                    aw1 = (1 - fwN) * (1 - fwD) * (jgw) / det;
                    as1 = (1 - fsN) * (1 - fsD) * (jgs) / det;
                    ab1 = (1 - fbN) * (1 - fbD) * (jgb) / det;
                    ac0 = - (ae1 + an1 + at1 + aw1 + as1 + ab1) - 2 * (jge * feD + jgw * fwD + jgn * fnD + jgs * fsD + jgt * ftD + jgb * fbD);

                    r.m[id3(i,j,k,r.size)] = b.m[id3(i,j,k,b.size)] - (ac0 * pc0 + ae1 * pe1 + an1 * pn1 + at1 * pt1 + aw1 * pw1 + as1 * ps1 + ab1 * pb1);
                }
            }
        }
    }
}

void MMAC::Poisson::pbicgstab(scalar_field<real_t> &p, scalar_field<real_t> &psi, Dom &dom) {
    Ctrl                   &c = dom.c;
    scalar_field<unsigned> &f = dom.f;

    real_t dummy_real = 0.0;
    int    dummy_int  = 0;
    real_t err        = 0;
    int    cnt        = 0;

    pcg_rho     = 0.0;
    pcg_rho_old = 1.0;
    pcg_alpha   = 1.0;
    pcg_beta    = 0.0;
    pcg_omega   = 1.0;
    pcg_it      = 0;

    calc_res(p, psi, pcg_r, dom);
    pcg_r0.copy(pcg_r);
    pcg_q.clear();

    c.poisson.it = 0;
    do {
        err = 0;
        cnt = 0;
        pcg_rho = dot(pcg_r, pcg_r0, dom);
        if (fabs(pcg_rho) < FLT_MIN) {
            c.poisson.err = fabs(pcg_rho);
            break;
        }

        if (c.poisson.it == 0) {
            pcg_p.copy(pcg_r);
        } else {
            pcg_beta = (pcg_rho / pcg_rho_old) * (pcg_alpha / pcg_omega);

            #pragma acc update device(pcg_beta, pcg_omega)
            #pragma acc kernels loop independent collapse(3) present(this[0:1], dom, f, pcg_p, pcg_r, pcg_q)
            for (int i = GUIDE; i < dom.size[0] - GUIDE; i ++) {
                for (int j = GUIDE; j < dom.size[1] - GUIDE; j ++) {
                    for (int k = GUIDE; k < dom.size[2] - GUIDE; k ++) {
                        if (Util::ibsee(f.m[id3(i,j,k,f.size)], Cell::Active, Util::Mask1)) {
                            pcg_p.m[id3(i,j,k,pcg_p.size)] = pcg_r.m[id3(i,j,k,pcg_r.size)] + pcg_beta * (pcg_p.m[id3(i,j,k,pcg_p.size)] - pcg_omega * pcg_q.m[id3(i,j,k,pcg_q.size)]);
                        }
                    }
                }
            }
        }

        pcg_p_.clear();
        sor_core(pcg_p_, pcg_p, dom, 1.2, 5, 0, true, dummy_int, dummy_real);
        calc_ax(pcg_p_, pcg_q, dom);
        pcg_alpha = pcg_rho / (dot(pcg_r0, pcg_q, dom));
        triad(pcg_r, pcg_q, pcg_s, -pcg_alpha, dom);
        pcg_s_.clear();
        sor_core(pcg_s_, pcg_s, dom, 1.2, 5, 0, true, dummy_int, dummy_real);
        calc_ax(pcg_s_, pcg_t, dom);
        pcg_omega = dot(pcg_t, pcg_s, dom) / dot(pcg_t, pcg_t, dom);

        #pragma acc update device(pcg_alpha, pcg_omega)
        #pragma acc kernels loop independent collapse(3) reduction(+:err, cnt) present(this[0:1], dom, f, p, pcg_p_, pcg_s_) copy(err, cnt)
        for (int i = GUIDE; i < dom.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < dom.size[1] - GUIDE; j ++) {
                for (int k = GUIDE; k < dom.size[2] - GUIDE; k ++) {
                    if (Util::ibsee(f.m[id3(i,j,k,f.size)], Cell::Active, Util::Mask1)) {
                        real_t r = pcg_alpha * pcg_p_.m[id3(i,j,k,pcg_p_.size)] + pcg_omega * pcg_s_.m[id3(i,j,k,pcg_s_.size)];
                        p.m[id3(i,j,k,p.size)] += r;
                        err += r * r;
                        cnt += 1;
                    }
                }
            }
        }
        BB::scalar_outer(p, dom);
        triad(pcg_s, pcg_t, pcg_r, -pcg_omega, dom);

        // #pragma acc kernels loop independent collapse(3) reduction(+:err, cnt) present(this[0:1], f, dom, pcg_r) copy(err, cnt)
        // for (int i = GUIDE; i < dom.size[0] - GUIDE; i ++) {
        //     for (int j = GUIDE; j < dom.size[1] - GUIDE; j ++) {
        //         for (int k = GUIDE; k < dom.size[2] - GUIDE; k ++) {
        //             if (Util::ibsee(f.m[id3(i,j,k,f.size)], Cell::Active, Util::Mask1)) {
        //                 real_t r = pcg_r.m[id3(i,j,k,pcg_r.size)];
        //                 err += r * r;
        //                 cnt += 1;
        //             }
        //         }
        //     }
        // }

        c.poisson.err = sqrt(err / cnt);
        pcg_rho_old = pcg_rho;
        pcg_it ++;
    } while (++c.poisson.it < c.poisson.maxit && c.poisson.err > c.poisson.epsilon);
}

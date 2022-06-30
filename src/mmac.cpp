#include <stdio.h>
#include "../include/scheme.h"
#include "../include/mmac.h"
#include "../include/flag.h"
#include "../include/boundary.h"

void MMAC::divergence_velocity(vector_field<real_t> &uu, scalar_field<real_t> &div, real_t &div_norm, Dom &dom) {
    scalar_field<real_t>   &ja = dom.ja;
    scalar_field<unsigned> &f  = dom.f;
    real_t sum = 0;
    real_t cnt = 0;

    #pragma acc kernels loop independent collapse(3) reduction(+:sum, cnt) present(uu, div, ja, f, dom) copy(sum, cnt)
    for (int i = GUIDE; i < dom.size[0] - GUIDE; i ++) {
        for (int j = GUIDE; j < dom.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < dom.size[2] - GUIDE; k ++) {
                if (Util::ibsee(f.m[id3(i,j,k,f.size)], Flag::Active, Util::Mask1)) {
                    real_t ufe, vfn, wft;
                    real_t ufw, vfs, wfb;
                    real_t det;
                    real_t dc0;

                    ufe = uu.m[id4(i  ,j  ,k  ,0,uu.size)];
                    vfn = uu.m[id4(i  ,j  ,k  ,1,uu.size)];
                    wft = uu.m[id4(i  ,j  ,k  ,2,uu.size)];
                    ufw = uu.m[id4(i-1,j  ,k  ,0,uu.size)];
                    vfs = uu.m[id4(i  ,j-1,k  ,1,uu.size)];
                    wfb = uu.m[id4(i  ,j  ,k-1,2,uu.size)];
                    det = ja.m[id3(i  ,j  ,k  ,  ja.size)];

                    dc0 = (ufe + vfn + wft - ufw - vfs - wfb) / det;
                    sum += dc0 * dc0 * det;
                    cnt += det;

                    div.m[id3(i,j,k,div.size)] = dc0;
                }
            }
        }
    }
    div_norm = sqrt(sum / cnt);
}

void MMAC::interpolate_velocity(vector_field<real_t> &u, vector_field<real_t> &uc, vector_field<real_t> &uu, Dom &dom) {
    scalar_field<unsigned> &f = dom.f;
    vector_field<real_t>   &x  = dom.x;
    vector_field<real_t>   &kx = dom.kx;
    scalar_field<real_t>   &ja = dom.ja;

    #pragma acc kernels loop independent collapse(3) present(u, uc, uu, f, x, kx, ja, dom)
    for (int i = GUIDE - 1; i <= dom.size[0] - GUIDE; i ++) {
        for (int j = GUIDE - 1; j <= dom.size[1] - GUIDE; j ++) {
            for (int k = GUIDE - 1; k <= dom.size[2] - GUIDE; k ++) {
                real_t k1, k2, k3, u1, u2, u3, de;

                u1 =  u.m[id4(i,j,k,0, u.size)];
                u2 =  u.m[id4(i,j,k,1, u.size)];
                u3 =  u.m[id4(i,j,k,2, u.size)];
                k1 = kx.m[id4(i,j,k,0,kx.size)];
                k2 = kx.m[id4(i,j,k,1,kx.size)];
                k3 = kx.m[id4(i,j,k,2,kx.size)];
                de = ja.m[id3(i,j,k,  ja.size)];

                uc.m[id4(i,j,k,0,uc.size)] = de * k1 * u1;
                uc.m[id4(i,j,k,1,uc.size)] = de * k2 * u2;
                uc.m[id4(i,j,k,2,uc.size)] = de * k3 * u3;
            }
        }
    }

    #pragma acc kernels loop independent collapse(3) present(u, uc, uu, f, x, kx, ja, dom)
    for (int i = GUIDE - 1; i < dom.size[0] - GUIDE; i ++) {
        for (int j = GUIDE - 1; j < dom.size[1] - GUIDE; j ++) {
            for (int k = GUIDE - 1; k < dom.size[2] - GUIDE; k ++) {
                unsigned int ff;
                unsigned int f3, m3, b3;
                real_t u0, u1, uf;
                real_t kk;
                real_t ref, dis;
                ff = f.m[id3(i,j,k,f.size)];

                if (j >= GUIDE && k >= GUIDE) {
                    m3 = Util::ibsee(ff, Flag::Me, Util::Mask1);
                    f3 = Util::ibsee(ff, Flag::Fe, Util::Mask8);
                    if (Util::ibsee(u.bflag[f3], BB::uu_locked, Util::Mask1)) {
                        uu.m[id4(i,j,k,0,uu.size)] = dom.uu.m[id4(i,j,k,0,dom.uu.size)];
                    } else {
                        b3 = Util::ibsee(u.bflag[f3], 0, Util::Mask8);
                        if (b3) {
                            u0 = u.m[id4(i  ,j,k,0,u.size)];
                            u1 = u.m[id4(i+1,j,k,0,u.size)];
                            BB::pre(m3, u0, u1, ref, dis);
                            uf = BB::eva(b3, ref, dis, u.b[id2(f3,0,u.bsize)]);
                            kk = 1 /    (x.m[id4(i+1,j,k,0,x.size)] -  x.m[id4(i,j,k,0,x.size)]);
                            uf = 0.5 * (ja.m[id3(i+1,j,k, ja.size)] + ja.m[id3(i,j,k, ja.size)]) * kk * uf;
                        } else {
                            u0 = uc.m[id4(i  ,j,k,0,uc.size)];
                            u1 = uc.m[id4(i+1,j,k,0,uc.size)];
                            uf = 0.5 * (u0 + u1);
                        }
                        uu.m[id4(i,j,k,0,uu.size)] = uf;
                    }
                }

                if (i >= GUIDE && k >= GUIDE) {
                    m3 = Util::ibsee(ff, Flag::Mn, Util::Mask1);
                    f3 = Util::ibsee(ff, Flag::Fn, Util::Mask8);
                    if (Util::ibsee(u.bflag[f3], BB::uu_locked, Util::Mask1)) {
                        uu.m[id4(i,j,k,1,uu.size)] = dom.uu.m[id4(i,j,k,1,dom.uu.size)];
                    } else {
                        b3 = Util::ibsee(u.bflag[f3], 0, Util::Mask8);
                        if (b3) {
                            u0 = u.m[id4(i,j  ,k,1,u.size)];
                            u1 = u.m[id4(i,j+1,k,1,u.size)];
                            BB::pre(m3, u0, u1, ref, dis);
                            uf = BB::eva(b3, ref, dis, u.b[id2(f3,1,u.bsize)]);
                            kk = 1 /    (x.m[id4(i,j+1,k,1,x.size)] -  x.m[id4(i,j,k,1,x.size)]);
                            uf = 0.5 * (ja.m[id3(i,j+1,k, ja.size)] + ja.m[id3(i,j,k, ja.size)]) * kk * uf;
                        } else {
                            u0 = uc.m[id4(i,j  ,k,1,uc.size)];
                            u1 = uc.m[id4(i,j+1,k,1,uc.size)];
                            uf = 0.5 * (u0 + u1);
                        }
                        uu.m[id4(i,j,k,1,uu.size)] = uf;
                    }
                }

                if (i >= GUIDE && j >= GUIDE) {
                    m3 = Util::ibsee(ff, Flag::Mt, Util::Mask1);
                    f3 = Util::ibsee(ff, Flag::Ft, Util::Mask8);
                    if (Util::ibsee(u.bflag[f3], BB::uu_locked, Util::Mask1)) {
                        uu.m[id4(i,j,k,2,uu.size)] = dom.uu.m[id4(i,j,k,2,dom.uu.size)];
                    } else {
                        b3 = Util::ibsee(u.bflag[f3], 0, Util::Mask8);
                        if (b3) {
                            u0 = u.m[id4(i,j,k  ,2,u.size)];
                            u1 = u.m[id4(i,j,k+1,2,u.size)];
                            BB::pre(m3, u0, u1, ref, dis);
                            uf = BB::eva(b3, ref, dis, u.b[id2(f3,2,u.bsize)]);
                            kk = 1 /    (x.m[id4(i,j,k+1,2,x.size)] -  x.m[id4(i,j,k,2,x.size)]);
                            uf = 0.5 * (ja.m[id3(i,j,k+1, ja.size)] + ja.m[id3(i,j,k, ja.size)]) * kk * uf;
                        } else {
                            u0 = uc.m[id4(i,j,k  ,2,uc.size)];
                            u1 = uc.m[id4(i,j,k+1,2,uc.size)];
                            uf = 0.5 * (u0 + u1);
                        }
                        uu.m[id4(i,j,k,2,uu.size)] = uf;
                    }
                }
            }
        }
    }
}

void MMAC::pseudo_velocity(vector_field<real_t> &u, vector_field<real_t> &uu, vector_field<real_t> &ua, scalar_field<real_t> &nut, Dom &dom) {
    Ctrl                   &c  = dom.c;
    scalar_field<unsigned> &f  = dom.f;
    vector_field<real_t>   &g  = dom.g;
    scalar_field<real_t>   &ja = dom.ja;

    real_t kappa = 1.0 / 3.0;
    real_t beta  = (3.0 - kappa) / (1.0 - kappa);
    real_t m1    = 1.0 - kappa;
    real_t m2    = 1.0 + kappa;

    #pragma acc kernels loop independent collapse(3) present(u, uu, ua, nut, f, g, ja, c, dom) copyin(kappa, beta, m1, m2)
    for (int i = GUIDE; i < dom.size[0] - GUIDE; i ++) {
        for (int j = GUIDE; j < dom.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < dom.size[2] - GUIDE; k++) {
                if (Util::ibsee(f.m[id3(i,j,k,f.size)], Flag::Active, Util::Mask1)) {
                    real_t uc0;
                    real_t ue1, un1, ut1;
                    real_t ue2, un2, ut2;
                    real_t uw1, us1, ub1;
                    real_t uw2, us2, ub2;
                    real_t ufe, vfn, wft;
                    real_t ufw, vfs, wfb;
                    real_t de1, dn1, dt1;
                    real_t dw1, ds1, db1;
                    real_t nc0;
                    real_t ne1, nn1, nt1;
                    real_t nw1, ns1, nb1;
                    real_t det;
                    real_t ref, dis;
                    real_t ad1, ad2, ad3;
                    real_t vi1, vi2, vi3;
                    real_t g1c, g2c, g3c;
                    real_t g1e, g2n, g3t;
                    real_t g1w, g2s, g3b;
                    unsigned f11, f12, f13, f14;
                    unsigned f21, f22, f23, f24;
                    unsigned f31, f32, f33, f34;
                    unsigned b11, b12, b13, b14;
                    unsigned b21, b22, b23, b24;
                    unsigned b31, b32, b33, b34;
                    unsigned m11, m12, m13, m14;
                    unsigned m21, m22, m23, m24;
                    unsigned m31, m32, m33, m34;
                    int epe, epn, ept;
                    int epw, eps, epb;
                    real_t u1, u2, u3, mag;

                    u1  = u.m[id4(i,j,k,0,u.size)];
                    u2  = u.m[id4(i,j,k,1,u.size)];
                    u3  = u.m[id4(i,j,k,2,u.size)];
                    mag = sqrt(u1 * u1 + u2 * u2 + u3 * u3);

                    ufe = uu.m[id4(i  ,j  ,k  ,0,uu.size)];
                    vfn = uu.m[id4(i  ,j  ,k  ,1,uu.size)];
                    wft = uu.m[id4(i  ,j  ,k  ,2,uu.size)];
                    ufw = uu.m[id4(i-1,j  ,k  ,0,uu.size)];
                    vfs = uu.m[id4(i  ,j-1,k  ,1,uu.size)];
                    wfb = uu.m[id4(i  ,j  ,k-1,2,uu.size)];
                    det = ja.m[id3(i  ,j  ,k  ,  ja.size)];
                    g1c =  g.m[id4(i  ,j  ,k  ,0, g.size)];
                    g2c =  g.m[id4(i  ,j  ,k  ,1, g.size)];
                    g3c =  g.m[id4(i  ,j  ,k  ,2, g.size)];
                    g1e =  g.m[id4(i+1,j  ,k  ,0, g.size)];
                    g2n =  g.m[id4(i  ,j+1,k  ,1, g.size)];
                    g3t =  g.m[id4(i  ,j  ,k+1,2, g.size)];
                    g1w =  g.m[id4(i-1,j  ,k  ,0, g.size)];
                    g2s =  g.m[id4(i  ,j-1,k  ,1, g.size)];
                    g3b =  g.m[id4(i  ,j  ,k-1,2, g.size)];
                    f11 =  f.m[id3(i-2,j  ,k  ,   f.size)];
                    f12 =  f.m[id3(i-1,j  ,k  ,   f.size)];
                    f13 =  f.m[id3(i  ,j  ,k  ,   f.size)];
                    f14 =  f.m[id3(i+1,j  ,k  ,   f.size)];
                    f21 =  f.m[id3(i  ,j-2,k  ,   f.size)];
                    f22 =  f.m[id3(i  ,j-1,k  ,   f.size)];
                    f23 =  f.m[id3(i  ,j  ,k  ,   f.size)];
                    f24 =  f.m[id3(i  ,j+1,k  ,   f.size)];
                    f31 =  f.m[id3(i  ,j  ,k-2,   f.size)];
                    f32 =  f.m[id3(i  ,j  ,k-1,   f.size)];
                    f33 =  f.m[id3(i  ,j  ,k  ,   f.size)];
                    f34 =  f.m[id3(i  ,j  ,k+1,   f.size)];
                    m11 = Util::ibsee(f11, Flag::Me, Util::Mask1);
                    m12 = Util::ibsee(f12, Flag::Me, Util::Mask1);
                    m13 = Util::ibsee(f13, Flag::Me, Util::Mask1);
                    m14 = Util::ibsee(f14, Flag::Me, Util::Mask1);
                    m21 = Util::ibsee(f21, Flag::Mn, Util::Mask1);
                    m22 = Util::ibsee(f22, Flag::Mn, Util::Mask1);
                    m23 = Util::ibsee(f23, Flag::Mn, Util::Mask1);
                    m24 = Util::ibsee(f24, Flag::Mn, Util::Mask1);
                    m31 = Util::ibsee(f31, Flag::Mt, Util::Mask1);
                    m32 = Util::ibsee(f32, Flag::Mt, Util::Mask1);
                    m33 = Util::ibsee(f33, Flag::Mt, Util::Mask1);
                    m34 = Util::ibsee(f34, Flag::Mt, Util::Mask1);
                    f11 = Util::ibsee(f11, Flag::Fe, Util::Mask8);
                    f12 = Util::ibsee(f12, Flag::Fe, Util::Mask8);
                    f13 = Util::ibsee(f13, Flag::Fe, Util::Mask8);
                    f14 = Util::ibsee(f14, Flag::Fe, Util::Mask8);
                    f21 = Util::ibsee(f21, Flag::Fn, Util::Mask8);
                    f22 = Util::ibsee(f22, Flag::Fn, Util::Mask8);
                    f23 = Util::ibsee(f23, Flag::Fn, Util::Mask8);
                    f24 = Util::ibsee(f24, Flag::Fn, Util::Mask8);
                    f31 = Util::ibsee(f31, Flag::Ft, Util::Mask8);
                    f32 = Util::ibsee(f32, Flag::Ft, Util::Mask8);
                    f33 = Util::ibsee(f33, Flag::Ft, Util::Mask8);
                    f34 = Util::ibsee(f34, Flag::Ft, Util::Mask8);

                    nc0 = nut.m[id3(i  ,j  ,k  ,nut.size)];
                    ne1 = nut.m[id3(i+1,j  ,k  ,nut.size)];
                    nn1 = nut.m[id3(i  ,j+1,k  ,nut.size)];
                    nt1 = nut.m[id3(i  ,j  ,k+1,nut.size)];
                    nw1 = nut.m[id3(i-1,j  ,k  ,nut.size)];
                    ns1 = nut.m[id3(i  ,j-1,k  ,nut.size)];
                    nb1 = nut.m[id3(i  ,j  ,k-1,nut.size)];
                    b12 = Util::ibsee(nut.bflag[f12], 0, Util::Mask8);
                    b13 = Util::ibsee(nut.bflag[f13], 0, Util::Mask8);
                    b22 = Util::ibsee(nut.bflag[f22], 0, Util::Mask8);
                    b23 = Util::ibsee(nut.bflag[f23], 0, Util::Mask8);
                    b32 = Util::ibsee(nut.bflag[f32], 0, Util::Mask8);
                    b33 = Util::ibsee(nut.bflag[f33], 0, Util::Mask8);
                    if (b13) {
                        BB::pre(m13, nc0, ne1, ref, dis);
                        ne1 = 2 * BB::eva(b13, ref, dis, nut.b[f13]) - nc0;
                    }
                    if (b23) {
                        BB::pre(m23, nc0, nn1, ref, dis);
                        nn1 = 2 * BB::eva(b23, ref, dis, nut.b[f23]) - nc0;
                    }
                    if (b33) {
                        BB::pre(m33, nc0, nt1, ref, dis);
                        nt1 = 2 * BB::eva(b33, ref, dis, nut.b[f33]) - nc0;
                    }
                    if (b12) {
                        BB::pre(m12, nw1, nc0, ref, dis);
                        nw1 = 2 * BB::eva(b12, ref, dis, nut.b[f12]) - nc0;
                    }
                    if (b22) {
                        BB::pre(m22, ns1, nc0, ref, dis);
                        ns1 = 2 * BB::eva(b22, ref, dis, nut.b[f22]) - nc0;
                    }
                    if (b32) {
                        BB::pre(m32, nb1, nc0, ref, dis);
                        nb1 = 2 * BB::eva(b32, ref, dis, nut.b[f32]) - nc0;
                    }

                    b11 = Util::ibsee(u.bflag[f11], 0, Util::Mask8);
                    b12 = Util::ibsee(u.bflag[f12], 0, Util::Mask8);
                    b13 = Util::ibsee(u.bflag[f13], 0, Util::Mask8);
                    b14 = Util::ibsee(u.bflag[f14], 0, Util::Mask8);
                    b21 = Util::ibsee(u.bflag[f21], 0, Util::Mask8);
                    b22 = Util::ibsee(u.bflag[f22], 0, Util::Mask8);
                    b23 = Util::ibsee(u.bflag[f23], 0, Util::Mask8);
                    b24 = Util::ibsee(u.bflag[f24], 0, Util::Mask8);
                    b31 = Util::ibsee(u.bflag[f31], 0, Util::Mask8);
                    b32 = Util::ibsee(u.bflag[f32], 0, Util::Mask8);
                    b33 = Util::ibsee(u.bflag[f33], 0, Util::Mask8);
                    b34 = Util::ibsee(u.bflag[f34], 0, Util::Mask8);
                    epe = (b13)? 0 : 1;
                    epn = (b23)? 0 : 1;
                    ept = (b33)? 0 : 1;
                    epw = (b12)? 0 : 1;
                    eps = (b22)? 0 : 1;
                    epb = (b32)? 0 : 1;
                    de1 = dn1 = dt1 = dw1 = ds1 = db1 = 0;

                    uc0 = u1;
                    uw2 = u.m[id4(i-2,j  ,k  ,0,u.size)];
                    uw1 = u.m[id4(i-1,j  ,k  ,0,u.size)];
                    ue1 = u.m[id4(i+1,j  ,k  ,0,u.size)];
                    ue2 = u.m[id4(i+2,j  ,k  ,0,u.size)];
                    us2 = u.m[id4(i  ,j-2,k  ,0,u.size)];
                    us1 = u.m[id4(i  ,j-1,k  ,0,u.size)];
                    un1 = u.m[id4(i  ,j+1,k  ,0,u.size)];
                    un2 = u.m[id4(i  ,j+2,k  ,0,u.size)];
                    ub2 = u.m[id4(i  ,j  ,k-2,0,u.size)];
                    ub1 = u.m[id4(i  ,j  ,k-1,0,u.size)];
                    ut1 = u.m[id4(i  ,j  ,k+1,0,u.size)];
                    ut2 = u.m[id4(i  ,j  ,k+2,0,u.size)];
                    if (b13) {
                        BB::pre(m13, uc0, ue1, ref, dis);
                        ue1 = 2 * BB::eva(b13, ref, dis, u.b[id2(f13,0,u.bsize)]) - uc0;
                    }
                    if (b23) {
                        BB::pre(m23, uc0, un1, ref, dis);
                        un1 = 2 * BB::eva(b23, ref, dis, u.b[id2(f23,0,u.bsize)]) - uc0;
                    }
                    if (b33) {
                        BB::pre(m33, uc0, ut1, ref, dis);
                        ut1 = 2 * BB::eva(b33, ref, dis, u.b[id2(f33,0,u.bsize)]) - uc0;
                    }
                    if (b12) {
                        BB::pre(m12, uw1, uc0, ref, dis);
                        uw1 = 2 * BB::eva(b12, ref, dis, u.b[id2(f12,0,u.bsize)]) - uc0;
                    }
                    if (b22) {
                        BB::pre(m22, us1, uc0, ref, dis);
                        us1 = 2 * BB::eva(b22, ref, dis, u.b[id2(f22,0,u.bsize)]) - uc0;
                    }
                    if (b32) {
                        BB::pre(m32, ub1, uc0, ref, dis);
                        ub1 = 2 * BB::eva(b32, ref, dis, u.b[id2(f32,0,u.bsize)]) - uc0;
                    }
                    if (b14) {
                        BB::pre(m14, ue1, ue2, ref, dis);
                        ue2 = 2 * BB::eva(b14, ref, dis, u.b[id2(f14,0,u.bsize)]) - ue1;
                    }
                    if (b24) {
                        BB::pre(m24, un1, un2, ref, dis);
                        un2 = 2 * BB::eva(b24, ref, dis, u.b[id2(f24,0,u.bsize)]) - un1;
                    }
                    if (b34) {
                        BB::pre(m34, ut1, ut2, ref, dis);
                        ut2 = 2 * BB::eva(b34, ref, dis, u.b[id2(f34,0,u.bsize)]) - ut1;
                    }
                    if (b11) {
                        BB::pre(m11, uw2, uw1, ref, dis);
                        uw2 = 2 * BB::eva(b11, ref, dis, u.b[id2(f11,0,u.bsize)]) - uw1;
                    }
                    if (b21) {
                        BB::pre(m21, us2, us1, ref, dis);
                        us2 = 2 * BB::eva(b21, ref, dis, u.b[id2(f21,0,u.bsize)]) - us1;
                    }
                    if (b31) {
                        BB::pre(m31, ub2, ub1, ref, dis);
                        ub2 = 2 * BB::eva(b31, ref, dis, u.b[id2(f31,0,u.bsize)]) - ub1;
                    }
                    ad1 = Scheme::ffvc_muscl(uc0, ue1, ue2, un1, un2, ut1, ut2, uw1, uw2, us1, us2, ub1, ub2, ufe, vfn, wft, ufw, vfs, wfb, det, epe, epn, ept, epw, eps, epb, kappa, beta, m1, m2);
                    vi1 = Scheme::viscosity(0, 0, 0, 0, 0, 0, mag, uc0, ue1, un1, ut1, uw1, us1, ub1, de1, dn1, dt1, dw1, ds1, db1, nc0, ne1, nn1, nt1, nw1, ns1, nb1, det, g1c, g2c, g3c, g1e, g2n, g3t, g1w, g2s, g3b, c.flow.ri);

                    uc0 = u2;
                    uw2 = u.m[id4(i-2,j  ,k  ,1,u.size)];
                    uw1 = u.m[id4(i-1,j  ,k  ,1,u.size)];
                    ue1 = u.m[id4(i+1,j  ,k  ,1,u.size)];
                    ue2 = u.m[id4(i+2,j  ,k  ,1,u.size)];
                    us2 = u.m[id4(i  ,j-2,k  ,1,u.size)];
                    us1 = u.m[id4(i  ,j-1,k  ,1,u.size)];
                    un1 = u.m[id4(i  ,j+1,k  ,1,u.size)];
                    un2 = u.m[id4(i  ,j+2,k  ,1,u.size)];
                    ub2 = u.m[id4(i  ,j  ,k-2,1,u.size)];
                    ub1 = u.m[id4(i  ,j  ,k-1,1,u.size)];
                    ut1 = u.m[id4(i  ,j  ,k+1,1,u.size)];
                    ut2 = u.m[id4(i  ,j  ,k+2,1,u.size)];
                    if (b13) {
                        BB::pre(m13, uc0, ue1, ref, dis);
                        ue1 = 2 * BB::eva(b13, ref, dis, u.b[id2(f13,1,u.bsize)]) - uc0;
                    }
                    if (b23) {
                        BB::pre(m23, uc0, un1, ref, dis);
                        un1 = 2 * BB::eva(b23, ref, dis, u.b[id2(f23,1,u.bsize)]) - uc0;
                    }
                    if (b33) {
                        BB::pre(m33, uc0, ut1, ref, dis);
                        ut1 = 2 * BB::eva(b33, ref, dis, u.b[id2(f33,1,u.bsize)]) - uc0;
                    }
                    if (b12) {
                        BB::pre(m12, uw1, uc0, ref, dis);
                        uw1 = 2 * BB::eva(b12, ref, dis, u.b[id2(f12,1,u.bsize)]) - uc0;
                    }
                    if (b22) {
                        BB::pre(m22, us1, uc0, ref, dis);
                        us1 = 2 * BB::eva(b22, ref, dis, u.b[id2(f22,1,u.bsize)]) - uc0;
                    }
                    if (b32) {
                        BB::pre(m32, ub1, uc0, ref, dis);
                        ub1 = 2 * BB::eva(b32, ref, dis, u.b[id2(f32,1,u.bsize)]) - uc0;
                    }
                    if (b14) {
                        BB::pre(m14, ue1, ue2, ref, dis);
                        ue2 = 2 * BB::eva(b14, ref, dis, u.b[id2(f14,1,u.bsize)]) - ue1;
                    }
                    if (b24) {
                        BB::pre(m24, un1, un2, ref, dis);
                        un2 = 2 * BB::eva(b24, ref, dis, u.b[id2(f24,1,u.bsize)]) - un1;
                    }
                    if (b34) {
                        BB::pre(m34, ut1, ut2, ref, dis);
                        ut2 = 2 * BB::eva(b34, ref, dis, u.b[id2(f34,1,u.bsize)]) - ut1;
                    }
                    if (b11) {
                        BB::pre(m11, uw2, uw1, ref, dis);
                        uw2 = 2 * BB::eva(b11, ref, dis, u.b[id2(f11,1,u.bsize)]) - uw1;
                    }
                    if (b21) {
                        BB::pre(m21, us2, us1, ref, dis);
                        us2 = 2 * BB::eva(b21, ref, dis, u.b[id2(f21,1,u.bsize)]) - us1;
                    }
                    if (b31) {
                        BB::pre(m31, ub2, ub1, ref, dis);
                        ub2 = 2 * BB::eva(b31, ref, dis, u.b[id2(f31,1,u.bsize)]) - ub1;
                    }
                    ad2 = Scheme::ffvc_muscl(uc0, ue1, ue2, un1, un2, ut1, ut2, uw1, uw2, us1, us2, ub1, ub2, ufe, vfn, wft, ufw, vfs, wfb, det, epe, epn, ept, epw, eps, epb, kappa, beta, m1, m2);
                    vi2 = Scheme::viscosity(0, 0, 0, 0, 0, 0, mag, uc0, ue1, un1, ut1, uw1, us1, ub1, de1, dn1, dt1, dw1, ds1, db1, nc0, ne1, nn1, nt1, nw1, ns1, nb1, det, g1c, g2c, g3c, g1e, g2n, g3t, g1w, g2s, g3b, c.flow.ri);

                    uc0 = u3;
                    uw2 = u.m[id4(i-2,j  ,k  ,2,u.size)];
                    uw1 = u.m[id4(i-1,j  ,k  ,2,u.size)];
                    ue1 = u.m[id4(i+1,j  ,k  ,2,u.size)];
                    ue2 = u.m[id4(i+2,j  ,k  ,2,u.size)];
                    us2 = u.m[id4(i  ,j-2,k  ,2,u.size)];
                    us1 = u.m[id4(i  ,j-1,k  ,2,u.size)];
                    un1 = u.m[id4(i  ,j+1,k  ,2,u.size)];
                    un2 = u.m[id4(i  ,j+2,k  ,2,u.size)];
                    ub2 = u.m[id4(i  ,j  ,k-2,2,u.size)];
                    ub1 = u.m[id4(i  ,j  ,k-1,2,u.size)];
                    ut1 = u.m[id4(i  ,j  ,k+1,2,u.size)];
                    ut2 = u.m[id4(i  ,j  ,k+2,2,u.size)];
                    if (b13) {
                        BB::pre(m13, uc0, ue1, ref, dis);
                        ue1 = 2 * BB::eva(b13, ref, dis, u.b[id2(f13,2,u.bsize)]) - uc0;
                    }
                    if (b23) {
                        BB::pre(m23, uc0, un1, ref, dis);
                        un1 = 2 * BB::eva(b23, ref, dis, u.b[id2(f23,2,u.bsize)]) - uc0;
                    }
                    if (b33) {
                        BB::pre(m33, uc0, ut1, ref, dis);
                        ut1 = 2 * BB::eva(b33, ref, dis, u.b[id2(f33,2,u.bsize)]) - uc0;
                    }
                    if (b12) {
                        BB::pre(m12, uw1, uc0, ref, dis);
                        uw1 = 2 * BB::eva(b12, ref, dis, u.b[id2(f12,2,u.bsize)]) - uc0;
                    }
                    if (b22) {
                        BB::pre(m22, us1, uc0, ref, dis);
                        us1 = 2 * BB::eva(b22, ref, dis, u.b[id2(f22,2,u.bsize)]) - uc0;
                    }
                    if (b32) {
                        BB::pre(m32, ub1, uc0, ref, dis);
                        ub1 = 2 * BB::eva(b32, ref, dis, u.b[id2(f32,2,u.bsize)]) - uc0;
                    }
                    if (b14) {
                        BB::pre(m14, ue1, ue2, ref, dis);
                        ue2 = 2 * BB::eva(b14, ref, dis, u.b[id2(f14,2,u.bsize)]) - ue1;
                    }
                    if (b24) {
                        BB::pre(m24, un1, un2, ref, dis);
                        un2 = 2 * BB::eva(b24, ref, dis, u.b[id2(f24,2,u.bsize)]) - un1;
                    }
                    if (b34) {
                        BB::pre(m34, ut1, ut2, ref, dis);
                        ut2 = 2 * BB::eva(b34, ref, dis, u.b[id2(f34,2,u.bsize)]) - ut1;
                    }
                    if (b11) {
                        BB::pre(m11, uw2, uw1, ref, dis);
                        uw2 = 2 * BB::eva(b11, ref, dis, u.b[id2(f11,2,u.bsize)]) - uw1;
                    }
                    if (b21) {
                        BB::pre(m21, us2, us1, ref, dis);
                        us2 = 2 * BB::eva(b21, ref, dis, u.b[id2(f21,2,u.bsize)]) - us1;
                    }
                    if (b31) {
                        BB::pre(m31, ub2, ub1, ref, dis);
                        ub2 = 2 * BB::eva(b31, ref, dis, u.b[id2(f31,2,u.bsize)]) - ub1;
                    }
                    ad3 = Scheme::ffvc_muscl(uc0, ue1, ue2, un1, un2, ut1, ut2, uw1, uw2, us1, us2, ub1, ub2, ufe, vfn, wft, ufw, vfs, wfb, det, epe, epn, ept, epw, eps, epb, kappa, beta, m1, m2);
                    vi3 = Scheme::viscosity(0, 0, 0, 0, 0, 0, mag, uc0, ue1, un1, ut1, uw1, us1, ub1, de1, dn1, dt1, dw1, ds1, db1, nc0, ne1, nn1, nt1, nw1, ns1, nb1, det, g1c, g2c, g3c, g1e, g2n, g3t, g1w, g2s, g3b, c.flow.ri);

                    ua.m[id4(i,j,k,0,ua.size)] = u1 + c.time.dt * (-ad1 + vi1);
                    ua.m[id4(i,j,k,1,ua.size)] = u2 + c.time.dt * (-ad2 + vi2);
                    ua.m[id4(i,j,k,2,ua.size)] = u3 + c.time.dt * (-ad3 + vi3);
                }
            }
        }
    }
}

void MMAC::correct_center_velocity(vector_field<real_t> &u, vector_field<real_t> &ua, scalar_field<real_t> &p, Dom &dom) {
    scalar_field<unsigned> &f  = dom.f;
    vector_field<real_t>   &kx = dom.kx;
    Ctrl                   &c  = dom.c;

    #pragma acc kernels loop independent collapse(3) present(u, ua, p, f, kx, c, dom)
    for (int i = GUIDE; i < dom.size[0] - GUIDE; i ++) {
        for (int j = GUIDE; j < dom.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < dom.size[2] - GUIDE; k ++) {
                if (Util::ibsee(f.m[id3(i,j,k,f.size)], Flag::Active, Util::Mask1)) {
                    unsigned int b12, b13;
                    unsigned int b22, b23;
                    unsigned int b32, b33;
                    unsigned int f12, f13;
                    unsigned int f22, f23;
                    unsigned int f32, f33;
                    unsigned int m12, m13;
                    unsigned int m22, m23;
                    unsigned int m32, m33;
                    double kx1, kx2, kx3;
                    double dp1, dp2, dp3;
                    double pc0;
                    double pe1, pn1, pt1;
                    double pw1, ps1, pb1;
                    double ref, dis;

                    kx1 = kx.m[id4(i  ,j  ,k  ,0,kx.size)];
                    kx2 = kx.m[id4(i  ,j  ,k  ,1,kx.size)];
                    kx3 = kx.m[id4(i  ,j  ,k  ,2,kx.size)];
                    pc0 =  p.m[id3(i  ,j  ,k  ,   p.size)];
                    pe1 =  p.m[id3(i+1,j  ,k  ,   p.size)];
                    pn1 =  p.m[id3(i  ,j+1,k  ,   p.size)];
                    pt1 =  p.m[id3(i  ,j  ,k+1,   p.size)];
                    pw1 =  p.m[id3(i-1,j  ,k  ,   p.size)];
                    ps1 =  p.m[id3(i  ,j-1,k  ,   p.size)];
                    pb1 =  p.m[id3(i  ,j  ,k-1,   p.size)];
                    f12 =  f.m[id3(i-1,j  ,k  ,   f.size)];
                    f13 =  f.m[id3(i  ,j  ,k  ,   f.size)];
                    f22 =  f.m[id3(i  ,j-1,k  ,   f.size)];
                    f23 =  f.m[id3(i  ,j  ,k  ,   f.size)];
                    f32 =  f.m[id3(i  ,j  ,k-1,   f.size)];
                    f33 =  f.m[id3(i  ,j  ,k  ,   f.size)];
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

                    dp1 = 0.5 * (pe1 - pw1) * kx1;
                    dp2 = 0.5 * (pn1 - ps1) * kx2;
                    dp3 = 0.5 * (pt1 - pb1) * kx3;

                    u.m[id4(i,j,k,0,u.size)] = ua.m[id4(i,j,k,0,ua.size)] - c.time.dt * dp1;
                    u.m[id4(i,j,k,1,u.size)] = ua.m[id4(i,j,k,1,ua.size)] - c.time.dt * dp2;
                    u.m[id4(i,j,k,2,u.size)] = ua.m[id4(i,j,k,2,ua.size)] - c.time.dt * dp3;
                }
            }
        }
    }
}

void MMAC::correct_face_velocity(vector_field<real_t> &u, vector_field<real_t> &uu, vector_field<real_t> &uua, scalar_field<real_t> &p, Dom &dom) {
    scalar_field<unsigned> &f  = dom.f;
    vector_field<real_t>   &x  = dom.x;
    vector_field<real_t>   &g  = dom.g;
    scalar_field<real_t>   &ja = dom.ja;
    Ctrl                   &c  = dom.c;

    #pragma acc kernels loop independent collapse(3) present(u, uu, uua, p, f, x, g, ja, c, dom)
    for (int i = GUIDE - 1; i < dom.size[0] - GUIDE; i ++) {
        for (int j = GUIDE - 1; j < dom.size[1] - GUIDE; j ++) {
            for (int k = GUIDE - 1; k < dom.size[2] - GUIDE; k ++) {
                unsigned int ff;
                unsigned int f3, m3, b3;
                double u0, u1, uf;
                double p0, p1;
                double g0, g1;
                double kk;
                double ref, dis;
                ff = f.m[id3(i,j,k,f.size)];

                if (j >= GUIDE && k >= GUIDE) {
                    m3 = Util::ibsee(ff, Flag::Me, Util::Mask1);
                    f3 = Util::ibsee(ff, Flag::Fe, Util::Mask8);
                    if (!Util::ibsee(u.bflag[f3], BB::uu_locked, Util::Mask1)) {
                        b3 = Util::ibsee(u.bflag[f3], 0, Util::Mask8);
                        if (b3) {
                            u0 = u.m[id4(i  ,j,k,0,u.size)];
                            u1 = u.m[id4(i+1,j,k,0,u.size)];
                            BB::pre(m3, u0, u1, ref, dis);
                            uf = BB::eva(b3, ref, dis, u.b[id2(f3,0,u.bsize)]);
                            kk = 1 /    (x.m[id4(i+1,j,k,0, x.size)] -  x.m[id4(i,j,k,0, x.size)]);
                            uf = 0.5 * (ja.m[id3(i+1,j,k,  ja.size)] + ja.m[id3(i,j,k,  ja.size)]) * kk * uf;
                        } else {
                            p0 = p.m[id3(i  ,j,k,  p.size)];
                            p1 = p.m[id3(i+1,j,k,  p.size)];
                            g0 = g.m[id4(i  ,j,k,0,g.size)];
                            g1 = g.m[id4(i+1,j,k,0,g.size)];
                            b3 = Util::ibsee(p.bflag[f3], 0, Util::Mask8);
                            if (b3) {
                                BB::pre(m3, p0, p1, ref, dis);
                                p1 = 2 * BB::eva(b3, ref, dis, p.b[f3]) - p0;
                            }
                            uf = uua.m[id4(i,j,k,0,uua.size)] - c.time.dt * 0.5 * (g0 + g1) * (p1 - p0);
                        }
                        uu.m[id4(i,j,k,0,uua.size)] = uf;
                    }
                }

                if (i >= GUIDE && k >= GUIDE) {
                    m3 = Util::ibsee(ff, Flag::Mn, Util::Mask1);
                    f3 = Util::ibsee(ff, Flag::Fn, Util::Mask8);
                    if (!Util::ibsee(u.bflag[f3], BB::uu_locked, Util::Mask1)) {
                        b3 = Util::ibsee(u.bflag[f3], 0, Util::Mask8);
                        if (b3) {
                            u0 = u.m[id4(i,j  ,k,1,u.size)];
                            u1 = u.m[id4(i,j+1,k,1,u.size)];
                            BB::pre(m3, u0, u1, ref, dis);
                            uf = BB::eva(b3, ref, dis, u.b[id2(f3,1,u.bsize)]);
                            kk = 1 /    (x.m[id4(i,j+1,k,1, x.size)] -  x.m[id4(i,j,k,1, x.size)]);
                            uf = 0.5 * (ja.m[id3(i,j+1,k,  ja.size)] + ja.m[id3(i,j,k,  ja.size)]) * kk * uf;
                        } else {
                            p0 = p.m[id3(i,j  ,k,  p.size)];
                            p1 = p.m[id3(i,j+1,k,  p.size)];
                            g0 = g.m[id4(i,j  ,k,1,g.size)];
                            g1 = g.m[id4(i,j+1,k,1,g.size)];
                            b3 = Util::ibsee(p.bflag[f3], 0, Util::Mask8);
                            if (b3) {
                                BB::pre(m3, p0, p1, ref, dis);
                                p1 = 2 * BB::eva(b3, ref, dis, p.b[f3]) - p0;
                            }
                            uf = uua.m[id4(i,j,k,1,uua.size)] - c.time.dt * 0.5 * (g0 + g1) * (p1 - p0);
                        }
                        uu.m[id4(i,j,k,1,uua.size)] = uf;
                    }
                }

                if (i >= GUIDE && j >= GUIDE) {
                    m3 = Util::ibsee(ff, Flag::Mt, Util::Mask1);
                    f3 = Util::ibsee(ff, Flag::Ft, Util::Mask8);
                    if (!Util::ibsee(u.bflag[f3], BB::uu_locked, Util::Mask1)) {
                        b3 = Util::ibsee(u.bflag[f3], 0, Util::Mask8);
                        if (b3) {
                            u0 = u.m[id4(i,j,k  ,2,u.size)];
                            u1 = u.m[id4(i,j,k+1,2,u.size)];
                            BB::pre(m3, u0, u1, ref, dis);
                            uf = BB::eva(b3, ref, dis, u.b[id2(f3,2,u.bsize)]);
                            kk = 1 /    (x.m[id4(i,j,k+1,2, x.size)] -  x.m[id4(i,j,k,2, x.size)]);
                            uf = 0.5 * (ja.m[id3(i,j,k+1,  ja.size)] + ja.m[id3(i,j,k,  ja.size)]) * kk * uf;
                        } else {
                            p0 = p.m[id3(i,j,k  ,  p.size)];
                            p1 = p.m[id3(i,j,k+1,  p.size)];
                            g0 = g.m[id4(i,j,k  ,2,g.size)];
                            g1 = g.m[id4(i,j,k+1,2,g.size)];
                            b3 = Util::ibsee(p.bflag[f3], 0, Util::Mask8);
                            if (b3) {
                                BB::pre(m3, p0, p1, ref, dis);
                                p1 = 2 * BB::eva(b3, ref, dis, p.b[f3]) - p0;
                            }
                            uf = uua.m[id4(i,j,k,2,uua.size)] - c.time.dt * 0.5 * (g0 + g1) * (p1 - p0);
                        }
                        uu.m[id4(i,j,k,2,uua.size)] = uf;
                    }
                }
            }
        }
    }
}

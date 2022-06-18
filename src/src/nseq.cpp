#include <math.h>
#include "nseq.h"
#include "util.h"
#include "flag.h"

NSEq::NSEq() {
    this->scheme = NSEq::Scheme::muscl;
    this->alpha = 0.0;
}

void NSEq::pseudo_velocity(D *dom, BC *bc) {
    unsigned int *f;
    double *u, *ua, *uu, *sgs, *x, *ja, *g;
    BD *b;
    f   = dom->f;
    u   = dom->u;
    ua  = dom->ua;
    uu  = dom->uu;
    sgs = dom->sgs;
    x   = dom->x;
    ja  = dom->ja;
    g   = dom->g;
    b   =  bc->b;
    unsigned int scheme = this->scheme;
    unsigned int alpha = this->alpha;
    double dt = dom->dt;
    double ri = dom->ri;
    int size[4] = {dom->size[0], dom->size[1], dom->size[2], 3};

    for (int i = D::GUIDE; i < size[0] - D::GUIDE; i ++) {
        for (int j = D::GUIDE; j < size[1] - D::GUIDE; j ++) {
            for (int k = D::GUIDE; k < size[2] - D::GUIDE; k++) {
                if (Util::ibsee(f[id3(i,j,k,size)], Flag::Active, Util::Mask1)) {
                    double uc0;
                    double ue1, un1, ut1;
                    double ue2, un2, ut2;
                    double uw1, us1, ub1;
                    double uw2, us2, ub2;
                    double ufe, vfn, wft;
                    double ufw, vfs, wfb;
                    double ute, utn, utt;
                    double utw, uts, utb;
                    double de1, dn1, dt1;
                    double dw1, ds1, db1;
                    double nc0;
                    double ne1, nn1, nt1;
                    double nw1, ns1, nb1;
                    double det;
                    double ref, dis;
                    double ad1, ad2, ad3;
                    double vi1, vi2, vi3;
                    double xc0, yc0, zc0;
                    double xe1, yn1, zt1;
                    double xw1, ys1, zb1;
                    double g1c, g2c, g3c;
                    double g1e, g2n, g3t;
                    double g1w, g2s, g3b;
                    unsigned int f11, f12, f13, f14;
                    unsigned int f21, f22, f23, f24;
                    unsigned int f31, f32, f33, f34;
                    unsigned int b11, b12, b13, b14;
                    unsigned int b21, b22, b23, b24;
                    unsigned int b31, b32, b33, b34;
                    unsigned int m11, m12, m13, m14;
                    unsigned int m21, m22, m23, m24;
                    unsigned int m31, m32, m33, m34;
                    int epe, epn, ept;
                    int epw, eps, epb;
                    double u1, u2, u3, mag;

                    u1  = u[id4(i,j,k,0,size)];
                    u2  = u[id4(i,j,k,1,size)];
                    u3  = u[id4(i,j,k,2,size)];
                    mag = sqrt(u1 * u1 + u2 * u2 + u3 * u3);

                    ufe = uu[id4(i  ,j  ,k  ,0,size)];
                    vfn = uu[id4(i  ,j  ,k  ,1,size)];
                    wft = uu[id4(i  ,j  ,k  ,2,size)];
                    ufw = uu[id4(i-1,j  ,k  ,0,size)];
                    vfs = uu[id4(i  ,j-1,k  ,1,size)];
                    wfb = uu[id4(i  ,j  ,k-1,2,size)];
                    det = ja[id3(i  ,j  ,k  ,  size)];
                    g1c =  g[id4(i  ,j  ,k  ,0,size)];
                    g2c =  g[id4(i  ,j  ,k  ,1,size)];
                    g3c =  g[id4(i  ,j  ,k  ,2,size)];
                    g1e =  g[id4(i+1,j  ,k  ,0,size)];
                    g2n =  g[id4(i  ,j+1,k  ,1,size)];
                    g3t =  g[id4(i  ,j  ,k+1,2,size)];
                    g1w =  g[id4(i-1,j  ,k  ,0,size)];
                    g2s =  g[id4(i  ,j-1,k  ,1,size)];
                    g3b =  g[id4(i  ,j  ,k-1,2,size)];
                    f11 =  f[id3(i-2,j  ,j  ,  size)];
                    f12 =  f[id3(i-1,j  ,k  ,  size)];
                    f13 =  f[id3(i  ,j  ,k  ,  size)];
                    f14 =  f[id3(i+1,j  ,k  ,  size)];
                    f21 =  f[id3(i  ,j-2,j  ,  size)];
                    f22 =  f[id3(i  ,j-1,k  ,  size)];
                    f23 =  f[id3(i  ,j  ,k  ,  size)];
                    f24 =  f[id3(i  ,j+1,k  ,  size)];
                    f31 =  f[id3(i  ,j  ,j-2,  size)];
                    f32 =  f[id3(i  ,j  ,k-1,  size)];
                    f33 =  f[id3(i  ,j  ,k  ,  size)];
                    f34 =  f[id3(i  ,j  ,k+1,  size)];
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
                    epe = (f13)? 0 : 1;
                    epn = (f23)? 0 : 1;
                    ept = (f33)? 0 : 1;
                    epw = (f12)? 0 : 1;
                    eps = (f22)? 0 : 1;
                    epb = (f32)? 0 : 1;
                    de1 = 0.0;
                    dn1 = 0.0;
                    dt1 = 0.0;
                    dw1 = 0.0;
                    ds1 = 0.0;
                    db1 = 0.0;

                    nc0 = sgs[id3(i  ,j  ,k  ,size)];
                    ne1 = sgs[id3(i+1,j  ,k  ,size)];
                    nn1 = sgs[id3(i  ,j+1,k  ,size)];
                    nt1 = sgs[id3(i  ,j  ,k+1,size)];
                    nw1 = sgs[id3(i-1,j  ,k  ,size)];
                    ns1 = sgs[id3(i  ,j-1,k  ,size)];
                    nb1 = sgs[id3(i  ,j  ,k-1,size)];
                    b12 = Util::ibsee(b[f12].flag, BD::Nb, Util::Mask2);
                    b13 = Util::ibsee(b[f13].flag, BD::Nb, Util::Mask2);
                    b22 = Util::ibsee(b[f22].flag, BD::Nb, Util::Mask2);
                    b23 = Util::ibsee(b[f23].flag, BD::Nb, Util::Mask2);
                    b32 = Util::ibsee(b[f32].flag, BD::Nb, Util::Mask2);
                    b33 = Util::ibsee(b[f33].flag, BD::Nb, Util::Mask2);
                    if (b13) {
                        BC::pre(m13, nc0, ne1, ref, dis);
                        ne1 = 2 * BC::eva(b13, ref, dis, b[f12].nut) - nc0;
                    }
                    if (b23) {
                        BC::pre(m23, nc0, nn1, ref, dis);
                        nn1 = 2 * BC::eva(b23, ref, dis, b[f23].nut) - nc0;
                    }
                    if (b33) {
                        BC::pre(m33, nc0, nt1, ref, dis);
                        nt1 = 2 * BC::eva(b33, ref, dis, b[f33].nut) - nc0;
                    }
                    if (b12) {
                        BC::pre(m12, nw1, nc0, ref, dis);
                        nw1 = 2 * BC::eva(b12, ref, dis, b[f12].nut) - nc0;
                    }
                    if (b22) {
                        BC::pre(m22, ns1, nc0, ref, dis);
                        ns1 = 2 * BC::eva(b22, ref, dis, b[f22].nut) - nc0;
                    }
                    if (b32) {
                        BC::pre(m32, nb1, nc0, ref, dis);
                        nb1 = 2 * BC::eva(b32, ref, dis, b[f32].nut) - nc0;
                    }

                    b11 = Util::ibsee(b[f11].flag, BD::Ub, Util::Mask2);
                    b12 = Util::ibsee(b[f12].flag, BD::Ub, Util::Mask2);
                    b13 = Util::ibsee(b[f13].flag, BD::Ub, Util::Mask2);
                    b14 = Util::ibsee(b[f14].flag, BD::Ub, Util::Mask2);
                    b21 = Util::ibsee(b[f21].flag, BD::Ub, Util::Mask2);
                    b22 = Util::ibsee(b[f22].flag, BD::Ub, Util::Mask2);
                    b23 = Util::ibsee(b[f23].flag, BD::Ub, Util::Mask2);
                    b24 = Util::ibsee(b[f24].flag, BD::Ub, Util::Mask2);
                    b31 = Util::ibsee(b[f31].flag, BD::Ub, Util::Mask2);
                    b32 = Util::ibsee(b[f32].flag, BD::Ub, Util::Mask2);
                    b33 = Util::ibsee(b[f33].flag, BD::Ub, Util::Mask2);
                    b34 = Util::ibsee(b[f34].flag, BD::Ub, Util::Mask2);

                    uc0 = u1;
                    uw2 = u[id4(i-2,j  ,k  ,0,size)];
                    uw1 = u[id4(i-1,j  ,k  ,0,size)];
                    ue1 = u[id4(i+1,j  ,k  ,0,size)];
                    ue2 = u[id4(i+2,j  ,k  ,0,size)];
                    us2 = u[id4(i  ,j-2,k  ,0,size)];
                    us1 = u[id4(i  ,j-1,k  ,0,size)];
                    un1 = u[id4(i  ,j+1,k  ,0,size)];
                    un2 = u[id4(i  ,j+2,k  ,0,size)];
                    ub2 = u[id4(i  ,j  ,k-2,0,size)];
                    ub1 = u[id4(i  ,j  ,k-1,0,size)];
                    ut1 = u[id4(i  ,j  ,k+1,0,size)];
                    ut2 = u[id4(i  ,j  ,k+2,0,size)];
                    if (b13) {
                        BC::pre(m13, uc0, ue1, ref, dis);
                        ue1 = 2 * BC::eva(b13, ref, dis, b[f12].u) - uc0;
                    }
                    if (b23) {
                        BC::pre(m23, uc0, un1, ref, dis);
                        un1 = 2 * BC::eva(b23, ref, dis, b[f23].u) - uc0;
                    }
                    if (b33) {
                        BC::pre(m33, uc0, ut1, ref, dis);
                        ut1 = 2 * BC::eva(b33, ref, dis, b[f33].u) - uc0;
                    }
                    if (b12) {
                        BC::pre(m12, uw1, uc0, ref, dis);
                        uw1 = 2 * BC::eva(b12, ref, dis, b[f12].u) - uc0;
                    }
                    if (b22) {
                        BC::pre(m22, us1, uc0, ref, dis);
                        us1 = 2 * BC::eva(b22, ref, dis, b[f22].u) - uc0;
                    }
                    if (b32) {
                        BC::pre(m32, ub1, uc0, ref, dis);
                        ub1 = 2 * BC::eva(b32, ref, dis, b[f32].u) - uc0;
                    }
                    ad1 = NSEq::ffvc_muscl(uc0, ue1, ue2, un1, un2, ut1, ut2, uw1, uw2, us1, us2, ub1, ub2, ufe, vfn, wft, ufw, vfs, wfb, det, epe, epn, ept, epw, eps, epb);
                    vi1 = NSEq::viscosity(f13, f23, f33, f12, f22, f32, mag, uc0, ue1, un1, ut1, uw1, us1, ub1, ute, utn, utt, utw, uts, utb, de1, dn1, dt1, dw1, ds1, db1, nc0, ne1, nn1, nt1, nw1, ns1, nb1, det, g1c, g2c, g3c, g1e, g2n, g3t, g1w, g2s, g3b, ri);

                    uc0 = u2;
                    uw2 = u[id4(i-2,j  ,k  ,1,size)];
                    uw1 = u[id4(i-1,j  ,k  ,1,size)];
                    ue1 = u[id4(i+1,j  ,k  ,1,size)];
                    ue2 = u[id4(i+2,j  ,k  ,1,size)];
                    us2 = u[id4(i  ,j-2,k  ,1,size)];
                    us1 = u[id4(i  ,j-1,k  ,1,size)];
                    un1 = u[id4(i  ,j+1,k  ,1,size)];
                    un2 = u[id4(i  ,j+2,k  ,1,size)];
                    ub2 = u[id4(i  ,j  ,k-2,1,size)];
                    ub1 = u[id4(i  ,j  ,k-1,1,size)];
                    ut1 = u[id4(i  ,j  ,k+1,1,size)];
                    ut2 = u[id4(i  ,j  ,k+2,1,size)];
                    if (b13) {
                        BC::pre(m13, uc0, ue1, ref, dis);
                        ue1 = 2 * BC::eva(b13, ref, dis, b[f12].v) - uc0;
                    }
                    if (b23) {
                        BC::pre(m23, uc0, un1, ref, dis);
                        un1 = 2 * BC::eva(b23, ref, dis, b[f23].v) - uc0;
                    }
                    if (b33) {
                        BC::pre(m33, uc0, ut1, ref, dis);
                        ut1 = 2 * BC::eva(b33, ref, dis, b[f33].v) - uc0;
                    }
                    if (b12) {
                        BC::pre(m12, uw1, uc0, ref, dis);
                        uw1 = 2 * BC::eva(b12, ref, dis, b[f12].v) - uc0;
                    }
                    if (b22) {
                        BC::pre(m22, us1, uc0, ref, dis);
                        us1 = 2 * BC::eva(b22, ref, dis, b[f22].v) - uc0;
                    }
                    if (b32) {
                        BC::pre(m32, ub1, uc0, ref, dis);
                        ub1 = 2 * BC::eva(b32, ref, dis, b[f32].v) - uc0;
                    }
                    ad2 = NSEq::ffvc_muscl(uc0, ue1, ue2, un1, un2, ut1, ut2, uw1, uw2, us1, us2, ub1, ub2, ufe, vfn, wft, ufw, vfs, wfb, det, epe, epn, ept, epw, eps, epb);
                    vi2 = NSEq::viscosity(f13, f23, f33, f12, f22, f32, mag, uc0, ue1, un1, ut1, uw1, us1, ub1, ute, utn, utt, utw, uts, utb, de1, dn1, dt1, dw1, ds1, db1, nc0, ne1, nn1, nt1, nw1, ns1, nb1, det, g1c, g2c, g3c, g1e, g2n, g3t, g1w, g2s, g3b, ri);

                    uc0 = u3;
                    uw2 = u[id4(i-2,j  ,k  ,2,size)];
                    uw1 = u[id4(i-1,j  ,k  ,2,size)];
                    ue1 = u[id4(i+1,j  ,k  ,2,size)];
                    ue2 = u[id4(i+2,j  ,k  ,2,size)];
                    us2 = u[id4(i  ,j-2,k  ,2,size)];
                    us1 = u[id4(i  ,j-1,k  ,2,size)];
                    un1 = u[id4(i  ,j+1,k  ,2,size)];
                    un2 = u[id4(i  ,j+2,k  ,2,size)];
                    ub2 = u[id4(i  ,j  ,k-2,2,size)];
                    ub1 = u[id4(i  ,j  ,k-1,2,size)];
                    ut1 = u[id4(i  ,j  ,k+1,2,size)];
                    ut2 = u[id4(i  ,j  ,k+2,2,size)];
                    if (b13) {
                        BC::pre(m13, uc0, ue1, ref, dis);
                        ue1 = 2 * BC::eva(b13, ref, dis, b[f12].w) - uc0;
                    }
                    if (b23) {
                        BC::pre(m23, uc0, un1, ref, dis);
                        un1 = 2 * BC::eva(b23, ref, dis, b[f23].w) - uc0;
                    }
                    if (b33) {
                        BC::pre(m33, uc0, ut1, ref, dis);
                        ut1 = 2 * BC::eva(b33, ref, dis, b[f33].w) - uc0;
                    }
                    if (b12) {
                        BC::pre(m12, uw1, uc0, ref, dis);
                        uw1 = 2 * BC::eva(b12, ref, dis, b[f12].w) - uc0;
                    }
                    if (b22) {
                        BC::pre(m22, us1, uc0, ref, dis);
                        us1 = 2 * BC::eva(b22, ref, dis, b[f22].w) - uc0;
                    }
                    if (b32) {
                        BC::pre(m32, ub1, uc0, ref, dis);
                        ub1 = 2 * BC::eva(b32, ref, dis, b[f32].w) - uc0;
                    }
                    ad3 = NSEq::ffvc_muscl(uc0, ue1, ue2, un1, un2, ut1, ut2, uw1, uw2, us1, us2, ub1, ub2, ufe, vfn, wft, ufw, vfs, wfb, det, epe, epn, ept, epw, eps, epb);
                    vi3 = NSEq::viscosity(f13, f23, f33, f12, f22, f32, mag, uc0, ue1, un1, ut1, uw1, us1, ub1, ute, utn, utt, utw, uts, utb, de1, dn1, dt1, dw1, ds1, db1, nc0, ne1, nn1, nt1, nw1, ns1, nb1, det, g1c, g2c, g3c, g1e, g2n, g3t, g1w, g2s, g3b, ri);

                    ua[id4(i,j,k,0,size)] = u1 - dt * (-ad1 + vi1);
                    ua[id4(i,j,k,1,size)] = u2 - dt * (-ad2 + vi2);
                    ua[id4(i,j,k,2,size)] = u3 - dt * (-ad3 + vi3);
                }
            }
        }
    }
}

void NSEq::contra_pseudo_velocity(D *dom, BC *bc) {
    unsigned int *f;
    double *ua, *uc, *uu, *x, *kx, *ja;
    BD *b;
    f  = dom->f;
    ua = dom->ua;
    uc = dom->uc;
    uu = dom->uu;
    x  = dom->x;
    kx = dom->kx;
    ja = dom->ja;
    b  =  bc->b;
    int size[4] = {dom->size[0], dom->size[1], dom->size[2], 3};

    for (int i = D::GUIDE - 1; i <= size[0] - D::GUIDE; i ++) {
        for (int j = D::GUIDE - 1; j <= size[1] - D::GUIDE; j ++) {
            for (int k = D::GUIDE - 1; k <= size[2] - D::GUIDE; k ++) {
                double k1, k2, k3, u1, u2, u3, de;

                u1 = ua[id4(i,j,k,0,size)];
                u2 = ua[id4(i,j,k,1,size)];
                u3 = ua[id4(i,j,k,2,size)];
                k1 = kx[id4(i,j,k,0,size)];
                k2 = kx[id4(i,j,k,1,size)];
                k3 = kx[id4(i,j,k,2,size)];
                de = ja[id3(i,j,k,  size)];

                uc[id4(i,j,k,0,size)] = de * k1 * u1;
                uc[id4(i,j,k,1,size)] = de * k2 * u2;
                uc[id4(i,j,k,2,size)] = de * k3 * u3;
            }
        }
    }

    for (int i = D::GUIDE - 1; i < size[0] - D::GUIDE; i ++) {
        for (int j = D::GUIDE - 1; j < size[1] - D::GUIDE; j ++) {
            for (int k = D::GUIDE - 1; k < size[2] - D::GUIDE; k ++) {
                unsigned int ff;
                unsigned int f3, m3, b3;
                double u0, u1, uf;
                double kk;
                double ref, dis;
                ff = f[id3(i,j,k,size)];

                m3 = Util::ibsee(ff, Flag::Me, Util::Mask1);
                f3 = Util::ibsee(ff, Flag::Fe, Util::Mask8);
                b3 = Util::ibsee(b[f3].flag, BD::Ub, Util::Mask2);
                if (b3) {
                    u0 = ua[id4(i  ,j,k,0,size)];
                    u1 = ua[id4(i+1,j,k,0,size)];
                    BC::pre(m3, u0, u1, ref, dis);
                    uf = BC::eva(b3, ref, dis, b[f3].u);
                    kk = 1 / (x[id4(i+1,j,k,0,size)] - x[id4(i,j,k,0,size)]);
                    uf = 0.5 * (ja[id3(i+1,j,k,size)] + ja[id3(i,j,k,size)]) * kk * uf;
                } else {
                    u0 = uc[id4(i  ,j,k,0,size)];
                    u1 = uc[id4(i+1,j,k,0,size)];
                    uf = 0.5 * (u0 + u1);
                }
                uu[id4(i,j,k,0,size)] = uf;

                m3 = Util::ibsee(ff, Flag::Mn, Util::Mask1);
                f3 = Util::ibsee(ff, Flag::Fn, Util::Mask8);
                b3 = Util::ibsee(b[f3].flag, BD::Ub, Util::Mask2);
                if (b3) {
                    u0 = ua[id4(i,j  ,k,1,size)];
                    u1 = ua[id4(i,j+1,k,1,size)];
                    BC::pre(m3, u0, u1, ref, dis);
                    uf = BC::eva(b3, ref, dis, b[f3].u);
                    kk = 1 / (x[id4(i,j+1,k,1,size)] - x[id4(i,j,k,1,size)]);
                    uf = 0.5 * (ja[id3(i,j+1,k,size)] + ja[id3(i,j,k,size)]) * kk * uf;
                } else {
                    u0 = uc[id4(i,j  ,k,1,size)];
                    u1 = uc[id4(i,j+1,k,1,size)];
                    uf = 0.5 * (u0 + u1);
                }
                uu[id4(i,j,k,1,size)] = uf;

                m3 = Util::ibsee(ff, Flag::Mt, Util::Mask1);
                f3 = Util::ibsee(ff, Flag::Ft, Util::Mask8);
                b3 = Util::ibsee(b[f3].flag, BD::Ub, Util::Mask2);
                if (b3) {
                    u0 = ua[id4(i,j,k  ,2,size)];
                    u1 = ua[id4(i,j,k+1,2,size)];
                    BC::pre(m3, u0, u1, ref, dis);
                    uf = BC::eva(b3, ref, dis, b[f3].u);
                    kk = 1 / (x[id4(i,j,k+1,2,size)] - x[id4(i,j,k,2,size)]);
                    uf = 0.5 * (ja[id3(i,j,k+1,size)] + ja[id3(i,j,k,size)]) * kk * uf;
                } else {
                    u0 = uc[id4(i,j,k  ,2,size)];
                    u1 = uc[id4(i,j,k+1,2,size)];
                    uf = 0.5 * (u0 + u1);
                }
                uu[id4(i,j,k,2,size)] = uf;
            }
        }
    }
}

void NSEq::correct_ceneter_velocity(D *dom, BC *bc) {
    unsigned int *f;
    double *u, *ua, *p, *kx;
    BD *b;
    f  = dom->f;
    u  = dom->u;
    ua = dom->ua;
    p  = dom->p;
    kx = dom->kx;
    b  =  bc->b;
    int size[4] = {dom->size[0], dom->size[1], dom->size[2], 3};
    double dt = dom->dt;

    for (int i = D::GUIDE; i < size[0] - D::GUIDE; i ++) {
        for (int j = D::GUIDE; j < size[1] - D::GUIDE; j ++) {
            for (int k = D::GUIDE; k < size[2] - D::GUIDE; k ++) {
                if (Util::ibsee(f[id3(i,j,k,size)], Flag::Active, Util::Mask1)) {
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

                    kx1 = kx[id4(i  ,j  ,k  ,0,size)];
                    kx2 = kx[id4(i  ,j  ,k  ,1,size)];
                    kx3 = kx[id4(i  ,j  ,k  ,2,size)];
                    pc0 =  p[id3(i  ,j  ,k  ,  size)];
                    pe1 =  p[id3(i+1,j  ,k  ,  size)];
                    pn1 =  p[id3(i  ,j+1,k  ,  size)];
                    pt1 =  p[id3(i  ,j  ,k+1,  size)];
                    pw1 =  p[id3(i-1,j  ,k  ,  size)];
                    ps1 =  p[id3(i  ,j-1,k  ,  size)];
                    pb1 =  p[id3(i  ,j  ,k-1,  size)];
                    f12 =  f[id3(i-1,j  ,k  ,  size)];
                    f13 =  f[id3(i  ,j  ,k  ,  size)];
                    f22 =  f[id3(i  ,j-1,k  ,  size)];
                    f23 =  f[id3(i  ,j  ,k  ,  size)];
                    f32 =  f[id3(i  ,j  ,k-1,  size)];
                    f33 =  f[id3(i  ,j  ,k  ,  size)];
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
                    b12 = Util::ibsee(b[f12].flag, BD::Pb, Util::Mask2);
                    b13 = Util::ibsee(b[f13].flag, BD::Pb, Util::Mask2);
                    b22 = Util::ibsee(b[f22].flag, BD::Pb, Util::Mask2);
                    b23 = Util::ibsee(b[f23].flag, BD::Pb, Util::Mask2);
                    b32 = Util::ibsee(b[f32].flag, BD::Pb, Util::Mask2);
                    b33 = Util::ibsee(b[f33].flag, BD::Pb, Util::Mask2);
                    if (b13) {
                        BC::pre(m13, pc0, pe1, ref, dis);
                        pe1 = 2 * BC::eva(b13, ref, dis, b[f13].p) - pc0;
                    }
                    if (b23) {
                        BC::pre(m23, pc0, pn1, ref, dis);
                        pn1 = 2 * BC::eva(b23, ref, dis, b[f23].p) - pc0;
                    }
                    if (b33) {
                        BC::pre(m33, pc0, pt1, ref, dis);
                        pt1 = 2 * BC::eva(b33, ref, dis, b[f33].p) - pc0;
                    }
                    if (b12) {
                        BC::pre(m12, pw1, pc0, ref, dis);
                        pw1 = 2 * BC::eva(b12, ref, dis, b[f12].p) - pc0;
                    }
                    if (b22) {
                        BC::pre(m22, ps1, pc0, ref, dis);
                        ps1 = 2 * BC::eva(b22, ref, dis, b[f22].p) - pc0;
                    }
                    if (b32) {
                        BC::pre(m32, pb1, pc0, ref, dis);
                        pb1 = 2 * BC::eva(b32, ref, dis, b[f32].p) - pc0;
                    }

                    dp1 = 0.5 * (pe1 - pw1) * kx1;
                    dp2 = 0.5 * (pn1 - ps1) * kx2;
                    dp3 = 0.5 * (pt1 - pb1) * kx3;

                    u[id4(i,j,k,0,size)] = ua[id4(i,j,k,0,size)] - dt * dp1;
                    u[id4(i,j,k,1,size)] = ua[id4(i,j,k,1,size)] - dt * dp2;
                    u[id4(i,j,k,2,size)] = ua[id4(i,j,k,2,size)] - dt * dp3;
                }
            }
        }
    }
}

void NSEq::correct_face_velocity(D *dom, BC *bc) {
    unsigned int *f;
    double *u, *uu, *uua, *p, *x, *g, *ja;
    BD *b;
    f   = dom->f;
    u   = dom->u;
    uu  = dom->uu;
    uua = dom->uua;
    p   = dom->p;
    x   = dom->x;
    g   = dom->g;
    ja  = dom->ja;
    b   =  bc->b;
    int size[4] = {dom->size[0], dom->size[1], dom->size[2], 3};
    double dt = dom->dt;

    for (int i = D::GUIDE - 1; i < size[0] - D::GUIDE; i ++) {
        for (int j = D::GUIDE - 1; j < size[1] - D::GUIDE; j ++) {
            for (int k = D::GUIDE - 1; k < size[2] - D::GUIDE; k ++) {
                unsigned int ff;
                unsigned int f3, m3, b3;
                double u0, u1, uf;
                double p0, p1;
                double g0, g1;
                double kk;
                double ref, dis;
                ff = f[id3(i,j,k,size)];

                m3 = Util::ibsee(ff, Flag::Me, Util::Mask1);
                f3 = Util::ibsee(ff, Flag::Fe, Util::Mask8);
                b3 = Util::ibsee(b[f3].flag, BD::Ub, Util::Mask2);
                if (b3) {
                    u0 = u[id4(i  ,j,k,0,size)];
                    u1 = u[id4(i+1,j,k,0,size)];
                    BC::pre(m3, u0, u1, ref, dis);
                    uf = BC::eva(b3, ref, dis, b[f3].u);
                    kk = 1 / (x[id4(i+1,j,k,0,size)] - x[id4(i,j,k,0,size)]);
                    uf = 0.5 * (ja[id3(i+1,j,k,size)] + ja[id3(i,j,k,size)]) * kk * uf;
                } else {
                    p0 = p[id3(i  ,j,k,  size)];
                    p1 = p[id3(i+1,j,k,  size)];
                    g0 = g[id4(i  ,j,k,0,size)];
                    g1 = g[id4(i+1,j,k,0,size)];
                    uf = uua[id4(i,j,k,0,size)] - dt * 0.5 * (g0 + g1) * (p1 - p0);
                }
                uu[id4(i,j,k,0,size)] = uf;

                m3 = Util::ibsee(ff, Flag::Mn, Util::Mask1);
                f3 = Util::ibsee(ff, Flag::Fn, Util::Mask8);
                b3 = Util::ibsee(b[f3].flag, BD::Ub, Util::Mask2);
                if (b3) {
                    u0 = u[id4(i,j  ,k,1,size)];
                    u1 = u[id4(i,j+1,k,1,size)];
                    BC::pre(m3, u0, u1, ref, dis);
                    uf = BC::eva(b3, ref, dis, b[f3].u);
                    kk = 1 / (x[id4(i,j+1,k,1,size)] - x[id4(i,j,k,1,size)]);
                    uf = 0.5 * (ja[id3(i,j+1,k,size)] + ja[id3(i,j,k,size)]) * kk * uf;
                } else {
                    p0 = p[id3(i,j  ,k,  size)];
                    p1 = p[id3(i,j+1,k,  size)];
                    g0 = g[id4(i,j  ,k,1,size)];
                    g1 = g[id4(i,j+1,k,1,size)];
                    uf = uua[id4(i,j,k,1,size)] - dt * 0.5 * (g0 + g1) * (p1 - p0);
                }
                uu[id4(i,j,k,1,size)] = uf;

                m3 = Util::ibsee(ff, Flag::Mt, Util::Mask1);
                f3 = Util::ibsee(ff, Flag::Ft, Util::Mask8);
                b3 = Util::ibsee(b[f3].flag, BD::Ub, Util::Mask2);
                if (b3) {
                    u0 = u[id4(i,j,k  ,2,size)];
                    u1 = u[id4(i,j,k+1,2,size)];
                    BC::pre(m3, u0, u1, ref, dis);
                    uf = BC::eva(b3, ref, dis, b[f3].u);
                    kk = 1 / (x[id4(i,j,k+1,2,size)] - x[id4(i,j,k,2,size)]);
                    uf = 0.5 * (ja[id3(i,j,k+1,size)] + ja[id3(i,j,k,size)]) * kk * uf;
                } else {
                    p0 = p[id3(i,j,k  ,  size)];
                    p1 = p[id3(i,j,k+1,  size)];
                    g0 = g[id4(i,j,k  ,2,size)];
                    g1 = g[id4(i,j,k+1,2,size)];
                    uf = uua[id4(i,j,k,2,size)] - dt * 0.5 * (g0 + g1) * (p1 - p0);
                }
                uu[id4(i,j,k,2,size)] = uf;
            }
        }
    }
}

double NSEq::ffvc_muscl(
    double uc0,
    double ue1,
    double ue2,
    double un1,
    double un2,
    double ut1,
    double ut2,
    double uw1,
    double uw2,
    double us1,
    double us2,
    double ub1,
    double ub2,
    double ufe,
    double vfn,
    double wft,
    double ufw,
    double vfs,
    double wfb,
    double det,
    int    epe,
    int    epn,
    int    ept,
    int    epw,
    int    eps,
    int    epb
) {
    double d1, d2, d3, d4;
    double s1, s2, s3, s4;
    double g1, g2, g3, g4, g5, g6;
    double u11, u10, u01, u00;
    double f1, f0;
    double adv;

    double kappa = 1.0 / 3.0;
    double beta  = (3.0 - kappa) / (1.0 - kappa);
    double m1    = 1.0 - kappa;
    double m2    = 1.0 + kappa;

    d4  = ue2 - ue1;
    d3  = ue1 - uc0;
    d2  = uc0 - uw1;
    d1  = uw1 - uw2;
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
    adv += (f1 - f0) / det;

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
    adv += (f1 - f0) / det;

    return adv;
}

double NSEq::viscosity(
    unsigned int f13,
    unsigned int f23,
    unsigned int f33,
    unsigned int f12,
    unsigned int f22,
    unsigned int f32,
    double       mag,
    double       uc0,
    double       ue1,
    double       un1,
    double       ut1,
    double       uw1,
    double       us1,
    double       ub1,
    double       ute,
    double       utn,
    double       utt,
    double       utw,
    double       uts,
    double       utb,
    double       de1,
    double       dn1,
    double       dt1,
    double       dw1,
    double       ds1,
    double       db1,
    double       nc0,
    double       ne1,
    double       nn1,
    double       nt1,
    double       nw1,
    double       ns1,
    double       nb1,
    double       det,
    double       g1c,
    double       g2c,
    double       g3c,
    double       g1e,
    double       g2n,
    double       g3t,
    double       g1w,
    double       g2s,
    double       g3b,
    double       ri
) {
    double ffe, ffn, fft, ffw, ffs, ffb;
    ffe = (ri + 0.5 * (nc0 + ne1)) * 0.5 * (g1e + g1c) * (ue1 - uc0);
    ffn = (ri + 0.5 * (nc0 + nn1)) * 0.5 * (g2n + g2c) * (un1 - uc0);
    fft = (ri + 0.5 * (nc0 + nt1)) * 0.5 * (g3t + g3c) * (ut1 - uc0);
    ffw = (ri + 0.5 * (nc0 + nw1)) * 0.5 * (g1c + g1w) * (uc0 - uw1);
    ffs = (ri + 0.5 * (nc0 + ns1)) * 0.5 * (g2c + g2s) * (uc0 - us1);
    ffb = (ri + 0.5 * (nc0 + nb1)) * 0.5 * (g3c + g3b) * (uc0 - ub1);

    return (ffe + ffn + fft - ffw - ffs - ffb) / det;
}
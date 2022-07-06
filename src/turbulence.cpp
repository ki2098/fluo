#include <math.h>
#include "../include/turbulence.h"
#include "../include/flag.h"
#include "../include/boundary.h"

void Turbulence::smagorinsky(Dom &dom) {
    Ctrl                   &c   = dom.c;
    scalar_field<unsigned> &f   = dom.f;
    vector_field<real_t>   &u   = dom.u;
    vector_field<real_t>   &kx  = dom.kx;
    scalar_field<real_t>   &ja  = dom.ja;
    scalar_field<real_t>   &nut = dom.nut;

    #pragma acc kernels loop independent collapse(3) present(f, u, kx, ja, nut, c, dom)
    for (int i = GUIDE; i < dom.size[0] - GUIDE; i ++) {
        for (int j = GUIDE; j < dom.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < dom.size[2] - GUIDE; k ++) {
                if (Util::ibsee(f.m[id3(i,j,k,f.size)], Cell::Active, Util::Mask1)) {
                    unsigned b12, b13;
                    unsigned b22, b23;
                    unsigned b32, b33;
                    unsigned f12, f13;
                    unsigned f22, f23;
                    unsigned f32, f33;
                    unsigned m12, m13;
                    unsigned m22, m23;
                    unsigned m32, m33;
                    real_t   kx1, kx2, kx3;
                    real_t   uc0;
                    real_t   ue1, un1, ut1;
                    real_t   uw1, us1, ub1;
                    real_t   ref, dis;
                    real_t   dux, duy, duz;
                    real_t   dvx, dvy, dvz;
                    real_t   dwx, dwy, dwz;
                    real_t   det;
                    real_t   d01, d02, d03, d04, d05, d06;
                    real_t   Del, D_u, l_0;

                    kx1 = kx.m[id4(i  ,j  ,k  ,0,kx.size)];
                    kx2 = kx.m[id4(i  ,j  ,k  ,1,kx.size)];
                    kx3 = kx.m[id4(i  ,j  ,k  ,2,kx.size)];
                    det = ja.m[id3(i  ,j  ,k  ,  ja.size)];
                    f12 =  f.m[id3(i-1,j  ,k  ,   f.size)];
                    f13 =  f.m[id3(i  ,j  ,k  ,   f.size)];
                    f22 =  f.m[id3(i  ,j-1,k  ,   f.size)];
                    f23 =  f.m[id3(i  ,j  ,k  ,   f.size)];
                    f32 =  f.m[id3(i  ,j  ,k-1,   f.size)];
                    f33 =  f.m[id3(i  ,j  ,k  ,   f.size)];
                    m12 = Util::ibsee(f12, Cell::Me, Util::Mask1);
                    m13 = Util::ibsee(f13, Cell::Me, Util::Mask1);
                    m22 = Util::ibsee(f22, Cell::Mn, Util::Mask1);
                    m23 = Util::ibsee(f23, Cell::Mn, Util::Mask1);
                    m32 = Util::ibsee(f32, Cell::Mt, Util::Mask1);
                    m33 = Util::ibsee(f33, Cell::Mt, Util::Mask1);
                    f12 = Util::ibsee(f12, Cell::Fe, Util::Mask8);
                    f13 = Util::ibsee(f13, Cell::Fe, Util::Mask8);
                    f22 = Util::ibsee(f22, Cell::Fn, Util::Mask8);
                    f23 = Util::ibsee(f23, Cell::Fn, Util::Mask8);
                    f32 = Util::ibsee(f32, Cell::Ft, Util::Mask8);
                    f33 = Util::ibsee(f33, Cell::Ft, Util::Mask8);
                    b12 = Util::ibsee(u.bflag[f12], 0, Util::Mask8);
                    b13 = Util::ibsee(u.bflag[f13], 0, Util::Mask8);
                    b22 = Util::ibsee(u.bflag[f22], 0, Util::Mask8);
                    b23 = Util::ibsee(u.bflag[f23], 0, Util::Mask8);
                    b32 = Util::ibsee(u.bflag[f32], 0, Util::Mask8);
                    b33 = Util::ibsee(u.bflag[f33], 0, Util::Mask8);

                    uc0 = u.m[id4(i  ,j  ,k  ,0,u.size)];
                    uw1 = u.m[id4(i-1,j  ,k  ,0,u.size)];
                    ue1 = u.m[id4(i+1,j  ,k  ,0,u.size)];
                    us1 = u.m[id4(i  ,j-1,k  ,0,u.size)];
                    un1 = u.m[id4(i  ,j+1,k  ,0,u.size)];
                    ub1 = u.m[id4(i  ,j  ,k-1,0,u.size)];
                    ut1 = u.m[id4(i  ,j  ,k+1,0,u.size)];
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
                    dux = kx1 * 0.5 * (ue1 - uw1);
                    duy = kx2 * 0.5 * (un1 - us1);
                    duz = kx3 * 0.5 * (ut1 - ub1);

                    uc0 = u.m[id4(i  ,j  ,k  ,1,u.size)];
                    uw1 = u.m[id4(i-1,j  ,k  ,1,u.size)];
                    ue1 = u.m[id4(i+1,j  ,k  ,1,u.size)];
                    us1 = u.m[id4(i  ,j-1,k  ,1,u.size)];
                    un1 = u.m[id4(i  ,j+1,k  ,1,u.size)];
                    ub1 = u.m[id4(i  ,j  ,k-1,1,u.size)];
                    ut1 = u.m[id4(i  ,j  ,k+1,1,u.size)];
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
                    dvx = kx1 * 0.5 * (ue1 - uw1);
                    dvy = kx2 * 0.5 * (un1 - us1);
                    dvz = kx3 * 0.5 * (ut1 - ub1);

                    uc0 = u.m[id4(i  ,j  ,k  ,2,u.size)];
                    uw1 = u.m[id4(i-1,j  ,k  ,2,u.size)];
                    ue1 = u.m[id4(i+1,j  ,k  ,2,u.size)];
                    us1 = u.m[id4(i  ,j-1,k  ,2,u.size)];
                    un1 = u.m[id4(i  ,j+1,k  ,2,u.size)];
                    ub1 = u.m[id4(i  ,j  ,k-1,2,u.size)];
                    ut1 = u.m[id4(i  ,j  ,k+1,2,u.size)];
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
                    dwx = kx1 * 0.5 * (ue1 - uw1);
                    dwy = kx2 * 0.5 * (un1 - us1);
                    dwz = kx3 * 0.5 * (ut1 - ub1);

                    d01 = 2 * dux * dux;
                    d02 = 2 * dvy * dvy;
                    d03 = 2 * dwz * dwz;
                    d04 = (dwy + dvz) * (dwy + dvz);
                    d05 = (duz + dwx) * (duz + dwx);
                    d06 = (duy + dvx) * (duy + dvx);
                    D_u = sqrt(d01 + d02 + d03 + d04 + d05 + d06);
                    Del = cbrt(det);
                    l_0 = c.turbulence.cs * Del;

                    nut.m[id3(i,j,k,nut.size)] = l_0 * l_0 * D_u;
                }
            }
        }
    }
}

void Turbulence::csm(Dom &dom) {
    Ctrl                   &c   = dom.c;
    scalar_field<unsigned> &f   = dom.f;
    vector_field<real_t>   &u   = dom.u;
    vector_field<real_t>   &kx  = dom.kx;
    scalar_field<real_t>   &ja  = dom.ja;
    scalar_field<real_t>   &nut = dom.nut;

    #pragma acc kernels loop independent collapse(3) present(f, u, kx, ja, nut, c, dom)
    for (int i = GUIDE; i < dom.size[0] - GUIDE; i ++) {
        for (int j = GUIDE; j < dom.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < dom.size[2] - GUIDE; k ++) {
                if (Util::ibsee(f.m[id3(i,j,k,f.size)], Cell::Active, Util::Mask1)) {
                    unsigned b12, b13;
                    unsigned b22, b23;
                    unsigned b32, b33;
                    unsigned f12, f13;
                    unsigned f22, f23;
                    unsigned f32, f33;
                    unsigned m12, m13;
                    unsigned m22, m23;
                    unsigned m32, m33;
                    real_t   kx1, kx2, kx3;
                    real_t   uc0;
                    real_t   ue1, un1, ut1;
                    real_t   uw1, us1, ub1;
                    real_t   ref, dis;
                    real_t   dux, duy, duz;
                    real_t   dvx, dvy, dvz;
                    real_t   dwx, dwy, dwz;
                    real_t   det;
                    real_t   d01, d02, d03, d04, d05, d06;
                    real_t   Del, D_u;
                    real_t   q, e, fcs, C;

                    kx1 = kx.m[id4(i  ,j  ,k  ,0,kx.size)];
                    kx2 = kx.m[id4(i  ,j  ,k  ,1,kx.size)];
                    kx3 = kx.m[id4(i  ,j  ,k  ,2,kx.size)];
                    det = ja.m[id3(i  ,j  ,k  ,  ja.size)];
                    f12 =  f.m[id3(i-1,j  ,k  ,   f.size)];
                    f13 =  f.m[id3(i  ,j  ,k  ,   f.size)];
                    f22 =  f.m[id3(i  ,j-1,k  ,   f.size)];
                    f23 =  f.m[id3(i  ,j  ,k  ,   f.size)];
                    f32 =  f.m[id3(i  ,j  ,k-1,   f.size)];
                    f33 =  f.m[id3(i  ,j  ,k  ,   f.size)];
                    m12 = Util::ibsee(f12, Cell::Me, Util::Mask1);
                    m13 = Util::ibsee(f13, Cell::Me, Util::Mask1);
                    m22 = Util::ibsee(f22, Cell::Mn, Util::Mask1);
                    m23 = Util::ibsee(f23, Cell::Mn, Util::Mask1);
                    m32 = Util::ibsee(f32, Cell::Mt, Util::Mask1);
                    m33 = Util::ibsee(f33, Cell::Mt, Util::Mask1);
                    f12 = Util::ibsee(f12, Cell::Fe, Util::Mask8);
                    f13 = Util::ibsee(f13, Cell::Fe, Util::Mask8);
                    f22 = Util::ibsee(f22, Cell::Fn, Util::Mask8);
                    f23 = Util::ibsee(f23, Cell::Fn, Util::Mask8);
                    f32 = Util::ibsee(f32, Cell::Ft, Util::Mask8);
                    f33 = Util::ibsee(f33, Cell::Ft, Util::Mask8);
                    b12 = Util::ibsee(u.bflag[f12], 0, Util::Mask8);
                    b13 = Util::ibsee(u.bflag[f13], 0, Util::Mask8);
                    b22 = Util::ibsee(u.bflag[f22], 0, Util::Mask8);
                    b23 = Util::ibsee(u.bflag[f23], 0, Util::Mask8);
                    b32 = Util::ibsee(u.bflag[f32], 0, Util::Mask8);
                    b33 = Util::ibsee(u.bflag[f33], 0, Util::Mask8);

                    uc0 = u.m[id4(i  ,j  ,k  ,0,u.size)];
                    uw1 = u.m[id4(i-1,j  ,k  ,0,u.size)];
                    ue1 = u.m[id4(i+1,j  ,k  ,0,u.size)];
                    us1 = u.m[id4(i  ,j-1,k  ,0,u.size)];
                    un1 = u.m[id4(i  ,j+1,k  ,0,u.size)];
                    ub1 = u.m[id4(i  ,j  ,k-1,0,u.size)];
                    ut1 = u.m[id4(i  ,j  ,k+1,0,u.size)];
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
                    dux = kx1 * 0.5 * (ue1 - uw1);
                    duy = kx2 * 0.5 * (un1 - us1);
                    duz = kx3 * 0.5 * (ut1 - ub1);

                    uc0 = u.m[id4(i  ,j  ,k  ,1,u.size)];
                    uw1 = u.m[id4(i-1,j  ,k  ,1,u.size)];
                    ue1 = u.m[id4(i+1,j  ,k  ,1,u.size)];
                    us1 = u.m[id4(i  ,j-1,k  ,1,u.size)];
                    un1 = u.m[id4(i  ,j+1,k  ,1,u.size)];
                    ub1 = u.m[id4(i  ,j  ,k-1,1,u.size)];
                    ut1 = u.m[id4(i  ,j  ,k+1,1,u.size)];
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
                    dvx = kx1 * 0.5 * (ue1 - uw1);
                    dvy = kx2 * 0.5 * (un1 - us1);
                    dvz = kx3 * 0.5 * (ut1 - ub1);

                    uc0 = u.m[id4(i  ,j  ,k  ,2,u.size)];
                    uw1 = u.m[id4(i-1,j  ,k  ,2,u.size)];
                    ue1 = u.m[id4(i+1,j  ,k  ,2,u.size)];
                    us1 = u.m[id4(i  ,j-1,k  ,2,u.size)];
                    un1 = u.m[id4(i  ,j+1,k  ,2,u.size)];
                    ub1 = u.m[id4(i  ,j  ,k-1,2,u.size)];
                    ut1 = u.m[id4(i  ,j  ,k+1,2,u.size)];
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
                    dwx = kx1 * 0.5 * (ue1 - uw1);
                    dwy = kx2 * 0.5 * (un1 - us1);
                    dwz = kx3 * 0.5 * (ut1 - ub1);

                    d01 = 2 * dux * dux;
                    d02 = 2 * dvy * dvy;
                    d03 = 2 * dwz * dwz;
                    d04 = (dwy + dvz) * (dwy + dvz);
                    d05 = (duz + dwx) * (duz + dwx);
                    d06 = (duy + dvx) * (duy + dvx);
                    D_u = sqrt(d01 + d02 + d03 + d04 + d05 + d06);
                    Del = cbrt(det);

                    e   = dux * dux + duy * duy + duz * duz;
                    e  += dvx * dvx + dvy * dvy + dvz * dvz;
                    e  += dwx * dwx + dwy * dwy + dwz * dwz;
                    e  *= 0.5;
                    q   = dux * dux + duy * dvx + duz * dwx;
                    q  += dvx * duy + dvy * dvy + dvz * dwy;
                    q  += dwx * duz + dwy * dvz + dwz * dwz;
                    q  *= (- 0.5);
                    fcs = (q + copysign(1E-12, q)) / (e + copysign(1E-12, e));

                    real_t afcs = fabs(fcs);
                    C = sqrt(afcs * afcs * afcs) * (1 - fcs) / 22.0;

                    nut.m[id3(i,j,k,nut.size)] = C * Del * Del * D_u;
                }
            }
        }
    }
}

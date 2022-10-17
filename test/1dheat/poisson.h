#ifndef _POISSON_H
#define _POISSON_H 1

#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include "mesh.h"
#include "param.h"

static void poisson_prepare_eq(Mat<real_t> &a, Mat<real_t> &b, Mesh &mesh, Dom &dom, stencil_t s) {
    if (s == stencil_t::d1s3) {
        #pragma acc kernels loop independent collapse(3) present(a, b, mesh, dom)
        for (int i = 0; i < dom.size[0]; i ++) {
            for (int j = 0; j < dom.size[1]; j ++) {
                for (int k = 0; k < dom.size[2]; k ++) {
                    if (i == 0) {
                        int    cur = id(i,j,k,dom.size);
                        real_t dxr = mesh.x.get(id(i+1,j,k,dom.size), 0) - mesh.x.get(cur, 0);
                        real_t dx  = mesh.h.get(cur, 0);
                        real_t cr  = kappa / (dx * dxr);
                        real_t cl  = 0;
                        real_t cc  = - (cr + 2 * kappa / (dx * dx));
                        real_t rc  = - 2 * kappa * TL / (dx * dx);
                        a.get(cur, 0) = cl;
                        a.get(cur, 1) = cc;
                        a.get(cur, 2) = cr;
                        b.get(cur)    = rc;
                    } else if (i == dom.size[0] - 1) {
                        int    cur = id(i,j,k,dom.size);
                        real_t dxl = mesh.x.get(cur, 0) - mesh.x.get(id(i-1,j,k,dom.size), 0);
                        real_t dx  = mesh.h.get(cur, 0);
                        real_t cr  = 0;
                        real_t cl  = kappa / (dx * dxl);
                        real_t cc  = - (cl + 2 * kappa / (dx * dx));
                        real_t rc  = - 2 * kappa * TR / (dx * dx);
                        a.get(cur, 0) = cl;
                        a.get(cur, 1) = cc;
                        a.get(cur, 2) = cr;
                        b.get(cur)    = rc;
                    } else {
                        int    cur = id(i,j,k,dom.size);
                        real_t dxr = mesh.x.get(id(i+1,j,k,dom.size), 0) - mesh.x.get(cur, 0);
                        real_t dxl = mesh.x.get(cur, 0) - mesh.x.get(id(i-1,j,k,dom.size), 0);
                        real_t dx  = mesh.h.get(cur, 0);
                        real_t cr  = kappa / (dx * dxr);
                        real_t cl  = kappa / (dx * dxl);
                        real_t cc  = - (cr + cl);
                        real_t rc  = 0;
                        a.get(cur, 0) = cl;
                        a.get(cur, 1) = cc;
                        a.get(cur, 2) = cr;
                        b.get(cur)    = rc;
                    }
                }
            }
        }
    }
}

static void scale_eq(Mat<real_t> &a, Mat<real_t> &b, Dom &dom) {
    real_t max_diag = 0;
    #pragma acc kernels loop independent collapse(3) reduction(max:max_diag) present(a, dom) copy(max_diag)
    for (int i = 0; i < dom.size[0]; i ++) {
        for (int j = 0; j < dom.size[1]; j ++) {
            for (int k = 0; k < dom.size[2]; k ++) {
                int    cur = id(i,j,k,dom.size);
                real_t cc  = a.get(cur, 1);
                if (fabs(cc) > max_diag) {
                    max_diag = fabs(cc);
                }
            }
        }
    }
    #pragma acc kernels loop independent collapse(3) present(a, b, dom) copyin(max_diag)
    for (int i = 0; i < dom.size[0]; i ++) {
        for (int j = 0; j < dom.size[1]; j ++) {
            for (int k = 0; k < dom.size[2]; k ++) {
                int    cur = id(i,j,k,dom.size);
                for (int m = 0; m < a.col; m ++) {
                    a.get(cur, m) = a.get(cur, m) / max_diag;
                }
                b.get(cur) = b.get(cur) / max_diag;
            }
        }
    }
}

static void calc_res(Mat<real_t> &a, Mat<real_t> &x, Mat<real_t> &b, Mat<real_t> &r, real_t &norm, Mesh &mesh, Dom &dom, stencil_t s) {
    real_t sum = 0;
    if (s == stencil_t::d1s3) {
        /* #pragma acc kernels loop independent collapse(2) reduction(+:sum) present(a, x, b, r, dom) copy(sum)
        for (int i = 0; i < dom.num; i ++) {
            for (int m = 0; m < x.col; m ++) {
                int    i0 = (i - 1 >= 0      )? (i - 1) : i;
                int    i1 = i;
                int    i2 = (i + 1 <  dom.num)? (i + 1) : i;
                real_t _0 = x.get(i0, m);
                real_t _1 = x.get(i1, m);
                real_t _2 = x.get(i2, m);
                real_t c0 = a.get(i1, 0);
                real_t c1 = a.get(i1, 1);
                real_t c2 = a.get(i1, 2);
                real_t r1 = b.get(i1) - (c0 * _0 + c1 * _1 + c2 * _2);
                r.get(i1, m) = r1;
                sum += r1 * r1;
            }
        } */
        #pragma acc kernels loop independent collapse(3) reduction(+:sum) present(a, x, b, r, dom, mesh) copy(sum)
        for (int i = 0; i < dom.size[0]; i ++) {
            for (int j = 0; j < dom.size[1]; j ++) {
                for (int k = 0; k < dom.size[2]; k ++) {
                    int    cur = id(i,j,k,dom.size);
                    real_t cr  = a.get(cur, 2);
                    real_t cc  = a.get(cur, 1);
                    real_t cl  = a.get(cur, 0);
                    int    ir  = mesh.map.get(cur, 3);
                    int    il  = mesh.map.get(cur, 0);
                    int    ic  = cur;
                    for (int m = 0; m < x.col; m ++) {
                        real_t _r = x.get(ir, m);
                        real_t _l = x.get(il, m);
                        real_t _c = x.get(ic, m);
                        real_t rc = b.get(ic, m) - (_r * cr + _c * cc + _l * cl);
                        r.get(ic, m) = rc;
                        sum += rc * rc;
                    }
                }
            }
        }
    }
    norm = sqrt(sum / x.num);
}

static void poisson_sor(Mat<real_t> &a, Mat<real_t> &t, Mat<real_t> &rhs, Mat<real_t> &res, real_t &norm, Mesh &mesh, Dom &dom) {
    int it = 0;
    do {
        /* #pragma acc kernels loop independent present(a, t, rhs, res, mesh, dom)
        for (int i = 0; i <  dom.num; i += 2) {
            int    i0 = (i - 1 >= 0      )? (i - 1) : i;
            int    i1 = i;
            int    i2 = (i + 1 <  dom.num)? (i + 1) : i;
            real_t _0 = t.get(i0);
            real_t _1 = t.get(i1);
            real_t _2 = t.get(i2);
            real_t c0 = a.get(i1, 0);
            real_t c1 = a.get(i1, 1);
            real_t c2 = a.get(i1, 2);
            real_t _d = (rhs.get(i1) - c0 * _0 - c1 * _1 - c2 * _2) / c1;
            t.get(i1) = _1 + omega * _d;
        }
        #pragma acc kernels loop independent present(a, t, rhs, res, mesh, dom)
        for (int i = 1; i <  dom.num; i += 2) {
            int    i0 = (i - 1 >= 0      )? (i - 1) : i;
            int    i1 = i;
            int    i2 = (i + 1 <  dom.num)? (i + 1) : i;
            real_t _0 = t.get(i0);
            real_t _1 = t.get(i1);
            real_t _2 = t.get(i2);
            real_t c0 = a.get(i1, 0);
            real_t c1 = a.get(i1, 1);
            real_t c2 = a.get(i1, 2);
            real_t _d = (rhs.get(i1) - c0 * _0 - c1 * _1 - c2 * _2) / c1;
            t.get(i1) = _1 + omega * _d;
        } */
        #pragma acc kernels loop independent collapse(3) present(a, t, rhs, res, mesh, dom)
        for (int i = 0; i < dom.size[0]; i ++) {
            for (int j = 0; j < dom.size[1]; j ++) {
                for (int k = 0; k < dom.size[2]; k ++) {
                    if ((i + j + k) % 2 == 0) {
                        int    cur = id(i,j,k,dom.size);
                        real_t cr  = a.get(cur, 2);
                        real_t cc  = a.get(cur, 1);
                        real_t cl  = a.get(cur, 0);
                        int    ir  = mesh.map.get(cur, 3);
                        int    il  = mesh.map.get(cur, 0);
                        int    ic  = cur;
                        real_t _r  = t.get(ir);
                        real_t _l  = t.get(il);
                        real_t _c  = t.get(ic);
                        real_t _d  = (rhs.get(ic) - (_r * cr + _c * cc + _l * cl)) / cc;
                        t.get(ic)  = _c + omega * _d;
                    }

                }
            }
        }
        #pragma acc kernels loop independent collapse(3) present(a, t, rhs, res, mesh, dom)
        for (int i = 0; i < dom.size[0]; i ++) {
            for (int j = 0; j < dom.size[1]; j ++) {
                for (int k = 0; k < dom.size[2]; k ++) {
                    if ((i + j + k) % 2 == 1) {
                        int    cur = id(i,j,k,dom.size);
                        real_t cr  = a.get(cur, 2);
                        real_t cc  = a.get(cur, 1);
                        real_t cl  = a.get(cur, 0);
                        int    ir  = mesh.map.get(cur, 3);
                        int    il  = mesh.map.get(cur, 0);
                        int    ic  = cur;
                        real_t _r  = t.get(ir);
                        real_t _l  = t.get(il);
                        real_t _c  = t.get(ic);
                        real_t _d  = (rhs.get(ic) - (_r * cr + _c * cc + _l * cl)) / cc;
                        t.get(ic)  = _c + omega * _d;
                    }

                }
            }
        }
        calc_res(a, t, rhs, res, norm, mesh, dom, stencil_t::d1s3);
        
        if (it % 10000 == 0) {
            t.update_self();
            FILE *fo;
            char fname[128];
            sprintf(fname, "temperature.csv.%d", it / 10000);
            fo = fopen(fname, "w+t");
            if (fo) {
                fprintf(fo, "x,y,z,t\n");
                real_t x, y, z, tp;
                for (int k = 0; k < dom.size[2]; k ++) {
                    for (int j = 0; j < dom.size[1]; j ++) {
                        for (int i = 0; i < dom.size[0]; i ++) {
                            int cur = id(i,j,k,dom.size);
                            x  = mesh.x.get(cur, 0);
                            y  = mesh.x.get(cur, 1);
                            z  = mesh.x.get(cur, 2);
                            tp = t.get(cur);
                            fprintf(fo, "%10.3e,%10.3e,%10.3e,%10.3e\n", x, y, z, tp);
                        }
                    }
                }
            }
        }

        printf("\r%6d %10.3e", ++it, norm);
    } while (norm > 1e-7);
    printf("\n");
}

#endif
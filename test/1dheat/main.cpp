#include "poisson.h"

Dom dom(N, 1, 1, 0, 0, 0);

int main() {
    Mesh mesh(dom);
    real_t h = L / N;
    real_t x = 0;
    for (int i = 0; i < dom.size[0]; i ++) {
        for (int j = 0; j < dom.size[1]; j ++) {
            for (int k = 0; k < dom.size[2]; k ++) {
                int cur = id(i,j,k,dom.size);
                mesh.h.get(cur, 0) = h;
                mesh.h.get(cur, 1) = h;
                mesh.h.get(cur, 2) = h;
                mesh.x.get(cur, 0) = x + 0.5 * h;
                mesh.x.get(cur, 1) = 0;
                mesh.x.get(cur, 2) = 0;
                mesh.v.get(cur)    = h * h * h;
                x += h;
            }
        }
    }
    for (int i = 0; i < dom.size[0]; i ++) {
        for (int j = 0; j < dom.size[1]; j ++) {
            for (int k = 0; k < dom.size[2]; k ++) {
                int cur = id(i,j,k,dom.size);
                mesh.map.get(cur, 0) = (i > 0            )? id(i-1,j,k,dom.size) : cur;
                mesh.map.get(cur, 1) = (j > 0            )? id(i,j-1,k,dom.size) : cur;
                mesh.map.get(cur, 2) = (k > 0            )? id(i,j,k-1,dom.size) : cur;
                mesh.map.get(cur, 3) = (i < dom.size[0]-1)? id(i+1,j,k,dom.size) : cur;
                mesh.map.get(cur, 4) = (j < dom.size[1]-1)? id(i,j+1,k,dom.size) : cur;
                mesh.map.get(cur, 5) = (k < dom.size[2]-1)? id(i,j,k+1,dom.size) : cur;
            }
        }
    }
    Mat<real_t> t(dom, 1);
    Mat<real_t> a(dom, 3);
    Mat<real_t> rhs(dom, 1);
    Mat<real_t> res(dom, 1);
    real_t norm;

    mesh.to_device();
    dom.to_device();
    t.to_device();
    a.to_device();
    rhs.to_device();
    res.to_device();

    poisson_prepare_eq(a, rhs, mesh, dom, stencil_t::d1s3);
    scale_eq(a, rhs, dom);
    poisson_sor(a, t, rhs, res, norm, mesh, dom);

    t.update_self();
    a.update_self();
    rhs.update_self();
    mesh.off_device();
    dom.off_device();
    t.off_device();
    a.off_device();
    rhs.off_device();
    res.off_device();

    // for (int k = 0; k < dom.size[2]; k ++) {
    //     for (int j = 0; j < dom.size[1]; j ++) {
    //         for (int i = 0; i < dom.size[0]; i ++) {
    //             int cur = id(i,j,k,dom.size);
    //             printf("%10.3e %10.3e %10.3e\n", a.get(cur, 0), a.get(cur, 1), a.get(cur, 2));
    //         }
    //     }
    // }

    for(int i = 0; i < dom.num; i ++) {
        printf("%10.3e ", rhs.get(i));
    }
    printf("\n");

    for(int i = 0; i < dom.num; i ++) {
        printf("%10.3e ", t.get(i));
    }
    printf("\n");

    FILE *fo;
    fo = fopen("final.csv", "w+t");
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

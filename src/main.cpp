#include <stdio.h>
#include "../include/flag.h"
#include "../include/fluo.h"

void param_out(Dom &dom) {
    FILE *fo;
    char fname[128];
    sprintf(fname, "param.%s.csv", dom.name);
    fo = fopen(fname, "w+t");
    if (fo) {
        fprintf(fo, "x,y,z,active,fe,fn,ft,me,mn,mt,k1,k2,k3,g1,g2,g3,j\n");
        for (int k = 0; k < dom.size[2]; k ++) {
            for (int j = 0; j < dom.size[1]; j ++) {
                for (int i = 0; i < dom.size[0]; i ++) {
                    unsigned flag  = dom.f.m[id3(i,j,k,dom.f.size)];
                    unsigned active = Util::ibsee(flag, Flag::Active, Util::Mask1);
                    unsigned fe = Util::ibsee(flag, Flag::Fe, Util::Mask8);
                    unsigned fn = Util::ibsee(flag, Flag::Fn, Util::Mask8);
                    unsigned ft = Util::ibsee(flag, Flag::Ft, Util::Mask8);
                    unsigned me = Util::ibsee(flag, Flag::Me, Util::Mask1);
                    unsigned mn = Util::ibsee(flag, Flag::Mn, Util::Mask1);
                    unsigned mt = Util::ibsee(flag, Flag::Mt, Util::Mask1);
                    real_t   x1 =  dom.x.m[id4(i,j,k,0, dom.x.size)];
                    real_t   x2 =  dom.x.m[id4(i,j,k,1, dom.x.size)];
                    real_t   x3 =  dom.x.m[id4(i,j,k,2, dom.x.size)];
                    real_t   k1 = dom.kx.m[id4(i,j,k,0,dom.kx.size)];
                    real_t   k2 = dom.kx.m[id4(i,j,k,1,dom.kx.size)];
                    real_t   k3 = dom.kx.m[id4(i,j,k,2,dom.kx.size)];
                    real_t   g1 =  dom.g.m[id4(i,j,k,0, dom.g.size)];
                    real_t   g2 =  dom.g.m[id4(i,j,k,1, dom.g.size)];
                    real_t   g3 =  dom.g.m[id4(i,j,k,2, dom.g.size)];
                    real_t   de = dom.ja.m[id3(i,j,k,  dom.ja.size)];
                    fprintf(fo, "%.5e,%.5e,%.5e,%u,%u,%u,%u,%u,%u,%u,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e\n", x1, x2, x3, active, fe, fn, ft, me, mn, mt, k1, k2, k3, g1, g2, g3, de);
                }
            }
        }
        fclose(fo);
    }
}

int main() {
    FLUO fluo("./setting/setting.json");
    fluo.init();
    fluo.show_info();
    

    // param_out(domain);

    return 0;
}

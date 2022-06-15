#include <stdio.h>
#include <iostream>
#include <bitset>
#include "jreader.h"
#include "util.h"
#include "flag.h"

using namespace std;

int main(void) {
    JReader jr((char*)"./setting/setting.json");
    BC bc;
    Driver driver;
    jr.setup_bc(&bc, &driver);
    printf("%d boundaries\n", bc.n);
    for (int i = 1; i <= bc.n; i ++) {
        printf("b%d:\n", i);
        cout << "\t" << bitset<32>(bc.b[i].flag) << endl;
        printf("\t%.3lf, %.3lf, %.3lf, %.3lf, %.3lf\n", bc.b[i].u, bc.b[i].v, bc.b[i].w, bc.b[i].p, bc.b[i].nt);
    }
    printf("Outer BC:\n");
    cout << "\t" << bitset<32>(bc.f0) << endl;
    cout << "\t" << bitset<32>(bc.f1) << endl;
    cout << "\t" << bitset<32>(bc.f2) << endl;
    cout << "\t" << bitset<32>(bc.f3) << endl;
    cout << "\t" << bitset<32>(bc.f4) << endl;
    cout << "\t" << bitset<32>(bc.f5) << endl;

    printf("Driver:\n");
    printf("\t%d %d %d %.3lf\n", driver.on, driver.direction, driver.type, driver.value);

    D dom;
    jr.setup_domain(&dom);
    printf("Domain:\n");
    printf("\tdata domain: (%d %d %d)\n", dom.size[0], dom.size[1], dom.size[2]);
    printf("\tinenr domain: (%d %d %d)\n", dom.size[0] - 2 * D::GUIDE, dom.size[1] - 2 * D::GUIDE, dom.size[2] - 2 * D::GUIDE);
    FILE *fo;
    fo = fopen("para.csv", "w+t");
    if (fo) {
        fprintf(fo, "x,y,z,active,fe,fn,ft,me,mn,mt\n");
        for (int k = 0; k < dom.size[2]; k ++) {
            for (int j = 0; j < dom.size[1]; j ++) {
                for (int i = 0; i < dom.size[0]; i ++) {
                    unsigned int flag  = dom.F[id3(i,j,k,dom.size)];
                    unsigned int active = Util::ibsee(flag, Flag::Active, Util::Mask1);
                    unsigned int fe = Util::ibsee(flag, Flag::Fe, Util::Mask8);
                    unsigned int fn = Util::ibsee(flag, Flag::Fn, Util::Mask8);
                    unsigned int ft = Util::ibsee(flag, Flag::Ft, Util::Mask8);
                    unsigned int me = Util::ibsee(flag, Flag::Me, Util::Mask1);
                    unsigned int mn = Util::ibsee(flag, Flag::Mn, Util::Mask1);
                    unsigned int mt = Util::ibsee(flag, Flag::Mt, Util::Mask1);
                    fprintf(fo, "%d,%d,%d,%u,%u,%u,%u,%u,%u,%u\n", i, j, k, active, fe, fn, ft, me, mn, mt);
                }
            }
        }
        fclose(fo);
    }

    return 0;
}
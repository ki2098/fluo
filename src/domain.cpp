#include "../include/domain.h"

Dom::Dom(const char *label) {
    name = label;
}

void Dom::init(int *field, int nbound) {
    size[0] = field[0];
    size[1] = field[1];
    size[2] = field[2];
    len     = size[0] * size[1] * size[2];
    bcount  = nbound;
    f.init(size, bcount);
    u.init(size, 3, bcount);
    uu.init(size, 3, bcount);
    p.init(size, bcount);
    nut.init(size, bcount);
    div.init(size, bcount);
    x.init(size, 3, bcount);
    kx.init(size, 3, bcount);
    g.init(size, 3, bcount);
    ja.init(size, bcount);
}

void Dom::monitor() {

}

void Dom::driver() {
    for (int x = 0; x < 3; x ++) {
        if (c.driver[x].type == Ctrl::Driver::Type::off) {
            continue;
        }
        if (c.driver[x].var == Ctrl::Driver::Var::P) {
            c.driver[x].dp = c.driver[x].value;
        } else if (c.driver[x].var == Ctrl::Driver::Var::U) {
            real_t m1   = 0;
            real_t m2   = 0;
            real_t area = 0;
            if (c.driver[x].type == Ctrl::Driver::Type::directional) {
                #pragma acc kernels loop independent collapse(3) reduction(+:m1, m2, area) present(this[0:1], u, c, ja, kx) copy(m1, m2, area)
                for (int i = GUIDE; i < size[0] - GUIDE; i ++) {
                    for (int j = GUIDE; j < size[1] - GUIDE; j ++) {
                        for (int k = GUIDE; k < size[2] - GUIDE; k ++) {
                            real_t det = ja.m[id3(i,j,k,ja.size)];
                            m1 += u.m[id4(i,j,k,x,u.size)] * det;
                            m2 += c.driver[x].value        * det;
                            if (x == 0 && i == GUIDE || x == 1 && j == GUIDE || x == 2 && k == GUIDE) {
                                area += det * kx.m[id4(i,j,k,x,kx.size)];
                            }
                        }
                    }
                }
                c.driver[x].dp = (m2 - m1) / (area * c.time.dt);
                #pragma acc update device(c.driver[x].dp)
            } else if (c.driver[x].type == Ctrl::Driver::Type::driver) {


                // just leave it unfinished for the moment

            }
        }
    }
}

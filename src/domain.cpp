#include "../include/domain.h"
#include "../include/flag.h"

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
    sz.init(size, 3, bcount);
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
                            if ((x == 0 && i == GUIDE) || (x == 1 && j == GUIDE) || (x == 2 && k == GUIDE)) {
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

void Dom::driver_monitor() {
    for (int x = 0; x < 3; x ++) {
        if (c.driver[x].type == Ctrl::Driver::Type::off) {
            c.driver[x].u_observed = 0;
        } else if (c.driver[x].type == Ctrl::Driver::Type::directional) {
            real_t sum = 0;
            real_t cnt = 0;
            #pragma acc kernels loop independent collapse(3) reduction(+:sum, cnt) present(this[0:1], u, ja) copy(sum, cnt)
            for (int i = GUIDE; i < size[0] - GUIDE; i ++) {
                for (int j = GUIDE; j < size[1] - GUIDE; j ++) {
                    for (int k = GUIDE; k < size[2] - GUIDE; k ++) {
                        real_t det = ja.m[id3(i,j,k,ja.size)];
                        sum += u.m[id4(i,j,k,x,u.size)] * det;
                        cnt += det;
                    }
                }
            }
            c.driver[x].u_observed = sum / cnt;
        } else if (c.driver[x].type == Ctrl::Driver::Type::driver) {

            // just leave it unfinished for the moment

        }
    }
}

void Dom::pressure_zero_average() {
    real_t sum = 0;
    int    cnt = 0;

    #pragma acc kernels loop independent collapse(3) reduction(+:sum, cnt) present(this[0:1], f, p) copy(sum, cnt)
    for (int i = GUIDE; i < p.size[0] - GUIDE; i ++) {
        for (int j = GUIDE; j < p.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < p.size[2] - GUIDE; k ++) {
                if (Util::ibsee(f.m[id3(i,j,k,f.size)], Cell::Active, Util::Mask1)) {
                    sum += p.m[id3(i,j,k,p.size)];
                    cnt += 1;
                }
            }
        }
    }

    real_t avg = sum / cnt;

    #pragma acc kernels loop independent collapse(3) present(this[0:1], f, p) copyin(avg)
    for (int i = GUIDE; i < p.size[0] - GUIDE; i ++) {
        for (int j = GUIDE; j < p.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < p.size[2] - GUIDE; k ++) {
                if (Util::ibsee(f.m[id3(i,j,k,f.size)], Cell::Active, Util::Mask1)) {
                    p.m[id3(i,j,k,p.size)] -= avg;;
                }
            }
        }
    }
}

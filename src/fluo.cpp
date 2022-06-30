#include <stdio.h>
#include "../include/fluo.h"
#include "../include/boundary.h"
#include "../include/flag.h"

void FLUO::show_info() {
    Ctrl &c = domain.c;

    printf("Main domain:\n");
    printf("\t%d %d %d\n", domain.size[0] - 2 * GUIDE, domain.size[1] - 2 * GUIDE, domain.size[2] - 2 * GUIDE);
    printf("\t%d boundaries\n", domain.bcount - 1);
    for (int i = 1; i < domain.bcount; i ++) {
        printf("\tb%d:\n", i);
        printf("\t\tu   %u (%lf %lf %lf)", Util::ibsee(domain.u.bflag[i], 0, Util::Mask8), domain.u.b[id2(i,0,domain.u.bsize)], domain.u.b[id2(i,1,domain.u.bsize)], domain.u.b[id2(i,2,domain.u.bsize)]);
        if (Util::ibsee(domain.u.bflag[i], BB::uu_locked, Util::Mask1)) {
            printf(" uu locked");
        }
        printf("\n");
        printf("\t\tp   %u %lf\n", Util::ibsee(domain.p.bflag[i], 0, Util::Mask8), domain.p.b[i]);
        printf("\t\tnut %u %lf\n", Util::ibsee(domain.nut.bflag[i], 0, Util::Mask8), domain.nut.b[i]);
    }

    const char *outer_face[6] = {"x-", "x+", "y-", "y+", "z-", "z+"};
    printf("\tOuter BC:\n");
    for (int i = 0; i < 6; i ++) {
        printf("\t%s\tu %2u  p %2u  nut %2u", outer_face[i], domain.u.obflag[i], domain.p.obflag[i], domain.nut.obflag[i]);
        if (domain.uob.type[i] == Dom::UOB::Type::designated) {
            printf("  outflow %lf", domain.uob.value[i]);
        } else if (domain.uob.type[i] == Dom::UOB::Type::average) {
            printf("  outflow average");
        } else if (domain.uob.type[i] == Dom::UOB::Type::minmax) {
            printf("  outflow minmax");
        }
        printf("\n");
    }

    const char *directions[3] = {"x", "y", "z"};
    const char *driver_type[3] = {"off", "directional", "driver"};
    for (int i = 0; i < 3; i ++) {
        if (c.driver[i].type == Ctrl::Driver::Type::off) {
            printf("Driver-%s:\n\toff\n", directions[i]);
        } else {
            printf("Driver-%s:\n", directions[i]);
            printf("\t%s\n", driver_type[(int)c.driver[i].type]);
            if (c.driver[i].var == Ctrl::Driver::Var::U) {
                printf("\tu %lf\n", c.driver[i].value);
            } else if (c.driver[i].var == Ctrl::Driver::Var::P) {
                printf("\tp gradient %.2e\n", c.driver[i].value);
            }
            if (c.driver[i].type == Ctrl::Driver::Type::driver) {
                printf("\tinflow from %s\n", outer_face[c.driver[i].inflow]);
                printf("\tlength %d\n", c.driver[i].length);
            }
        }
    }

    printf("Poisson:\n");
    const char *lstype[4] = {"jacobi", "sor", "bicgstab", "pbicgstab"};
    printf("\t%s\n", lstype[(int)c.poisson.type]);
    if (c.poisson.type == Ctrl::LS::Type::sor) {
        printf("\tomega %lf\n", c.poisson.omega);
    }
    printf("\tresidual tolerance %.2e\n", c.poisson.epsilon);
    printf("\tmax iteration %d\n", c.poisson.maxit);

    printf("Flow:\n");
    printf("\tRe %lf\n", c.flow.re);
    printf("\tdiv tolerance %.2e\n", c.flow.tdiv);
    printf("\tmax iteration %d\n", c.flow.maxit);
    const char *scheme[3] = {"1st order upwind", "3rd order upwind", "muscl"};
    printf("\tconvective scheme %s\n", scheme[(int)c.flow.scheme]);

    if (c.monitor.type == Ctrl::Monitor::Type::off) {
        printf("Monitor:\n\toff\n");
    } else {
        printf("Monitor:\n");
        printf("\t%d\n", (int)c.monitor.type);
        if (c.monitor.type == Ctrl::Monitor::Type::T) {
            printf("\toutput every %.2e t\n", c.monitor.t_interval);
        } else if (c.monitor.type == Ctrl::Monitor::Type::Dt) {
            printf("\toutput every %d dt\n", c.monitor.dt_interval);
        }
    }
    
    printf("Time:\n");
    printf("\tdt %.2e", c.time.dt);
    printf("\ttarget t %lf\n", c.time.dt * c.time.ndt);
    printf("\ttarget number of dt %d\n", c.time.ndt);
}

void FLUO::param_out() {
    FILE *fo;
    char fname[128];
    sprintf(fname, "param.%s.csv", domain.name);
    fo = fopen(fname, "w+t");
    if (fo) {
        fprintf(fo, "x,y,z,active,fe,fn,ft,me,mn,mt,k1,k2,k3,g1,g2,g3,j,sz1,sz2,sz3\n");
        for (int k = 0; k < domain.size[2]; k ++) {
            for (int j = 0; j < domain.size[1]; j ++) {
                for (int i = 0; i < domain.size[0]; i ++) {
                    unsigned flag  = domain.f.m[id3(i,j,k,domain.f.size)];
                    unsigned active = Util::ibsee(flag, Flag::Active, Util::Mask1);
                    unsigned fe = Util::ibsee(flag, Flag::Fe, Util::Mask8);
                    unsigned fn = Util::ibsee(flag, Flag::Fn, Util::Mask8);
                    unsigned ft = Util::ibsee(flag, Flag::Ft, Util::Mask8);
                    unsigned me = Util::ibsee(flag, Flag::Me, Util::Mask1);
                    unsigned mn = Util::ibsee(flag, Flag::Mn, Util::Mask1);
                    unsigned mt = Util::ibsee(flag, Flag::Mt, Util::Mask1);
                    real_t   x1 =  domain.x.m[id4(i,j,k,0, domain.x.size)];
                    real_t   x2 =  domain.x.m[id4(i,j,k,1, domain.x.size)];
                    real_t   x3 =  domain.x.m[id4(i,j,k,2, domain.x.size)];
                    real_t   k1 = domain.kx.m[id4(i,j,k,0,domain.kx.size)];
                    real_t   k2 = domain.kx.m[id4(i,j,k,1,domain.kx.size)];
                    real_t   k3 = domain.kx.m[id4(i,j,k,2,domain.kx.size)];
                    real_t   g1 =  domain.g.m[id4(i,j,k,0, domain.g.size)];
                    real_t   g2 =  domain.g.m[id4(i,j,k,1, domain.g.size)];
                    real_t   g3 =  domain.g.m[id4(i,j,k,2, domain.g.size)];
                    real_t   de = domain.ja.m[id3(i,j,k,  domain.ja.size)];
                    real_t  sz1 = domain.sz.m[id4(i,j,k,0,domain.sz.size)];
                    real_t  sz2 = domain.sz.m[id4(i,j,k,1,domain.sz.size)];
                    real_t  sz3 = domain.sz.m[id4(i,j,k,2,domain.sz.size)];
                    fprintf(fo, "%.5e,%.5e,%.5e,%u,%u,%u,%u,%u,%u,%u,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e,%.5e\n", x1, x2, x3, active, fe, fn, ft, me, mn, mt, k1, k2, k3, g1, g2, g3, de, sz1, sz2, sz3);
                }
            }
        }
        fclose(fo);
    }
}

void FLUO::var_out(const char* fname) {
    domain.x.update_self();
    domain.u.update_self();
    domain.uu.update_self();
    domain.nut.update_self();
    domain.div.update_self();
    mmac.diva.update_self();
    domain.p.update_self();
    FILE *fo;
    fo = fopen(fname, "w+t");
    if (fo == NULL) {
        printf("\nERROR when opening file\n");
        fflush(stdout);
    }
    else {
        fprintf(fo, "x,y,z,u,v,w,uu,vv,ww,p,nue,dvr,dva\n");
        for (int k = 0; k < domain.size[2]; k ++) {
            for (int j = 0; j < domain.size[1]; j ++) {
                for (int i = 0; i < domain.size[0]; i ++) {
                    real_t   x1 =   domain.x.m[id4(i,j,k,0,  domain.x.size)];
                    real_t   x2 =   domain.x.m[id4(i,j,k,1,  domain.x.size)];
                    real_t   x3 =   domain.x.m[id4(i,j,k,2,  domain.x.size)];
                    real_t   u1 =   domain.u.m[id4(i,j,k,0,  domain.u.size)];
                    real_t   u2 =   domain.u.m[id4(i,j,k,1,  domain.u.size)];
                    real_t   u3 =   domain.u.m[id4(i,j,k,2,  domain.u.size)];
                    real_t  uu1 =  domain.uu.m[id4(i,j,k,0, domain.uu.size)];
                    real_t  uu2 =  domain.uu.m[id4(i,j,k,1, domain.uu.size)];
                    real_t  uu3 =  domain.uu.m[id4(i,j,k,2, domain.uu.size)];
                    real_t  nut = domain.nut.m[id3(i,j,k,  domain.nut.size)];
                    real_t  div = domain.div.m[id3(i,j,k,  domain.div.size)];
                    real_t diva =  mmac.diva.m[id3(i,j,k,   mmac.diva.size)];
                    real_t    p =   domain.p.m[id3(i,j,k,    domain.p.size)];
                    fprintf(fo, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", x1, x2, x3, u1, u2, u3, uu1, uu2, uu3, p, nut, div, diva);
                }
            }
        }
        fclose(fo);
    }
}

void FLUO::fractional_step() {
    Ctrl &c = domain.c;

    mmac.ua.init(domain.u);
    mmac.uc.init(domain.u);
    mmac.uua.init(domain.uu);
    mmac.diva.init(domain.div);
    if (c.poisson.type == Ctrl::LS::Type::jacobi) {
        poisson.pd.init(domain.p);
    }

    domain.to_device();
    c.to_device();
    poisson.to_device();
    mmac.to_device();

    domain.f.to_device();
    domain.u.to_device();
    domain.uu.to_device();
    domain.p.to_device();
    domain.nut.to_device();
    domain.div.to_device();
    domain.x.to_device();
    domain.kx.to_device();
    domain.g.to_device();
    domain.ja.to_device();
    domain.sz.to_device();

    mmac.ua.to_device();
    mmac.uc.to_device();
    mmac.uua.to_device();
    mmac.diva.to_device();
    if (c.poisson.type == Ctrl::LS::Type::jacobi) {
        poisson.pd.to_device();
    }

    BB::vector_outer(domain.u, domain);
    BB::scalar_outer(domain.p, domain);
    domain.driver();
    mmac.interpolate_velocity(domain.u, mmac.uc, domain.uu, domain);

    real_t diva;
    int n_file = 0;

    char fname[128];
    sprintf(fname, "./data/var.csv.%d", n_file++);
    var_out(fname);

    for (c.time.idt = 1; c.time.idt <= c.time.ndt; domain.tick()) {
        mmac.pseudo_velocity(domain.u, domain.uu, mmac.ua, domain.nut, domain);
        BB::vector_outer(mmac.ua, domain);
        BB::vector_outflow(domain.u, domain);
        BB::velocity_outflow_correction(domain);
        mmac.interpolate_velocity(mmac.ua, mmac.uc, mmac.uua, domain);
        mmac.divergence_velocity(mmac.uua, mmac.diva, diva, domain);
        domain.driver();

        c.flow.it = 0;
        do {
            if (c.poisson.type == Ctrl::LS::Type::jacobi) {
                poisson.jacobi(domain.p, mmac.diva, domain);
            } else if (c.poisson.type == Ctrl::LS::Type::sor) {
                poisson.sor(domain.p, mmac.diva, domain);
            }
            mmac.correct_center_velocity(domain.u, mmac.ua, domain.p, domain);
            BB::vector_outer(domain.u, domain);
            mmac.correct_face_velocity(domain.u, domain.uu, mmac.uua, domain.p, domain);
            mmac.divergence_velocity(domain.uu, domain.div, c.flow.div, domain);

            printf("\r(%6d %6d),p(%4d %7.2e),d(%7.2e %7.2e)", c.time.idt, c.flow.it+1, c.poisson.it, c.poisson.res, diva, c.flow.div);
            fflush(stdout);
        } while (++c.flow.it < c.flow.maxit && c.flow.div > c.flow.tdiv);

        if (c.flow.it >= c.flow.maxit && c.flow.div > c.flow.tdiv) {
            goto Terminate;
        }

        domain.pressure_zero_average();
        BB::scalar_outer(domain.p, domain);

        if (c.time.idt % c.monitor.dt_interval == 0) {
            sprintf(fname, "./data/var.csv.%d", n_file++);
            var_out(fname);
        }
    }
Terminate:
    printf("\n");
    var_out("./data/final.csv");

    domain.f.end_device();
    domain.u.end_device();
    domain.uu.end_device();
    domain.p.end_device();
    domain.nut.end_device();
    domain.div.end_device();
    domain.x.end_device();
    domain.kx.end_device();
    domain.g.end_device();
    domain.ja.end_device();
    domain.sz.end_device();

    mmac.ua.end_device();
    mmac.uc.end_device();
    mmac.uua.end_device();
    mmac.diva.end_device();
    if (c.poisson.type == Ctrl::LS::Type::jacobi) {
        poisson.pd.end_device();
    }

    
    c.end_device();
    poisson.end_device();
    mmac.end_device();
    domain.end_device();
}

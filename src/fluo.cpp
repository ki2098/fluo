#include <stdio.h>
#include "../include/fluo.h"
#include "../include/boundary.h"

void FLUO::show_info() {
    Ctrl &c = domain.c;

    printf("Main domain:\n");
    printf("\t%d %d %d\n", domain.size[0] - 2 * GUIDE, domain.size[1] - 2 * GUIDE, domain.size[2] - 2 * GUIDE);
    printf("\t%d boundaries\n", domain.bcount - 1);
    for (int i = 1; i < domain.bcount; i ++) {
        printf("\tb%d:\n", i);
        printf("\t\tu   %u (%lf %lf %lf)\n", domain.u.bflag[i], domain.u.b[id2(i,0,domain.u.bsize)], domain.u.b[id2(i,1,domain.u.bsize)], domain.u.b[id2(i,2,domain.u.bsize)]);
        printf("\t\tp   %u %lf\n", domain.p.bflag[i], domain.p.b[i]);
        printf("\t\tnut %u %lf\n", domain.nut.bflag[i], domain.nut.b[i]);
    }

    const char *outer_face[6] = {"x-", "x+", "y-", "y+", "z-", "z+"};
    printf("\tOuter BC:\n");
    for (int i = 0; i < 6; i ++) {
        printf("\t%s\tu %2u  p %2u  nut %2u\n", outer_face[i], domain.u.obflag[i], domain.p.obflag[i], domain.nut.obflag[i]);
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

void FLUO::fractional_step() {
    Ctrl &c = domain.c;

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

    mmac.ua.init(domain.u);
    mmac.ua.to_device();
    mmac.uc.init(domain.u);
    mmac.uc.to_device();
    mmac.uua.init(domain.uu);
    mmac.uua.to_device();
    mmac.diva.init(domain.div);
    mmac.diva.to_device();

    if (c.poisson.type == Ctrl::LS::Type::jacobi) {
        poisson.pd.init(domain.p);
        poisson.pd.to_device();
    }

    BB::vector_outer(domain.u, domain);
    BB::scalar_outer(domain.p, domain);
    domain.driver();
    mmac.interpolate_velocity(domain.u, mmac.uc, domain.uu, domain);

    real_t diva;

    for (c.time.idt = 1; c.time.idt <= c.time.ndt; domain.tick()) {
        mmac.pseudo_velocity(domain.u, domain.uu, mmac.ua, domain.nut, domain);
        BB::vector_outer(mmac.ua, domain);
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

            printf("\r(%6d %6d),p(%4d %7.2e),d(%7.2e %7.2e)", c.time.idt, c.flow.it, c.poisson.it, c.poisson.res, diva, c.flow.div);
        } while (c.flow.div > c.flow.tdiv && ++c.flow.it < c.flow.maxit);
    }

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

    mmac.ua.end_device();
    mmac.uc.end_device();
    mmac.uua.end_device();
    mmac.diva.end_device();

    if (c.poisson.type == Ctrl::LS::Type::jacobi) {
        poisson.pd.end_device();
    }
}

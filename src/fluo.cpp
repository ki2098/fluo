#include <stdio.h>
#include "../include/fluo.h"
#include "../include/boundary.h"
#include "../include/flag.h"

void FLUO::time_sum() {
    Dom                  &dom   = domain;
    scalar_field<real_t> &p_avg = dom.p_avg;
    scalar_field<real_t> &p     = dom.p;
    vector_field<real_t> &u_avg = dom.u_avg;
    vector_field<real_t> &u     = dom.u;
    Ctrl                 &c     = dom.c;
    #pragma acc kernels loop independent collapse(3) present(p_avg, p, u_avg, u, dom)
    for (int i = 0; i < dom.size[0]; i ++) {
        for (int j = 0; j < dom.size[1]; j ++) {
            for (int k = 0; k < dom.size[2]; k ++) {
                p_avg.m[id3(i,j,k,  p_avg.size)] += p.m[id3(i,j,k,  p.size)];
                u_avg.m[id4(i,j,k,0,u_avg.size)] += u.m[id4(i,j,k,0,u.size)];
                u_avg.m[id4(i,j,k,1,u_avg.size)] += u.m[id4(i,j,k,1,u.size)];
                u_avg.m[id4(i,j,k,2,u_avg.size)] += u.m[id4(i,j,k,2,u.size)];
            }
        }
    }
    c.statistics.avg_steps += 1;
}

void FLUO::time_average() {
    Dom                  &dom   = domain;
    scalar_field<real_t> &p_avg = dom.p_avg;
    vector_field<real_t> &u_avg = dom.u_avg;
    Ctrl                 &c     = dom.c;
    #pragma acc update device(c.statistics.avg_steps)
    #pragma acc kernels loop independent collapse(3) present(p_avg, u_avg, c, dom)
    for (int i = 0; i < dom.size[0]; i ++) {
        for (int j = 0; j < dom.size[1]; j ++) {
            for (int k = 0; k < dom.size[2]; k ++) {
                p_avg.m[id3(i,j,k,  p_avg.size)] /= c.statistics.avg_steps;
                u_avg.m[id4(i,j,k,0,u_avg.size)] /= c.statistics.avg_steps;
                u_avg.m[id4(i,j,k,1,u_avg.size)] /= c.statistics.avg_steps;
                u_avg.m[id4(i,j,k,2,u_avg.size)] /= c.statistics.avg_steps;
            }
        }
    }
}

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
        if (Util::ibsee(domain.u.bflag[i], BB::wall_func, Util::Mask1)) {
            printf(" wall function");
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
    
    printf("Turbulence:\n");
    if (c.turbulence.model == Ctrl::Turbulence::Model::off) {
        printf("\toff\n");
    } else if (c.turbulence.model == Ctrl::Turbulence::Model::smagorinsky) {
        printf("\tsmagorinsky %.2lf\n", c.turbulence.cs);
    } else if (c.turbulence.model == Ctrl::Turbulence::Model::csm) {
        printf("\tcsm\n");
    }

    printf("Time:\n");
    printf("\tdt %.2e", c.time.dt);
    printf("\ttarget t %lf\n", c.time.dt * c.time.ndt);
    printf("\ttarget number of dt %d\n", c.time.ndt);

    printf("Time average:\n");
    if (c.statistics.type == Ctrl::Statistics::Type::off) {
        printf("\toff\n");
    } else {
        printf("\taverage from %.1lf\n", c.statistics.avg_from);
    }
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
                    unsigned active = Util::ibsee(flag, Cell::Active, Util::Mask1);
                    unsigned fe = Util::ibsee(flag, Cell::Fe, Util::Mask8);
                    unsigned fn = Util::ibsee(flag, Cell::Fn, Util::Mask8);
                    unsigned ft = Util::ibsee(flag, Cell::Ft, Util::Mask8);
                    unsigned me = Util::ibsee(flag, Cell::Me, Util::Mask1);
                    unsigned mn = Util::ibsee(flag, Cell::Mn, Util::Mask1);
                    unsigned mt = Util::ibsee(flag, Cell::Mt, Util::Mask1);
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
    domain.nut.update_self();
    domain.div.update_self();
    mmac.diva.update_self();
    domain.p.update_self();
    FILE *fo;
    fo = fopen(fname, "w+t");
    if (fo == NULL) {
        printf("\nERROR when opening file\n");
        fflush(stdout);
    } else {
        fprintf(fo, "x,y,z,u,v,w,p,nue,dvr,dva\n");
        for (int k = 0; k < domain.size[2]; k ++) {
            for (int j = 0; j < domain.size[1]; j ++) {
                for (int i = 0; i < domain.size[0]; i ++) {
                    real_t   x1 =   domain.x.m[id4(i,j,k,0,  domain.x.size)];
                    real_t   x2 =   domain.x.m[id4(i,j,k,1,  domain.x.size)];
                    real_t   x3 =   domain.x.m[id4(i,j,k,2,  domain.x.size)];
                    real_t   u1 =   domain.u.m[id4(i,j,k,0,  domain.u.size)];
                    real_t   u2 =   domain.u.m[id4(i,j,k,1,  domain.u.size)];
                    real_t   u3 =   domain.u.m[id4(i,j,k,2,  domain.u.size)];
                    real_t  nut = domain.nut.m[id3(i,j,k,  domain.nut.size)];
                    real_t  div = domain.div.m[id3(i,j,k,  domain.div.size)];
                    real_t diva =  mmac.diva.m[id3(i,j,k,   mmac.diva.size)];
                    real_t    p =   domain.p.m[id3(i,j,k,    domain.p.size)];
                    fprintf(fo, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%.5e,%.5e,%lf\n", x1, x2, x3, u1, u2, u3, p, nut, div, diva);
                }
            }
        }
        fclose(fo);
    }
}

void FLUO::time_avg_out(const char *fname) {
    domain.x.update_self();
    domain.u_avg.update_self();
    domain.p_avg.update_self();
    FILE *fo;
    fo = fopen(fname, "w+t");
    if (fo == NULL) {
        printf("\nERROR when opening file\n");
        fflush(stdout);
    } else {
        fprintf(fo, "x,y,z,u,v,w,p\n");
        for (int k = 0; k < domain.size[2]; k ++) {
            for (int j = 0; j < domain.size[1]; j ++) {
                for (int i = 0; i < domain.size[0]; i ++) {
                    real_t   x1 =   domain.x.m[id4(i,j,k,0,  domain.x.size)];
                    real_t   x2 =   domain.x.m[id4(i,j,k,1,  domain.x.size)];
                    real_t   x3 =   domain.x.m[id4(i,j,k,2,  domain.x.size)];
                    real_t   u1 =   domain.u_avg.m[id4(i,j,k,0,  domain.u.size)];
                    real_t   u2 =   domain.u_avg.m[id4(i,j,k,1,  domain.u.size)];
                    real_t   u3 =   domain.u_avg.m[id4(i,j,k,2,  domain.u.size)];
                    real_t    p =   domain.p_avg.m[id3(i,j,k,    domain.p.size)];
                    fprintf(fo, "%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", x1, x2, x3, u1, u2, u3, p);
                }
            }
        }
        fclose(fo);
    }
}

void FLUO::fractional_step() {
    Ctrl          &c       = domain.c;
    MMAC::Poisson &poisson = mmac.poisson;

    mmac.ua.init(domain.u);
    mmac.uc.init(domain.u);
    mmac.uua.init(domain.uu);
    mmac.diva.init(domain.div);
    mmac.psi.init(domain.div);
    if (c.poisson.type == Ctrl::LS::Type::jacobi) {
        poisson.pd.init(domain.p);
    }
    if (c.poisson.subtype == Ctrl::LS::Type::jacobi) {
        poisson.pd.init(domain.p);
        BB::disable_driver(poisson.pd);
    }
    if (c.poisson.type == Ctrl::LS::Type::pbicgstab) {
        poisson.pcg_r.init(domain.p);
        poisson.pcg_r0.init(domain.p);
        poisson.pcg_p.init(domain.p);
        poisson.pcg_p_.init(domain.p);
        poisson.pcg_q.init(domain.p);
        poisson.pcg_s.init(domain.p);
        poisson.pcg_s_.init(domain.p);
        poisson.pcg_t.init(domain.p);
        poisson.pcg_t_.init(domain.p);
        BB::disable_driver(poisson.pcg_r);
        BB::disable_driver(poisson.pcg_r0);
        BB::disable_driver(poisson.pcg_p);
        BB::disable_driver(poisson.pcg_p_);
        BB::disable_driver(poisson.pcg_q);
        BB::disable_driver(poisson.pcg_s);
        BB::disable_driver(poisson.pcg_s_);
        BB::disable_driver(poisson.pcg_t);
        BB::disable_driver(poisson.pcg_t_);
    }
    if (c.statistics.type == Ctrl::Statistics::Type::on) {
        domain.u_avg.init(domain.u);
        domain.p_avg.init(domain.p);
        c.statistics.avg_steps = 0;
    }

    domain.to_device();
    c.to_device();
    mmac.to_device();
    poisson.to_device();

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
    mmac.psi.to_device();
    if (c.poisson.type == Ctrl::LS::Type::jacobi || c.poisson.subtype == Ctrl::LS::Type::jacobi) {
        poisson.pd.to_device();
    }
    if (c.poisson.type == Ctrl::LS::Type::pbicgstab) {
        poisson.pcg_r.to_device();
        poisson.pcg_r0.to_device();
        poisson.pcg_p.to_device();
        poisson.pcg_p_.to_device();
        poisson.pcg_q.to_device();
        poisson.pcg_s.to_device();
        poisson.pcg_s_.to_device();
        poisson.pcg_t.to_device();
        poisson.pcg_t_.to_device();
    }
    if (c.statistics.type == Ctrl::Statistics::Type::on) {
        domain.u_avg.to_device();
        domain.p_avg.to_device();
    }

    BB::vector_outer(domain.u, domain);
    domain.driver();
    BB::scalar_outer(domain.p, domain);
    mmac.interpolate_velocity(domain.u, mmac.uc, domain.uu, domain);
    if (c.turbulence.model == Ctrl::Turbulence::Model::csm) {
        turb.csm(domain);
    } else if (c.turbulence.model == Ctrl::Turbulence::Model::smagorinsky) {
        turb.smagorinsky(domain);
    }
    BB::scalar_outer(domain.nut, domain);

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
        mmac.calc_rhs(mmac.diva, domain);
        domain.driver();

        c.flow.it = 0;
        do {
            if (c.poisson.type == Ctrl::LS::Type::jacobi) {
                poisson.jacobi(domain.p, mmac.psi, domain);
            } else if (c.poisson.type == Ctrl::LS::Type::sor) {
                poisson.sor(domain.p, mmac.psi, domain);
            } else if (c.poisson.type == Ctrl::LS::Type::pbicgstab) {
                poisson.pbicgstab(domain.p, mmac.psi, domain);
            }
            mmac.correct_center_velocity(domain.u, mmac.ua, domain.p, domain);
            BB::vector_outer(domain.u, domain);
            mmac.correct_face_velocity(domain.u, domain.uu, mmac.uua, domain.p, domain);
            mmac.divergence_velocity(domain.uu, domain.div, c.flow.div, domain);

            domain.driver_monitor();

            printf("\r(%6d %4d)p(%4d %7.2e)d(%7.2e %7.2e)(%6.3lf %6.3lf %6.3lf)", c.time.idt, c.flow.it+1, c.poisson.it, c.poisson.err, diva, c.flow.div, c.driver[0].u_observed, c.driver[1].u_observed, c.driver[2].u_observed);
            fflush(stdout);
        } while (++c.flow.it < c.flow.maxit && c.flow.div > c.flow.tdiv);

        if (c.turbulence.model == Ctrl::Turbulence::Model::csm) {
            turb.csm(domain);
        } else if (c.turbulence.model == Ctrl::Turbulence::Model::smagorinsky) {
            turb.smagorinsky(domain);
        }
        BB::scalar_outflow(domain.nut, domain);
        BB::scalar_outer(domain.nut, domain);

        domain.pressure_zero_average();
        BB::scalar_outer(domain.p, domain);

        if (c.time.idt % c.monitor.dt_interval == 0) {
            sprintf(fname, "./data/var.csv.%d", n_file++);
            var_out(fname);
        }

        if (c.statistics.type == Ctrl::Statistics::Type::on && c.time.idt * c.time.dt > c.statistics.avg_from) {
            time_sum();
        }

        if (c.flow.it >= c.flow.maxit && c.flow.div > c.flow.tdiv) {
            goto Terminate;
        }
    }
Terminate:
    printf("\n");
    var_out("./data/final.csv");

    if (c.statistics.type == Ctrl::Statistics::Type::on) {
        time_average();
        printf("time average of %.1lf t and %d steps", c.time.idt * c.time.dt - c.statistics.avg_from, c.statistics.avg_steps);
        time_avg_out("./data/time_average.csv");
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
    domain.sz.end_device();

    mmac.ua.end_device();
    mmac.uc.end_device();
    mmac.uua.end_device();
    mmac.diva.end_device();
    mmac.psi.end_device();
    if (c.poisson.type == Ctrl::LS::Type::jacobi || c.poisson.subtype == Ctrl::LS::Type::jacobi) {
        poisson.pd.end_device();
    }
    if (c.poisson.type == Ctrl::LS::Type::pbicgstab) {
        poisson.pcg_r.end_device();
        poisson.pcg_r0.end_device();
        poisson.pcg_p.end_device();
        poisson.pcg_p_.end_device();
        poisson.pcg_q.end_device();
        poisson.pcg_s.end_device();
        poisson.pcg_s_.end_device();
        poisson.pcg_t.end_device();
        poisson.pcg_t_.end_device();
    }
    if (c.statistics.type == Ctrl::Statistics::Type::on) {
        domain.u_avg.end_device();
        domain.p_avg.end_device();
    }
    
    c.end_device();
    poisson.end_device();
    mmac.end_device();
    domain.end_device();
}

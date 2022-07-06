#ifndef _CONTROL_H
#define _CONTROL_H 1

#include "basic.h"

class Ctrl {
public:
    struct Time {
    public:
        real_t dt;
        int    ndt;
        int    idt;
    };

    struct LS {
    public:
        enum class Type {jacobi, sor, bicgstab, pbicgstab};
        Ctrl::LS::Type type;
        Ctrl::LS::Type subtype;
        real_t omega;
        real_t epsilon;
        real_t err;
        int    maxit;
        int    it;
    };

    struct Flow {
    public:
        enum class Scheme {upwind1, upwind3, muscl};
        Ctrl::Flow::Scheme scheme;
        real_t alpha;
        real_t re;
        real_t ri;
        real_t tdiv;
        real_t div;
        int    maxit;
        int    it;
    };

    struct Monitor {
    public:
        enum class Type {off, T, Dt};
        Ctrl::Monitor::Type type;
        real_t t_interval;
        int dt_interval;
    };

    struct Driver {
    public:
        enum class Type {off, directional, driver};
        enum class Var {U, P};
        enum class Face {minus, plus};
        Ctrl::Driver::Type type;
        Ctrl::Driver::Var var;
        real_t value;
        real_t dp;
        real_t u_observed;
        int inflow;
        int length;
    };

    struct Turbulence {
    public:
        enum class Model {off, smagorinsky, csm};
        Ctrl::Turbulence::Model model;
        real_t cs;
    };

    struct Statistics {
    public:
        enum class Type {off, on};
        Ctrl::Statistics::Type type;
        real_t avg_from;
        int avg_steps;
    };

public:
    Ctrl::Time       time;
    Ctrl::LS         poisson;
    Ctrl::Flow       flow;
    Ctrl::Monitor    monitor;
    Ctrl::Driver     driver[3];
    Ctrl::Turbulence turbulence;
    Ctrl::Statistics statistics;
public:
    Ctrl() {
        monitor.type     = Ctrl::Monitor::Type::off;
        driver[0].type   = Ctrl::Driver::Type::off;
        driver[1].type   = Ctrl::Driver::Type::off;
        driver[2].type   = Ctrl::Driver::Type::off;
        turbulence.model = Ctrl::Turbulence::Model::off;
        statistics.type  = Ctrl::Statistics::Type::off;
    }
    void to_device() {
        #pragma acc enter data copyin(this[0:1])
    }
    void end_device() {
        #pragma acc exit data delete(this[0:1])
    }
    void update_self() {
        #pragma acc update self(this[0:1])
    }
    void update_device() {
        #pragma acc update device(this[0:1])
    }
};

#endif

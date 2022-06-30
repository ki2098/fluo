#ifndef _FLUO_H
#define _FLUO_H 1

#include "jreader.h"
#include "domain.h"
#include "poisson.h"
#include "mmac.h"

class FLUO {
public:
    JReader jreader;
    Dom     domain;
    Poisson poisson;
    MMAC    mmac;
public:
    FLUO(const char *fname) : jreader(fname), domain("domain") {}

    void init(bool load_mesh = true) {
        jreader.load_domain(domain, load_mesh);
    }

    void fractional_step();
    void show_info();
    void param_out();
    void var_out(const char* fname);
};

#endif

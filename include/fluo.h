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

    void init() {
        jreader.load_domain(domain);
    }

    void fractional_step();
    void show_info();
};

#endif

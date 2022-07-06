#ifndef _FLUO_H
#define _FLUO_H 1

#include "jreader.h"
#include "domain.h"
#include "mmac.h"
#include "turbulence.h"

class FLUO {
public:
    JReader              jreader;
    Dom                  domain;
    MMAC                 mmac;
    Turbulence           turb;
public:
    FLUO(const char *fname) : jreader(fname), domain("domain") {}

    void init(bool load_mesh = true) {
        jreader.load_domain(domain, load_mesh);
    }

    void fractional_step();
    void show_info();
    void param_out();
    void var_out(const char* fname);
    void time_sum();
    void time_average();
    void time_avg_out(const char *fname);
};

#endif

#ifndef _BOUNDARY_H
#define _BOUNDARY_H 1

#include "scalarfield.h"
#include "vectorfield.h"
#include "domain.h"

class BB {
public:

static const unsigned symmetric   = 0U;
static const unsigned cyclic      = 1U;
static const unsigned directional = 2U;
static const unsigned semicyclic  = 3U;
static const unsigned driver      = 4U;
static const unsigned outflow     = 5U;

static const unsigned dirichlet   = 0U;
static const unsigned neumann     = 1U;

static void pre(unsigned m, real_t _0, real_t _1, real_t &ref, real_t &dis) {
    ref = (m)? _1 : _0;
    dis = 0.5 - m;
}

static real_t eva(unsigned flag, real_t ref, real_t dis, real_t value) {
    if (Util::ibsee(flag, BB::dirichlet, Util::Mask1)) {
        return value;
    } else if (Util::ibsee(flag, BB::neumann, Util::Mask1)) {
        return ref + dis * value;
    }
    return 0.0;
}

static void scalar_outer(scalar_field<real_t> &fld, Dom &dom);
static void vector_outer(vector_field<real_t> &fld, Dom &dom);

};
#endif

#ifndef _BOUNDARY_H
#define _BOUNDARY_H 1

#include "scalarfield.h"
#include "vectorfield.h"
#include "domain.h"

class BB {
public:
    static const unsigned symmetric   = 0U;
    static const unsigned slip        = 1U;
    static const unsigned cyclic      = 2U;
    static const unsigned directional = 3U;
    static const unsigned semicyclic  = 4U;
    static const unsigned driver      = 5U;
    static const unsigned outflow     = 6U;

    static const unsigned dirichlet   = 0U;
    static const unsigned neumann     = 1U;
    static const unsigned uu_locked   = 8U;
    static const unsigned wall_func   = 9U;
public:
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
    static void disable_driver(scalar_field<real_t> &field) {
        for (int i = 0; i < 6; i ++) {
            if (Util::ibsee(field.obflag[i], BB::directional, Util::Mask1)) {
                field.obflag[i] = Util::ibset(field.obflag[i], BB::directional, Util::Mask1, 0U);
                field.obflag[i] = Util::ibset(field.obflag[i], BB::cyclic, Util::Mask1, 1U);
            }
            if (Util::ibsee(field.obflag[i], BB::driver, Util::Mask1)) {
                field.obflag[i] = Util::ibset(field.obflag[i], BB::driver, Util::Mask1, 0U);
                field.obflag[i] = Util::ibset(field.obflag[i], BB::semicyclic, Util::Mask1, 1U);
            }
        }
    }
    static void disable_driver(vector_field<real_t> &field) {
        for (int i = 0; i < 6; i ++) {
            if (Util::ibsee(field.obflag[i], BB::directional, Util::Mask1)) {
                field.obflag[i] = Util::ibset(field.obflag[i], BB::directional, Util::Mask1, 0U);
                field.obflag[i] = Util::ibset(field.obflag[i], BB::cyclic, Util::Mask1, 1U);
            }
            if (Util::ibsee(field.obflag[i], BB::driver, Util::Mask1)) {
                field.obflag[i] = Util::ibset(field.obflag[i], BB::driver, Util::Mask1, 0U);
                field.obflag[i] = Util::ibset(field.obflag[i], BB::semicyclic, Util::Mask1, 1U);
            }
        }
    }
public:
    static void scalar_outer(scalar_field<real_t> &fld, Dom &dom);
    static void vector_outer(vector_field<real_t> &fld, Dom &dom);
    static void scalar_outflow(scalar_field<real_t> &fld, Dom &dom);
    static void vector_outflow(vector_field<real_t> &fld, Dom &dom);
    static void velocity_outflow_correction(Dom &dom);
    static void calc_u_ob(Dom &dom);
};
#endif

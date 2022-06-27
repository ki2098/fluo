#ifndef _MMAC_H
#define _MMAC_H 1

#include "basic.h"
#include "domain.h"

class MMAC {
public:
    vector_field<real_t> ua, uc, up, uua, uup;
    scalar_field<real_t> pp, diva, divp;
public:
    void pseudo_velocity(vector_field<real_t> &u, vector_field<real_t> &uu, vector_field<real_t> &ua, scalar_field<real_t> &nut, Dom &dom);
    void correct_center_velocity(vector_field<real_t> &u, vector_field<real_t> &ua, scalar_field<real_t> &p, Dom &dom);
    void correct_face_velocity(vector_field<real_t> &u, vector_field<real_t> &uu, vector_field<real_t> &uua, scalar_field<real_t> &p, Dom &dom);
    void interpolate_velocity(vector_field<real_t> &u, vector_field<real_t> &uc, vector_field<real_t> &uu, Dom &dom);
    void divergence_velocity(vector_field<real_t> &uu, scalar_field<real_t> &div, real_t &div_norm, Dom &dom);
public:
    void to_device() {
        #pragma acc enter data copyin(this[0:1])
    }
    void end_device() {
        #pragma acc exit data delete(this[0:1])
    }
};

#endif

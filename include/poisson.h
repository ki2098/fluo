#ifndef _POISSON_H
#define _POISSON_H 1

#include "domain.h"
#include "scalarfield.h"
#include "vectorfield.h"

class Poisson {
public:
    scalar_field<real_t> pd;
public:
    void sor(scalar_field<real_t> &p, scalar_field<real_t> &div, Dom &dom);
    void jacobi(scalar_field<real_t> &p, scalar_field<real_t> &div, Dom &dom);
public:
    void to_device() {
        #pragma acc enter data copyin(this[0:1])
    }
    void end_device() {
        #pragma acc exit data delete(this[0:1])
    }
};

#endif

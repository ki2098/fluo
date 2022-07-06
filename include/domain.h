#ifndef _DOMAIN_H
#define _DOMAIN_H 1

#include "basic.h"
#include "util.h"
#include "control.h"
#include "scalarfield.h"
#include "vectorfield.h"

class Dom {
public:
    int    size[3];
    int    len;
    int    bcount;
    Ctrl   c;
    struct UOB {
        enum class Type {off, average, minmax, designated};
        Type   type[6];
        real_t value[6];
        UOB() {
            type[0] = type[1] = type[2] = type[3] = type[4] = type[5] = Type::off;
        }
    } uob;
public:
    scalar_field<unsigned> f;
    vector_field<real_t>   u;
    vector_field<real_t>   uu;
    scalar_field<real_t>   p;
    scalar_field<real_t>   nut;
    scalar_field<real_t>   div;
    vector_field<real_t>   x;
    vector_field<real_t>   kx;
    vector_field<real_t>   g;
    scalar_field<real_t>   ja;
    vector_field<real_t>   sz;
    vector_field<real_t>   u_avg;
    scalar_field<real_t>   p_avg;
public:
    const char *name;
public:
    Dom(const char *label);
    void init(int *field, int nbound);
    void driver();
    void driver_monitor();
    void pressure_zero_average();
    void tick() {c.time.idt ++;}
public:
    void to_device() {
        #pragma acc enter data copyin(this[0:1])
    }
    void end_device() {
        #pragma acc exit data delete(this[0:1])
    }
};

#endif

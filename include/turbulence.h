#ifndef _TURBULENCE_H
#define _TURBULENCE_H 1

#include "domain.h"

class Turbulence {
public:
    void smagorinsky(Dom &dom);
    void csm(Dom &dom);
};

#endif

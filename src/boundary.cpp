#include <stdlib.h>
#include "boundary.h"
#include "util.h"

BC::BC() {
    this->f0 = 0U;
    this->f1 = 0U;
    this->f2 = 0U;
    this->f3 = 0U;
    this->f4 = 0U;
    this->f5 = 0U;
    this->n = 0;
    this->b = NULL;
}

BC::~BC() {
    delete[] this->b;
}

void BC::pre(unsigned int m, double _0, double _1, double &ref, double &dis) {
    ref = (m)? _1 : _0;
    dis = 0.5 - m;
}

double BC::eva(unsigned int flag, double ref, double dis, double value) {
    unsigned int drc = Util::ibsee(flag, BD::Bd, Util::Mask1);
    unsigned int neu = Util::ibsee(flag, BD::Bn, Util::Mask1);
    if (drc) {
        return value;
    } else if (neu) {
        return ref + dis * value;
    }
    return 0.0;
}
#include <stdlib.h>
#include "domain.h"

D::D() {
    this->size[0] = 0;
    this->size[1] = 0;
    this->size[2] = 0;
    this->F = NULL;
}

D::~D() {
    delete[] this->F;
}
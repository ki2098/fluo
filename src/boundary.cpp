#include <stdlib.h>
#include "boundary.h"

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
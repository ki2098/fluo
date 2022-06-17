#include <stdlib.h>
#include "domain.h"

D::D() {
    this->size[0] = 0;
    this->size[1] = 0;
    this->size[2] = 0;
    this->F   = NULL;
    this->U   = NULL;
    this->UU  = NULL;
    this->UC  = NULL;
    this->UA  = NULL;
    this->UUA = NULL;
    this->UR  = NULL;
    this->P   = NULL;
    this->PD  = NULL;
    this->PR  = NULL;
    this->SGS = NULL;
    this->X   = NULL;
    this->KX  = NULL;
    this->G   = NULL;
    this->J   = NULL;
    this->monitor = 0U;
}

D::~D() {
    delete[] this->F;
    delete[] this->U;
    delete[] this->UU;
    delete[] this->UC;
    delete[] this->UA;
    delete[] this->UUA;
    delete[] this->UR;
    delete[] this->P;
    delete[] this->PD;
    delete[] this->PR;
    delete[] this->SGS;
    delete[] this->X;
    delete[] this->KX;
    delete[] this->G;
    delete[] this->J;
}
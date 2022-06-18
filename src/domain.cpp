#include <stdlib.h>
#include "domain.h"

D::D() {
    this->size[0] = 0;
    this->size[1] = 0;
    this->size[2] = 0;
    this->f   = NULL;
    this->u   = NULL;
    this->uu  = NULL;
    this->uc  = NULL;
    this->ua  = NULL;
    this->uua = NULL;
    this->ur  = NULL;
    this->p   = NULL;
    this->pd  = NULL;
    this->pr  = NULL;
    this->sgs = NULL;
    this->x   = NULL;
    this->kx  = NULL;
    this->g   = NULL;
    this->ja  = NULL;
    this->monitor = 0U;
}

D::~D() {
    delete[] this->f;
    delete[] this->u;
    delete[] this->uu;
    delete[] this->uc;
    delete[] this->ua;
    delete[] this->uua;
    delete[] this->ur;
    delete[] this->p;
    delete[] this->pd;
    delete[] this->pr;
    delete[] this->sgs;
    delete[] this->x;
    delete[] this->kx;
    delete[] this->g;
    delete[] this->ja;
}
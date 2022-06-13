#include <string.h>
#include "alloc.h"

double* Alloc::double_3d(int* size) {
    int nn = size[0] * size[1] * size[2];

    double* var = new double[nn];
    memset(var, 0, sizeof(double) * nn);

    return var;
}

double* Alloc::double_4d(int* size) {
    int nn = size[0] * size[1] * size[2] * size[3];

    double* var = new double[nn];
    memset(var, 0, sizeof(double) * nn);
    
    return var;
}

double* Alloc::double_5d(int* size) {
    int nn = size[0] * size[1] * size[2] * size[3] * size[4];

    double* var = new double[nn];
    memset(var, 0, sizeof(double) * nn);
    
    return var;
}

unsigned int* Alloc::uint_3d(int* size) {
    int nn = size[0] * size[1] * size[2];

    unsigned int* var = new unsigned int[nn];
    memset(var, 0, sizeof(unsigned int) * nn);

    return var;
}

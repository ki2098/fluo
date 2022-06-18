#ifndef _ALLOC_H
#define _ALLOC_H 1

class Alloc {
public:
    static double* double_3d(int* size);
    static double* double_4d(int* size);
    static double* double_5d(int* size);
    static unsigned int* uint_3d(int* size);
};

#endif
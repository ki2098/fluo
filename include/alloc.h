#ifndef _ALLOC_H
#define _ALLOC_H 1

#include <string.h>

class Alloc {
public:

template<class T> static T *alloc_array(int len) {
    T *var = new T[len];
    memset(var, 0, sizeof(T) * len);
    return var;
}

template<class T> static T *alloc_scalar(int *size) {
    int len = size[0] * size[1] * size[2];
    T *var = new T[len];
    memset(var, 0, sizeof(T) * len);
    return var;
}

template<class T> static T *alloc_vector(int *size) {
    int len = size[0] * size[1] * size[2] * size[3];
    T *var = new T[len];
    memset(var, 0, sizeof(T) * len);
    return var;
}

};

#endif

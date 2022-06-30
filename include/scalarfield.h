#ifndef _SCALARFIELD_H
#define _SCALARFIELD_H 1

#include "alloc.h"
#include "basic.h"

template<class T> class scalar_field {
public:
    int size[3];
    int len;
    T  *m;
public:
    int       bcount;
    int       blen;
    T        *b;
    unsigned *bflag;
public:
    unsigned obflag[6];
public:
    bool inited;
public:
    scalar_field();
    ~scalar_field();
public:
    void init(int *field, int n);
    void init(scalar_field<T> &field);
public:
    void to_device() {
        #pragma acc enter data copyin(this[0:1], m[0:len], b[0:blen], bflag[0:bcount])
    }
    void end_device() {
        #pragma acc exit data delete(m[0:len], b[0:blen], bflag[0:bcount], this[0:1])
    }
    void update_device() {
        #pragma acc update device(m[0:len])
    }
    void update_self() {
        #pragma acc update self(m[0:len])
    }
};

template<class T> scalar_field<T>::scalar_field() {
    size[0] = 0;
    size[1] = 0;
    size[2] = 0;
    len     = 0;
    m       = nullptr;

    bcount  = 0;
    blen    = 0;
    b       = nullptr;
    bflag   = nullptr;

    for (int i = 0; i < 6; i ++) {
        obflag[i] = 0U;
    }

    inited = false;
}

template<class T> scalar_field<T>::~scalar_field() {
    delete[] m;
    delete[] b;
    delete[] bflag;
}

template<class T> void scalar_field<T>::init(int *field, int n) {
    if (!inited) {
        size[0]   = field[0];
        size[1]   = field[1];
        size[2]   = field[2];
        len       = size[0] * size[1] * size[2];
        m         = Alloc::alloc_array<T>(len);

        bcount    = n;
        blen      = n;
        b         = Alloc::alloc_array<T>(blen);
        bflag     = Alloc::alloc_array<unsigned>(bcount);

        for (int i = 0; i < 6; i ++) {
            obflag[i] = 0U;
        }
        
        inited    = true;
    }
}

template<class T> void scalar_field<T>::init(scalar_field<T> &field) {
    if (!inited) {
        size[0]   = field.size[0];
        size[1]   = field.size[1];
        size[2]   = field.size[2];
        len       = field.len;
        m         = Alloc::alloc_array<T>(len);

        bcount    = field.bcount;
        blen      = field.blen;
        b         = Alloc::alloc_array<T>(blen);
        for (int i = 0; i < blen; i ++) {
            b[i] = field.b[i];
        }
        bflag     = Alloc::alloc_array<unsigned>(bcount);
        for (int i = 0; i < bcount; i ++) {
            bflag[i] = field.bflag[i];
        }

        for (int i = 0; i < 6; i ++) {
            obflag[i] = field.obflag[i];
        }
        
        inited    = true;
    }
}

#endif

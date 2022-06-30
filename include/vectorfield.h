#ifndef _VECTORFIELD_H
#define _VECTORFIELD_H 1

#include "alloc.h"
#include "basic.h"

template<class T> class vector_field {
public:
    int size[4];
    int len;
    T  *m;
public:
    int       bcount;
    int       bsize[2];
    int       blen;
    T        *b;
    unsigned *bflag;
public:
    unsigned obflag[6];
public:
    bool inited;
public:
    vector_field();
    ~vector_field();
public:
    void init(int *field, int dim, int n);
    void init(vector_field<T> &field);
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

template<class T> vector_field<T>::vector_field() {
    size[0]   = 0;
    size[1]   = 0;
    size[2]   = 0;
    size[3]   = 0;
    len       = 0;
    m         = nullptr;

    bcount    = 0;
    bsize[0]  = 0;
    bsize[1]  = 0;
    blen      = 0;
    b         = nullptr;
    bflag     = nullptr;

    for (int i = 0; i < 6; i ++) {
        obflag[i] = 0U;
    }

    inited    = false;
}

template<class T> vector_field<T>::~vector_field() {
    delete[] m;
    delete[] b;
    delete[] bflag;
}

template<class T> void vector_field<T>::init(int *field, int dim, int n) {
    if (!inited) {
        size[0]   = field[0];
        size[1]   = field[1];
        size[2]   = field[2];
        size[3]   = dim;
        len       = size[0] * size[1] * size[2] * size[3];
        m         = Alloc::alloc_array<T>(len);

        bcount    = n;
        bsize[0]  = n;
        bsize[1]  = dim;
        blen      = bsize[0] * bsize[1];
        b         = Alloc::alloc_array<T>(blen);
        bflag     = Alloc::alloc_array<unsigned>(bcount);

        for (int i = 0; i < 6; i ++) {
            obflag[i] = 0U;
        }

        inited    = true;
    }
}

template<class T> void vector_field<T>::init(vector_field<T> &field) {
    if (!inited) {
        size[0]  = field.size[0];
        size[1]  = field.size[1];
        size[2]  = field.size[2];
        size[3]  = field.size[3];
        len      = field.len;
        m        = Alloc::alloc_vector<T>(size);

        bcount   = field.bcount;
        bsize[0] = field.bsize[0];
        bsize[1] = field.bsize[1];
        blen     = field.blen;
        b        = Alloc::alloc_array<T>(blen);
        for (int i = 0; i < blen; i ++) {
            b[i] = field.b[i];
        }
        bflag    = Alloc::alloc_array<unsigned>(bcount);
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

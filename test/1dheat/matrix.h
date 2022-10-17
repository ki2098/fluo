#ifndef _MATRIX_H
#define _MATRIX_H 1

#include <string.h>
#include <math.h>
#include "type.h"
#include "dom.h"
#include "util.h"

template<class T>
class Mat {
public:
    T  *mat;
    int row;
    int col;
    int num;
public:
    Mat(int row, int col);
    Mat(Dom &dom, int col);
    T &get(int i, int j);
    T &get(int i);
    void to_device();
    void off_device();
    void update_device();
    void update_self();
    ~Mat();
};

template<class T>
Mat<T>::Mat(int row, int col) : row(row), col(col) {
    num = row * col;
    mat = new T[num];
    memset(mat, 0, num * sizeof(T));
}

template<class T>
Mat<T>::Mat(Dom &dom, int col) : row(dom.num), col(col) {
    num = row * col;
    mat = new T[num];
    memset(mat, 0, num * sizeof(T));
}

template<class T>
inline T &Mat<T>::get(int i, int j) {
    return mat[i * col + j];
}

template<class T>
inline T &Mat<T>::get(int i) {
    return mat[i];
}

template<class T>
void Mat<T>::to_device() {
    #pragma acc enter data copyin(this[0:1], mat[0:num])
}

template<class T>
void Mat<T>::off_device() {
    #pragma acc exit data delete(mat[0:num], this[0:1])
}

template<class T>
void Mat<T>::update_device() {
    #pragma acc update device(mat[0:num])
}

template<class T>
void Mat<T>::update_self() {
    #pragma acc update self(mat[0:num])
}

template<class T>
Mat<T>::~Mat() {
    delete[] mat;
}

#endif
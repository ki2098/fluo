#ifndef _LS_H
#define _LS_H 1

class LS {
public:
    int type;
    enum Type{jacobi, sor, bicgstab, pbicgstab};
    double omega;
};

#endif
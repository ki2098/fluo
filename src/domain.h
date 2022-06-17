#ifndef _DOMAIN_H
#define _DOMAIN_H 1

class D {
public:
    static const int GUIDE = 2;
public:
    int size[3];
    unsigned int *F;
    double       *U;
    double       *UU;
    double       *UC;
    double       *UA;
    double       *UUA;
    double       *UR;
    double       *P;
    double       *PD;
    double       *PR;
    double       *SGS;
    double       *X;
    double       *KX;
    double       *G;
    double       *J;
    double        re;
    double        ri;
    double        dt;
    int           ntime;
    double        tdiv;
    unsigned int  monitor;
    int           monitor_interval;
public:
    D();
    ~D();
};

#endif
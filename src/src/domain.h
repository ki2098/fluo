#ifndef _DOMAIN_H
#define _DOMAIN_H 1

class D {
public:
    static const int GUIDE = 2;
public:
    int size[3];
    unsigned int *f;
    double       *u;
    double       *uu;
    double       *uc;
    double       *ua;
    double       *uua;
    double       *ur;
    double       *p;
    double       *pd;
    double       *pr;
    double       *sgs;
    double       *x;
    double       *kx;
    double       *g;
    double       *ja;
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
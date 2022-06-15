#ifndef _BOUNDARY_H
#define _BOUNDARY_H 1

class BD {
public:
    static const unsigned int Ud = 0U;
    static const unsigned int Un = 1U;
    static const unsigned int Pd = 6U;
    static const unsigned int Pn = 7U;
    static const unsigned int Nd = 8U;
    static const unsigned int Nn = 9U;
    static const unsigned int Ub = Ud;
    static const unsigned int Pb = Pd;
    static const unsigned int Nb = Nd;
    static const unsigned int Bd = 0U;
    static const unsigned int Bn = 1U;
public:
    static const unsigned int WALLFLux = 20U;
    static const unsigned int WallNS   = 21U;
    static const unsigned int WallFunc = 22U;
public:
    unsigned int flag;
    double u;
    double v; 
    double w; 
    double p; 
    double nt; 
};

class BC {
public:
    unsigned int f0, f1, f2, f3, f4, f5;
    enum OBC{US, UP, UO, PS, PP, NS, NP};
public:
    int n;
    BD* b;
public:
    BC();
    ~BC();
};

#endif
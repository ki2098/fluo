#ifndef _DRIVER_H 
#define _DRIVER_H 1

class Driver{
public:
    unsigned int on;
    unsigned int direction;
    unsigned int type;
    double       value;
    double       dp;
public:
    enum Type{U, P};
};

#endif
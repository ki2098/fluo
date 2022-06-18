#ifndef _UTIL_H
#define _UTIL_H 1

#define id3(i,j,k,size) ((i)*size[1]*size[2] + (j)*size[2] + (k))
#define id4(i,j,k,m,size) ((i)*size[1]*size[2]*size[3] + (j)*size[2]*size[3] + (k)*size[3] + (m))
#define id5(i,j,k,m,n,size) ((i)*size[1]*size[2]*size[3]*size[4] + (j)*size[2]*size[3]*size[4] + (k)*size[3]*size[4] + (m)*size[4] +(n))

class Util {
public:
    static const unsigned int Mask1 = 1U;
    static const unsigned int Mask2 = 3U;
    static const unsigned int Mask4 = 15U;
    static const unsigned int Mask8 = 255U;
public:
    static unsigned int ibsee(unsigned int bits, unsigned int i, unsigned int mask);
    static unsigned int ibset(unsigned int bits, unsigned int i, unsigned int mask, unsigned int value);
    static double max(double a, double b);
    static double min(double a, double b);
};

#endif

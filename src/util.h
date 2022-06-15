#ifndef _UTIL_H
#define _UTIL_H 1

#define id3(i,j,k,len) ((i)*len[1]*len[2] + (j)*len[2]+(k))
#define id4(i,j,k,m,len) ((i)*len[1]*len[2]*len[3] + (j)*len[2]*len[3] + (k)*len[3] + (m))
#define id5(i,j,k,m,n,len) ((i)*len[1]*len[2]*len[3]*len[4] + (j)*len[2]*len[3]*len[4] + (k)*len[3]*len[4] + (m)*len[4] +(n))

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

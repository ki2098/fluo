#ifndef _UTIL_H
#define _UTIL_H 1

#define id2(i,j,size) ((i)*(size[1]) + (j))
#define id3(i,j,k,size) ((i)*(size[1])*(size[2]) + (j)*(size[2]) + (k))
#define id4(i,j,k,m,size) ((i)*(size[1])*(size[2])*(size[3]) + (j)*(size[2])*(size[3]) + (k)*(size[3]) + (m))

class Util {
public:

static const unsigned Mask1 = 1U;
static const unsigned Mask2 = 3U;
static const unsigned Mask4 = 15U;
static const unsigned Mask8 = 255U;

static unsigned ibsee(unsigned int bits, unsigned int i, unsigned int mask) {
    bits = (bits >> i) & mask;
    return bits;
}

static unsigned ibset(unsigned int bits, unsigned int i, unsigned int mask, unsigned int value) {
    bits = bits & ~(mask << i);
    bits = bits | (value << i);
    return bits;
}

static double max(double a, double b) {
    return ((a > b)? a : b);
}

static double min(double a, double b) {
    return ((a < b)? a : b);
}

};

#endif

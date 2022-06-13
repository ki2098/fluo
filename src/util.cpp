#include "util.h"

unsigned int Util::ibsee(unsigned int bits, unsigned int i, unsigned int mask) {
    bits = (bits >> i) & mask;
    return bits;
}

unsigned int Util::ibset(unsigned int bits, unsigned int i, unsigned int mask, unsigned int value) {
    bits = bits & ~(mask << i);
    bits = bits | (value << i);
    return bits;
}

double Util::max(double a, double b) {
    return ((a > b)? a : b);
}

double Util::min(double a, double b) {
    return ((a < b)? a : b);
}

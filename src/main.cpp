#include <stdio.h>
#include "setup.h"
#include "basic.h"

int main(void) {
    Setup ss((char*)"./setup/setup.json");
    printf("%d %d %d\n", ss.N[0], ss.N[1], ss.N[2]);

    int size[3] = {
        ss.N[0] + 2 * GUIDE,
        ss.N[1] + 2 * GUIDE,
        ss.N[2] + 2 * GUIDE
    };

    unsigned int* F = ss.setup_flag();
    printf("%d %d %d\n", F[idx3d(0,0,0,size)], F[idx3d(GUIDE,GUIDE,GUIDE,size)], F[idx3d(92,52,12,size)]);
}
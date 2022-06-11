#include <stdio.h>
#include "setup.h"

int main(void) {
    Setup ss((char*)"./setup/setup.json");
    printf("%d %d %d\n", ss.N[0], ss.N[1], ss.N[2]);
}
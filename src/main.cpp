#include <stdio.h>
#include <string.h>
#include <time.h>
#include "../include/flag.h"
#include "../include/fluo.h"

int main(int argc, char **argv) {
    FLUO fluo("./setting/setting.json");
    if (argc == 1 || argc == 2 && strcmp(argv[1], "full") == 0) {
        fluo.init();
        fluo.show_info();
        fluo.param_out();
    } else if (argc == 2 && strcmp(argv[1], "setting") == 0) {
        fluo.init(false);
        fluo.show_info();
        return 0;
    } else if (argc == 2 && strcmp(argv[1], "mesh") == 0) {
        fluo.init();
        fluo.show_info();
        fluo.param_out();
        return 0;
    } else if (argc == 2 && strcmp(argv[1], "0") == 0) {
        fluo.init();
        fluo.show_info();
        fluo.param_out();
        fluo.domain.c.time.ndt = 0;
    }
    
    clock_t start, end;
    start = clock();

    fluo.fractional_step();
    
    end = clock();
    printf("cpu time %f sec\n", double(end - start) / CLOCKS_PER_SEC);

    return 0;
}

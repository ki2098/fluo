#include <stdio.h>
#include "setup.h"
#include "alloc.h"
#include "basic.h"



Setup::Setup(char* fname) {
    this->doc = yyjson_read_file(fname, 0, NULL, NULL);
    this->root = yyjson_doc_get_root(this->doc);

    yyjson_val *domain = yyjson_obj_get(this->root, "domain");
    this->N[0] = yyjson_get_uint(yyjson_arr_get(domain, 0));
    this->N[1] = yyjson_get_uint(yyjson_arr_get(domain, 1));
    this->N[2] = yyjson_get_uint(yyjson_arr_get(domain, 2));
}

Setup::~Setup() {
    yyjson_doc_free(this->doc);
}

unsigned int* Setup::setup_flag(void) {
    int nn[3] = {
        this->N[0] + 2 * GUIDE,
        this->N[1] + 2 * GUIDE,
        this->N[2] + 2 * GUIDE,
    };

    unsigned int* F = Alloc::uint_3d(nn);

    printf("F allocated\n");

    for (int i = GUIDE; i < this->N[0] + GUIDE; i ++) {
        for (int j = GUIDE; j < this->N[1] + GUIDE; j ++) {
            for (int k = GUIDE; k < this->N[2] + GUIDE; k ++) {
                F[idx3d(i,j,k,nn)] = 1;
            }
        }
    }

    yyjson_val *boundary = yyjson_obj_get(this->root, "boundary");
    yyjson_val *b;
    yyjson_arr_iter iter;
    yyjson_arr_iter_init(boundary, &iter);
    while(b = yyjson_arr_iter_next(&iter)) {
        if(yyjson_equals_str(yyjson_obj_get(b, "type"), "obstacle")) {
            int origin[3];
            int size[3];

            yyjson_val* arr = yyjson_obj_get(b, "origin");
            origin[0] = yyjson_get_uint(yyjson_arr_get(arr, 0));
            origin[1] = yyjson_get_uint(yyjson_arr_get(arr, 1));
            origin[2] = yyjson_get_uint(yyjson_arr_get(arr, 2));

            arr = yyjson_obj_get(b, "size");
            size[0] = yyjson_get_uint(yyjson_arr_get(arr, 0));
            size[1] = yyjson_get_uint(yyjson_arr_get(arr, 1));
            size[2] = yyjson_get_uint(yyjson_arr_get(arr, 2));

            for (int i = origin[0]; i < origin[0] + size[0]; i ++) {
                for (int j = origin[1]; j < origin[1] + size[1]; j ++) {
                    for (int k = origin[2]; k < origin[2] + size[2]; k ++) {
                        F[idx3d(i,j,k,nn)] = 0;
                    }
                }
            }
        }
    }

    printf("F setted\n");

    return F;
}
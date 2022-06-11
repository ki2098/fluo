#include "setup.h"
#include "alloc.h"
#include "basic.h"

Setup::Setup(char* fname) {
    this->doc = yyjson_read_file(fname, 0, NULL, NULL);
    this->root = yyjson_doc_get_root(this->doc);

    yyjson_val *domain = yyjson_obj_get(this->root, "domain");
    yyjson_val *nn;
    yyjson_arr_iter iter;
    int id = 0;
    yyjson_arr_iter_init(domain, &iter);
    while (nn = yyjson_arr_iter_next(&iter)) {
        this->N[id] = yyjson_get_sint(nn);
        id ++;
    }
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

    
}
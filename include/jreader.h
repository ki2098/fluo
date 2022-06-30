#ifndef _JREADER_H
#define _JREADER_H 1

#include "../lib/yyjson.h"
#include "domain.h"

class JReader {
public:
    yyjson_doc *doc;
    yyjson_val *root;
public:
    JReader(const char *fname);
    ~JReader();
    void load_domain(Dom &dom, bool load_mesh = true);
private:
    void load_bcfi(Dom &dom, vector_field<real_t> &var, unsigned flag, yyjson_val *topo);
    void load_bcfi(Dom &dom, scalar_field<real_t> &var, unsigned flag, yyjson_val *topo);
    void load_uobi(Dom &dom, Dom::UOB::Type type, yyjson_val *topo);
    void load_uobi(Dom &dom, real_t value, yyjson_val *topo);
    int search_boundary(yyjson_val *boundary, const char *label);
};

#endif

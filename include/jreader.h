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
};

#endif

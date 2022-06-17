#ifndef _JREADER_H
#define _JREADER_H 1

#include "../lib/yyjson.h"
#include "boundary.h"
#include "driver.h"
#include "domain.h"
#include "ls.h"

class JReader {
public:
    yyjson_doc *doc;
    yyjson_val *root;
public:
    JReader(char *fname);
    ~JReader();
    void load_bc(BC *bc, Driver *driver);
    void load_domain(D *dom, LS *ls);
    void load_mesh(D *dom);
};

#endif
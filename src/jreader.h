#ifndef _JREADER_H
#define _JREADER_H 1

#include "../lib/yyjson.h"
#include "boundary.h"
#include "driver.h"
#include "domain.h"

class JReader {
public:
    yyjson_doc* doc;
    yyjson_val* root;
public:
    JReader(char* fname);
    ~JReader();
    void setup_bc(BC* bc, Driver* driver);
    void setup_domain(D* dom);
};

#endif
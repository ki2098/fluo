#ifndef _SETUP_H
#define _SETUP_H 1

#include "../lib/yyjson.h"

class Setup{
public:
    yyjson_doc *doc;
    yyjson_val *root;
    int N[3];
public:
    Setup(char* fname);
    ~Setup();
    unsigned int* setup_flag(void);
};

#endif

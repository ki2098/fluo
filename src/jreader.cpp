#include <stdio.h>
#include "jreader.h"
#include "util.h"
#include "alloc.h"
#include "flag.h"

void load_bcfi(BC* bc, unsigned int flag, yyjson_val* topo) {
    yyjson_val *io = yyjson_obj_get(topo, "io");
    if (yyjson_equals_str(io, "outer")) {
        yyjson_val *faces = yyjson_obj_get(topo, "face");
        yyjson_arr_iter iter;
        yyjson_arr_iter_init(faces, &iter);
        yyjson_val *face;
        while (face = yyjson_arr_iter_next(&iter)) {
            if (yyjson_equals_str(face, "x-")) {
                bc->f0 = Util::ibset(bc->f0, flag, Util::Mask1, 1);
            } else if (yyjson_equals_str(face, "x+")) {
                bc->f1 = Util::ibset(bc->f1, flag, Util::Mask1, 1);
            } else if (yyjson_equals_str(face, "y-")) {
                bc->f2 = Util::ibset(bc->f2, flag, Util::Mask1, 1);
            } else if (yyjson_equals_str(face, "y+")) {
                bc->f3 = Util::ibset(bc->f3, flag, Util::Mask1, 1);
            } else if (yyjson_equals_str(face, "z-")) {
                bc->f4 = Util::ibset(bc->f4, flag, Util::Mask1, 1);
            } else if (yyjson_equals_str(face, "z+")) {
                bc->f5 = Util::ibset(bc->f5, flag, Util::Mask1, 1);
            }
        }
    }
}

JReader::JReader(char* fname) {
    this->doc = yyjson_read_file(fname, 0, NULL, NULL);
    this->root = yyjson_doc_get_root(this->doc);
}

JReader::~JReader() {
    yyjson_doc_free(this->doc);
}

void JReader::setup_bc(BC* bc, Driver* driver) {
    yyjson_val *boundary = yyjson_obj_get(this->root, "boundary");
    bc->n = yyjson_arr_size(boundary);
    bc->b = new BD[bc->n + 1];
    yyjson_val *bound;
    yyjson_arr_iter iter;
    yyjson_arr_iter_init(boundary, &iter);
    int index = 1;
    while (bound = yyjson_arr_iter_next(&iter)) {
        yyjson_val *var, *type, *value;
        unsigned int flag = 0U;
        var = yyjson_obj_get(bound, "U");
        if (var) {
            type = yyjson_obj_get(var, "type");
            if (yyjson_equals_str(type, "fixedValue")) {
                flag = Util::ibset(flag, BD::Ud, Util::Mask1, 1U);
            } else if (yyjson_equals_str(type, "fixedGradiet")) {
                flag = Util::ibset(flag, BD::Un, Util::Mask1, 1U);
            } else if (yyjson_equals_str(type, "symmetric")) {
                load_bcfi(bc, BC::OBC::US, yyjson_obj_get(bound, "topo"));
            } else if (yyjson_equals_str(type, "periodic")) {
                load_bcfi(bc, BC::OBC::UP, yyjson_obj_get(bound, "topo"));
            } else if (yyjson_equals_str(type, "outlofw")) {
                load_bcfi(bc, BC::OBC::UO, yyjson_obj_get(bound, "topo"));
            }
            value = yyjson_obj_get(var, "value");
            bc->b[index].u = yyjson_get_real(yyjson_arr_get(value, 0));
            bc->b[index].v = yyjson_get_real(yyjson_arr_get(value, 1));
            bc->b[index].w = yyjson_get_real(yyjson_arr_get(value, 2));
        }
        var = yyjson_obj_get(bound, "P");
        if (var) {
            type = yyjson_obj_get(var, "type");
            if (yyjson_equals_str(type, "fixedValue")) {
                flag = Util::ibset(flag, BD::Pd, Util::Mask1, 1U);
            } else if (yyjson_equals_str(type, "fixedGradient")) {
                flag = Util::ibset(flag, BD::Pn, Util::Mask1, 1U);
            } else if (yyjson_equals_str(type, "symmetric")) {
                load_bcfi(bc, BC::OBC::PS, yyjson_obj_get(bound, "topo"));
            } else if (yyjson_equals_str(type, "periodic")) {
                load_bcfi(bc, BC::OBC::PP, yyjson_obj_get(bound, "topo"));
            } else if (yyjson_equals_str(type, "directional")) {
                driver->on = 1U;
                load_bcfi(bc, BC::OBC::PP, yyjson_obj_get(bound, "topo"));
                yyjson_val *reference = yyjson_obj_get(var,"U");
                int refvalue;
                if (reference) {
                    driver->type = Driver::Type::U;
                    refvalue = yyjson_get_real(reference);
                }
                reference = yyjson_obj_get(var, "P");
                if (reference) {
                    driver->type = Driver::Type::P;
                    refvalue = yyjson_get_real(reference);
                }
                yyjson_val* direction = yyjson_obj_get(var, "direction");
                if (yyjson_equals_str(direction, "x+")) {
                    driver->direction = 0U;
                    driver->value = refvalue;
                } else if (yyjson_equals_str(direction, "x-")) {
                    driver->direction = 0U;
                    driver->value = - refvalue;
                } else if (yyjson_equals_str(direction, "y+")) {
                    driver->direction = 1U;
                    driver->value = refvalue;
                } else if (yyjson_equals_str(direction, "y-")) {
                    driver->direction = 1U;
                    driver->value = - refvalue;
                } else if (yyjson_equals_str(direction, "z+")) {
                    driver->direction = 2U;
                    driver->value = refvalue;
                } else if (yyjson_equals_str(direction, "z-")) {
                    driver->direction = 2U;
                    driver->value = - refvalue;
                }
                
            }
            value = yyjson_obj_get(var, "value");
            bc->b[index].p = yyjson_get_real(value);
        }
        var = yyjson_obj_get(bound, "Nt");
        if (var) {
            type = yyjson_obj_get(var, "type");
            if (yyjson_equals_str(type, "fixedValue")) {
                flag = Util::ibset(flag, BD::Nd, Util::Mask1, 1U);
            } else if (yyjson_equals_str(type, "fixedGradient")) {
                flag = Util::ibset(flag, BD::Nn, Util::Mask1, 1U);
            } else if (yyjson_equals_str(type, "symmetric")) {
                load_bcfi(bc, BC::OBC::NS, yyjson_obj_get(bound, "topo"));
            } else if (yyjson_equals_str(type, "periodic")) {
                load_bcfi(bc, BC::OBC::NP, yyjson_obj_get(bound, "topo"));
            }
            value = yyjson_obj_get(var, "value");
            bc->b[index].nt = yyjson_get_real(value);
        }
        bc->b[index].flag = flag;
        index ++;
    }
}

void JReader::setup_domain(D *dom) {
    yyjson_val *domain = yyjson_obj_get(this->root, "domain");
    dom->size[0] = yyjson_get_int(yyjson_arr_get(domain, 0)) + 2 * D::GUIDE;
    dom->size[1] = yyjson_get_int(yyjson_arr_get(domain, 1)) + 2 * D::GUIDE;
    dom->size[2] = yyjson_get_int(yyjson_arr_get(domain, 2)) + 2 * D::GUIDE;
    dom->F = Alloc::uint_3d(dom->size);
    int *size = dom->size;
    unsigned int *F = dom->F;

    for (int i = D::GUIDE; i < size[0] - D::GUIDE; i ++) {
        for (int j = D::GUIDE; j < size[1] - D::GUIDE; j ++) {
            for (int k = D::GUIDE; k < size[2] - D::GUIDE; k ++) {
                F[id3(i,j,k,size)] = Util::ibset(F[id3(i,j,k,size)], Flag::Active, Util::Mask1, 1);
            }
        }
    }

    yyjson_val *boundary = yyjson_obj_get(this->root, "boundary");
    yyjson_val *bound;
    yyjson_arr_iter iter;
    yyjson_arr_iter_init(boundary, &iter);
    int index = 1;
    while (bound = yyjson_arr_iter_next(&iter)) {
        yyjson_val *topo = yyjson_obj_get(bound, "topo");
        yyjson_val *io = yyjson_obj_get(topo, "io");
        if (yyjson_equals_str(io, "outer")) {
            yyjson_val *faces = yyjson_obj_get(topo, "face");
            yyjson_arr_iter iter;
            yyjson_arr_iter_init(faces, &iter);
            yyjson_val *face;
            while (face = yyjson_arr_iter_next(&iter)) {
                if (yyjson_equals_str(face, "x+")) {
                    for (int j = D::GUIDE; j < size[1] - D::GUIDE; j ++) {
                        for (int k = D::GUIDE; k < size[2] - D::GUIDE; k ++) {
                            int i = size[0] - D::GUIDE - 1;
                            unsigned int flag = F[id3(i,j,k,size)];
                            flag = Util::ibset(flag, Flag::Fe, Util::Mask8, index);
                            flag = Util::ibset(flag, Flag::Me, Util::Mask1, 0);
                            F[id3(i,j,k,size)] = flag;
                        }
                    }
                } else if (yyjson_equals_str(face, "x-")) {
                    for (int j = D::GUIDE; j <size[1] - D::GUIDE; j ++) {
                        for (int k = D::GUIDE; k < size[2] - D::GUIDE; k ++) {
                            int i = D::GUIDE - 1;
                            unsigned int flag = F[id3(i,j,k,size)];
                            flag = Util::ibset(flag, Flag::Fe, Util::Mask8, index);
                            flag = Util::ibset(flag, Flag::Me, Util::Mask1, 1);
                            F[id3(i,j,k,size)] = flag;
                        }
                    }
                } else if (yyjson_equals_str(face, "y+")) {
                    for (int i = D::GUIDE; i < size[0] - D::GUIDE; i ++) {
                        for (int k = D::GUIDE; k < size[2] - D::GUIDE; k ++) {
                            int j = size[1] - D::GUIDE - 1;
                            unsigned int flag = F[id3(i,j,k,size)];
                            flag = Util::ibset(flag, Flag::Fn, Util::Mask8, index);
                            flag = Util::ibset(flag, Flag::Mn, Util::Mask1, 0);
                            F[id3(i,j,k,size)] = flag;
                        }
                    }
                } else if (yyjson_equals_str(face, "y-")) {
                    for (int i = D::GUIDE; i < size[0] - D::GUIDE; i ++) {
                        for (int k = D::GUIDE; k < size[2] - D::GUIDE; k ++) {
                            int j = D::GUIDE - 1;
                            unsigned int flag = F[id3(i,j,k,size)];
                            flag = Util::ibset(flag, Flag::Fn, Util::Mask8, index);
                            flag = Util::ibset(flag, Flag::Mn, Util::Mask1, 1);
                            F[id3(i,j,k,size)] = flag;
                        }
                    }
                } else if (yyjson_equals_str(face, "z+")) {
                    for (int i = D::GUIDE; i < size[0] - D::GUIDE; i ++) {
                        for (int j = D::GUIDE; j < size[1] - D::GUIDE; j ++) {
                            int k = size[2] - D::GUIDE - 1;
                            unsigned int flag = F[id3(i,j,k,size)];
                            flag = Util::ibset(flag, Flag::Ft, Util::Mask8, index);
                            flag = Util::ibset(flag, Flag::Mt, Util::Mask1, 0);
                            F[id3(i,j,k,size)] = flag;
                        }
                    }
                } else if (yyjson_equals_str(face, "z-")) {
                    for (int i = D::GUIDE; i < size[0] - D::GUIDE; i ++) {
                        for (int j = D::GUIDE; j < size[1] - D::GUIDE; j ++) {
                            int k = D::GUIDE - 1;
                            unsigned int flag = F[id3(i,j,k,size)];
                            flag = Util::ibset(flag, Flag::Ft, Util::Mask8, index);
                            flag = Util::ibset(flag, Flag::Mt, Util::Mask1, 1);
                            F[id3(i,j,k,size)] = flag;
                        }
                    }
                }
            }
        } else if (yyjson_equals_str(io, "inner")) {
            yyjson_val *shape = yyjson_obj_get(topo, "shape");
            if (yyjson_equals_str(shape, "cuboid")) {
                yyjson_val *position = yyjson_obj_get(topo, "position");
                int p1[3], p2[3];
                p1[0] = yyjson_get_int(yyjson_arr_get(yyjson_arr_get(position, 0), 0));
                p1[1] = yyjson_get_int(yyjson_arr_get(yyjson_arr_get(position, 0), 1));
                p1[2] = yyjson_get_int(yyjson_arr_get(yyjson_arr_get(position, 0), 2));
                p2[0] = yyjson_get_int(yyjson_arr_get(yyjson_arr_get(position, 1), 0));
                p2[1] = yyjson_get_int(yyjson_arr_get(yyjson_arr_get(position, 1), 1));
                p2[2] = yyjson_get_int(yyjson_arr_get(yyjson_arr_get(position, 1), 2));
                for (int i = p1[0] + D::GUIDE; i <= p2[0] + D::GUIDE; i ++) {
                    for (int j = p1[1] + D::GUIDE; j <= p2[1] + D::GUIDE; j ++) {
                        for (int k = p1[2] + D::GUIDE; k <= p2[2] + D::GUIDE; k ++) {
                            F[id3(i,j,k,size)] = Util::ibset(F[id3(i,j,k,size)], Flag::Active, Util::Mask1, 0);
                        }
                    }
                }
                if (yyjson_obj_get(bound, "U") || yyjson_obj_get(bound, "P")) {
                    for (int j = p1[1]  + D::GUIDE; j <= p2[1] + D::GUIDE; j ++) {
                        for (int k = p1[2] + D::GUIDE; k <= p2[2] + D::GUIDE; k ++) {
                            int i = p1[0] - 1 + D::GUIDE;
                            unsigned int flag = F[id3(i,j,k,size)];
                            flag = Util::ibset(flag, Flag::Fe, Util::Mask8, index);
                            flag = Util::ibset(flag, Flag::Me, Util::Mask1, 0);
                            F[id3(i,j,k,size)] = flag;

                            i = p2[0] + D::GUIDE;
                            flag = F[id3(i,j,k,size)];
                            flag = Util::ibset(flag, Flag::Fe, Util::Mask8, index);
                            flag = Util::ibset(flag, Flag::Me, Util::Mask1, 1);
                            F[id3(i,j,k,size)] = flag;
                        }
                    }
                    for (int i = p1[0] + D::GUIDE; i <= p2[0] + D::GUIDE; i ++) {
                        for (int k = p1[2] + D::GUIDE; k <= p2[2] + D::GUIDE; k ++) {
                            int j = p1[1] - 1 + D::GUIDE;
                            unsigned int flag = F[id3(i,j,k,size)];
                            flag = Util::ibset(flag, Flag::Fn, Util::Mask8, index);
                            flag = Util::ibset(flag, Flag::Mn, Util::Mask1, 0);
                            F[id3(i,j,k,size)] = flag;

                            j = p2[1] + D::GUIDE;
                            flag = F[id3(i,j,k,size)];
                            flag = Util::ibset(flag, Flag::Fn, Util::Mask8, index);
                            flag = Util::ibset(flag, Flag::Mn, Util::Mask1, 1);
                            F[id3(i,j,k,size)] = flag;
                        }
                    }
                    for (int i = p1[0] + D::GUIDE; i <= p2[0] + D::GUIDE; i ++) {
                        for (int j = p1[1] + D::GUIDE; j <= p2[1] + D::GUIDE; j ++) {
                            int k = p1[2] - 1 + D::GUIDE;
                            unsigned int flag = F[id3(i,j,k,size)];
                            flag = Util::ibset(flag, Flag::Ft, Util::Mask8, index);
                            flag = Util::ibset(flag, Flag::Mt, Util::Mask1, 0);
                            F[id3(i,j,k,size)] = flag;

                            k = p2[2] + D::GUIDE;
                            flag = F[id3(i,j,k,size)];
                            flag = Util::ibset(flag, Flag::Ft, Util::Mask8, index);
                            flag = Util::ibset(flag, Flag::Mt, Util::Mask1, 1);
                            F[id3(i,j,k,size)] = flag;
                        }
                    }
                }
            }
        }
        index ++;
    }
}
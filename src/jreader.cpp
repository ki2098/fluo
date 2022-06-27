#include <stdio.h>
#include "../include/jreader.h"
#include "../include/basic.h"
#include "../include/boundary.h"
#include "../include/flag.h"

void load_bcfi(vector_field<real_t> &var, unsigned flag, yyjson_val *topo) {
    yyjson_val *io = yyjson_obj_get(topo, "io");
    if (yyjson_equals_str(io, "outer")) {
        yyjson_val *faces = yyjson_obj_get(topo, "face");
        if (yyjson_is_arr(faces)) {
            yyjson_arr_iter iter;
            yyjson_arr_iter_init(faces, &iter);
            yyjson_val *face;
            while (face = yyjson_arr_iter_next(&iter)) {
                if (yyjson_equals_str(face, "x-")) {
                    var.obflag[0] = Util::ibset(var.obflag[0], flag, Util::Mask1, 1U);
                } else if (yyjson_equals_str(face, "x+")) {
                    var.obflag[1] = Util::ibset(var.obflag[1], flag, Util::Mask1, 1U);
                } else if (yyjson_equals_str(face, "y-")) {
                    var.obflag[2] = Util::ibset(var.obflag[2], flag, Util::Mask1, 1U);
                } else if (yyjson_equals_str(face, "y+")) {
                    var.obflag[3] = Util::ibset(var.obflag[3], flag, Util::Mask1, 1U);
                } else if (yyjson_equals_str(face, "z-")) {
                    var.obflag[4] = Util::ibset(var.obflag[4], flag, Util::Mask1, 1U);
                } else if (yyjson_equals_str(face, "z+")) {
                    var.obflag[5] = Util::ibset(var.obflag[5], flag, Util::Mask1, 1U);
                }
            }
        } else {
            yyjson_val *face = faces;
            if (yyjson_equals_str(face, "x-")) {
                var.obflag[0] = Util::ibset(var.obflag[0], flag, Util::Mask1, 1U);
            } else if (yyjson_equals_str(face, "x+")) {
                var.obflag[1] = Util::ibset(var.obflag[1], flag, Util::Mask1, 1U);
            } else if (yyjson_equals_str(face, "y-")) {
                var.obflag[2] = Util::ibset(var.obflag[2], flag, Util::Mask1, 1U);
            } else if (yyjson_equals_str(face, "y+")) {
                var.obflag[3] = Util::ibset(var.obflag[3], flag, Util::Mask1, 1U);
            } else if (yyjson_equals_str(face, "z-")) {
                var.obflag[4] = Util::ibset(var.obflag[4], flag, Util::Mask1, 1U);
            } else if (yyjson_equals_str(face, "z+")) {
                var.obflag[5] = Util::ibset(var.obflag[5], flag, Util::Mask1, 1U);
            }
        }
    }
}

void load_bcfi(scalar_field<real_t> &var, unsigned flag, yyjson_val *topo) {
    yyjson_val *io = yyjson_obj_get(topo, "io");
    if (yyjson_equals_str(io, "outer")) {
        yyjson_val *faces = yyjson_obj_get(topo, "face");
        if (yyjson_is_arr(faces)) {
            yyjson_arr_iter iter;
            yyjson_arr_iter_init(faces, &iter);
            yyjson_val *face;
            while (face = yyjson_arr_iter_next(&iter)) {
                if (yyjson_equals_str(face, "x-")) {
                    var.obflag[0] = Util::ibset(var.obflag[0], flag, Util::Mask1, 1U);
                } else if (yyjson_equals_str(face, "x+")) {
                    var.obflag[1] = Util::ibset(var.obflag[1], flag, Util::Mask1, 1U);
                } else if (yyjson_equals_str(face, "y-")) {
                    var.obflag[2] = Util::ibset(var.obflag[2], flag, Util::Mask1, 1U);
                } else if (yyjson_equals_str(face, "y+")) {
                    var.obflag[3] = Util::ibset(var.obflag[3], flag, Util::Mask1, 1U);
                } else if (yyjson_equals_str(face, "z-")) {
                    var.obflag[4] = Util::ibset(var.obflag[4], flag, Util::Mask1, 1U);
                } else if (yyjson_equals_str(face, "z+")) {
                    var.obflag[5] = Util::ibset(var.obflag[5], flag, Util::Mask1, 1U);
                }
            }
        } else {
            yyjson_val *face = faces;
            if (yyjson_equals_str(face, "x-")) {
                var.obflag[0] = Util::ibset(var.obflag[0], flag, Util::Mask1, 1U);
            } else if (yyjson_equals_str(face, "x+")) {
                var.obflag[1] = Util::ibset(var.obflag[1], flag, Util::Mask1, 1U);
            } else if (yyjson_equals_str(face, "y-")) {
                var.obflag[2] = Util::ibset(var.obflag[2], flag, Util::Mask1, 1U);
            } else if (yyjson_equals_str(face, "y+")) {
                var.obflag[3] = Util::ibset(var.obflag[3], flag, Util::Mask1, 1U);
            } else if (yyjson_equals_str(face, "z-")) {
                var.obflag[4] = Util::ibset(var.obflag[4], flag, Util::Mask1, 1U);
            } else if (yyjson_equals_str(face, "z+")) {
                var.obflag[5] = Util::ibset(var.obflag[5], flag, Util::Mask1, 1U);
            }
        }
    }
}

int search_boundary(yyjson_val *boundary, const char *label) {
    yyjson_arr_iter iter;
    yyjson_arr_iter_init(boundary, &iter);
    int index = 1;
    yyjson_val *bound;
    while (bound = yyjson_arr_iter_next(&iter)) {
        if (yyjson_equals_str(yyjson_obj_get(bound, "label"), label)) {
            return index;
        }
        index ++;
    }
    return 0;
}

JReader::JReader(const char *fname) {
    doc = yyjson_read_file(fname, 0, NULL, NULL);
    root = yyjson_doc_get_root(doc);
}

JReader::~JReader() {
    yyjson_doc_free(doc);
}

void JReader::load_domain(Dom &dom, bool load_mesh) {
    yyjson_val *domain = yyjson_obj_get(root, dom.name);
    int size[4];
    size[0] = yyjson_get_int(yyjson_arr_get(domain, 0)) + 2 * GUIDE;
    size[1] = yyjson_get_int(yyjson_arr_get(domain, 1)) + 2 * GUIDE;
    size[2] = yyjson_get_int(yyjson_arr_get(domain, 2)) + 2 * GUIDE;
    size[3] = 3;
    yyjson_val *boundary = yyjson_obj_get(root, "boundary");
    int bcount = yyjson_arr_size(boundary);
    dom.init(size, bcount + 1);

    Ctrl &c                   = dom.c;
    vector_field<real_t> &u   = dom.u;
    scalar_field<real_t> &p   = dom.p;
    scalar_field<real_t> &nut = dom.nut;

    yyjson_arr_iter iter;
    yyjson_arr_iter_init(boundary, &iter);
    int index = 1;
    while (yyjson_val *bound = yyjson_arr_iter_next(&iter)) {
        yyjson_val *var, *type, *value;
        var = yyjson_obj_get(bound, "u");
        if (var) {
            
            unsigned int flag = 0U;
            type = yyjson_obj_get(var, "type");
            if (yyjson_equals_str(type, "fixedValue")) {
                flag = Util::ibset(flag, BB::dirichlet, Util::Mask1, 1U);
            } else if (yyjson_equals_str(type, "fixedGradiet")) {
                flag = Util::ibset(flag, BB::neumann, Util::Mask1, 1U);
            } else if (yyjson_equals_str(type, "symmetric")) {
                load_bcfi(u, BB::symmetric, yyjson_obj_get(bound, "topo"));
            } else if (yyjson_equals_str(type, "cyclic")) {
                load_bcfi(u, BB::cyclic, yyjson_obj_get(bound, "topo"));
            } else if (yyjson_equals_str(type, "directional")) {
                load_bcfi(u, BB::cyclic, yyjson_obj_get(bound, "topo"));
            } else if (yyjson_equals_str(type, "outflow")) {
                load_bcfi(u, BB::outflow, yyjson_obj_get(bound, "topo"));
            } else if (yyjson_equals_str(type, "driver")) {
                load_bcfi(u, BB::semicyclic, yyjson_obj_get(bound, "topo"));
            }
            value = yyjson_obj_get(var, "value");
            u.b[id2(index,0,u.bsize)] = yyjson_get_real(yyjson_arr_get(value, 0));
            u.b[id2(index,1,u.bsize)] = yyjson_get_real(yyjson_arr_get(value, 1));
            u.b[id2(index,2,u.bsize)] = yyjson_get_real(yyjson_arr_get(value, 2));
            u.bflag[index] = flag;
        }
        var = yyjson_obj_get(bound, "p");
        if (var) {
            
            unsigned int flag = 0U;
            type = yyjson_obj_get(var, "type");
            if (yyjson_equals_str(type, "fixedValue")) {
                flag = Util::ibset(flag, BB::dirichlet, Util::Mask1, 1U);
            } else if (yyjson_equals_str(type, "fixedGradient")) {
                flag = Util::ibset(flag, BB::neumann, Util::Mask1, 1U);
            } else if (yyjson_equals_str(type, "symmetric")) {
                load_bcfi(p, BB::symmetric, yyjson_obj_get(bound, "topo"));
            } else if (yyjson_equals_str(type, "cyclic")) {
                load_bcfi(p, BB::cyclic, yyjson_obj_get(bound, "topo"));
            } else if (yyjson_equals_str(type, "directional")) {
                load_bcfi(p, BB::directional, yyjson_obj_get(bound, "topo"));
            } else if (yyjson_equals_str(type, "driver")) {
                load_bcfi(p, BB::driver, yyjson_obj_get(bound, "topo"));
            }
            value = yyjson_obj_get(var, "value");
            p.b[index] = yyjson_get_real(value);
            p.bflag[index] = flag;
        }
        var = yyjson_obj_get(bound, "nut");
        if (var) {
            
            unsigned int flag = 0U;
            type = yyjson_obj_get(var, "type");
            if (yyjson_equals_str(type, "fixedValue")) {
                flag = Util::ibset(flag, BB::dirichlet, Util::Mask1, 1U);
            } else if (yyjson_equals_str(type, "fixedGradient")) {
                flag = Util::ibset(flag, BB::neumann, Util::Mask1, 1U);
            } else if (yyjson_equals_str(type, "symmetric")) {
                load_bcfi(nut, BB::symmetric, yyjson_obj_get(bound, "topo"));
            } else if (yyjson_equals_str(type, "cyclic")) {
                load_bcfi(nut, BB::cyclic, yyjson_obj_get(bound, "topo"));
            } else if (yyjson_equals_str(type, "directional")) {
                load_bcfi(nut, BB::cyclic, yyjson_obj_get(bound, "topo"));
            } else if (yyjson_equals_str(type, "driver")) {
                load_bcfi(nut, BB::semicyclic, yyjson_obj_get(bound, "topo"));
            }
            value = yyjson_obj_get(var, "value");
            nut.b[index] = yyjson_get_real(value);
            nut.bflag[index] = flag;
        }
        index ++;
    }

    yyjson_val *driver = yyjson_obj_get(root, "driver");
    const char *directions[3] = {"x", "y", "z"};
    const char *outer_face[6] = {"x-", "x+", "y-", "y+", "z-", "z+"};
    yyjson_val *direction;
    for (int i = 0; i < 3; i ++) {
        direction = yyjson_obj_get(driver, directions[i]);
        if (!direction) {
            continue;
        }
        if (yyjson_equals_str(yyjson_obj_get(direction, "type"), "directional")) {
            c.driver[i].type = Ctrl::Driver::Type::directional;
            if (yyjson_val *var = yyjson_obj_get(direction, "u")) {
                c.driver[i].var = Ctrl::Driver::Var::U;
                c.driver[i].value = yyjson_get_real(var);
            } else if (yyjson_val *var = yyjson_obj_get(direction, "p")) {
                c.driver[i].var = Ctrl::Driver::Var::P;
                c.driver[i].value = yyjson_get_real(var);
            }
        } else if (yyjson_equals_str(yyjson_obj_get(direction, "type"), "driver")) {
            c.driver[i].type = Ctrl::Driver::Type::driver;
            if (yyjson_val *var = yyjson_obj_get(direction, "u")) {
                c.driver[i].var = Ctrl::Driver::Var::U;
                c.driver[i].value = yyjson_get_real(var);
            } else if (yyjson_val *var = yyjson_obj_get(direction, "p")) {
                c.driver[i].var = Ctrl::Driver::Var::P;
                c.driver[i].value = yyjson_get_real(var);
            }
            c.driver[i].length = yyjson_get_int(yyjson_obj_get(direction, "length"));
            if (yyjson_equals_str(yyjson_obj_get(direction, "inflow"), outer_face[i*2 + 0])) {
                c.driver[i].inflow = i*2 + 0;
            } else if (yyjson_equals_str(yyjson_obj_get(direction, "inflow"), outer_face[i*2 + 1])) {
                c.driver[i].inflow = i*2 + 1;
            }
        }
    }

    scalar_field<unsigned> &f = dom.f;
    for (int i = GUIDE; i < dom.size[0] - GUIDE; i ++) {
        for (int j = GUIDE; j < dom.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < dom.size[2] - GUIDE; k ++) {
                f.m[id3(i,j,k,f.size)] = Util::ibset(f.m[id3(i,j,k,f.size)], Flag::Active, Util::Mask1, 1);
            }
        }
    }
    yyjson_arr_iter_init(boundary, &iter);
    index = 1;
    while (yyjson_val *bound = yyjson_arr_iter_next(&iter)) {
        yyjson_val *topo = yyjson_obj_get(bound, "topo");
        yyjson_val *io = yyjson_obj_get(topo, "io");
        if (yyjson_equals_str(io, "outer")) {
            yyjson_val *faces = yyjson_obj_get(topo, "face");
            yyjson_arr_iter fiter;
            yyjson_arr_iter_init(faces, &fiter);
            yyjson_val *face;
            while (face = yyjson_arr_iter_next(&fiter)) {
                if (yyjson_equals_str(face, "x+")) {
                    for (int j = GUIDE; j < dom.size[1] - GUIDE; j ++) {
                        for (int k = GUIDE; k < dom.size[2] - GUIDE; k ++) {
                            int i = dom.size[0] - GUIDE - 1;
                            unsigned int flag = f.m[id3(i,j,k,f.size)];
                            flag = Util::ibset(flag, Flag::Fe, Util::Mask8, index);
                            flag = Util::ibset(flag, Flag::Me, Util::Mask1, 0);
                            f.m[id3(i,j,k,f.size)] = flag;
                        }
                    }
                } else if (yyjson_equals_str(face, "x-")) {
                    for (int j = GUIDE; j < dom.size[1] - GUIDE; j ++) {
                        for (int k = GUIDE; k < dom.size[2] - GUIDE; k ++) {
                            int i = GUIDE - 1;
                            unsigned int flag = f.m[id3(i,j,k,f.size)];
                            flag = Util::ibset(flag, Flag::Fe, Util::Mask8, index);
                            flag = Util::ibset(flag, Flag::Me, Util::Mask1, 1);
                            f.m[id3(i,j,k,f.size)] = flag;
                        }
                    }
                } else if (yyjson_equals_str(face, "y+")) {
                    for (int i = GUIDE; i < dom.size[0] - GUIDE; i ++) {
                        for (int k = GUIDE; k < dom.size[2] - GUIDE; k ++) {
                            int j = dom.size[1] - GUIDE - 1;
                            unsigned int flag = f.m[id3(i,j,k,f.size)];
                            flag = Util::ibset(flag, Flag::Fn, Util::Mask8, index);
                            flag = Util::ibset(flag, Flag::Mn, Util::Mask1, 0);
                            f.m[id3(i,j,k,f.size)] = flag;
                        }
                    }
                } else if (yyjson_equals_str(face, "y-")) {
                    for (int i = GUIDE; i < dom.size[0] - GUIDE; i ++) {
                        for (int k = GUIDE; k < dom.size[2] - GUIDE; k ++) {
                            int j = GUIDE - 1;
                            unsigned int flag = f.m[id3(i,j,k,f.size)];
                            flag = Util::ibset(flag, Flag::Fn, Util::Mask8, index);
                            flag = Util::ibset(flag, Flag::Mn, Util::Mask1, 1);
                            f.m[id3(i,j,k,f.size)] = flag;
                        }
                    }
                } else if (yyjson_equals_str(face, "z+")) {
                    for (int i = GUIDE; i < dom.size[0] - GUIDE; i ++) {
                        for (int j = GUIDE; j < dom.size[1] - GUIDE; j ++) {
                            int k = dom.size[2] - GUIDE - 1;
                            unsigned int flag = f.m[id3(i,j,k,f.size)];
                            flag = Util::ibset(flag, Flag::Ft, Util::Mask8, index);
                            flag = Util::ibset(flag, Flag::Mt, Util::Mask1, 0);
                            f.m[id3(i,j,k,f.size)] = flag;
                        }
                    }
                } else if (yyjson_equals_str(face, "z-")) {
                    for (int i = GUIDE; i < dom.size[0] - GUIDE; i ++) {
                        for (int j = GUIDE; j < dom.size[1] - GUIDE; j ++) {
                            int k = GUIDE - 1;
                            unsigned int flag = f.m[id3(i,j,k,f.size)];
                            flag = Util::ibset(flag, Flag::Ft, Util::Mask8, index);
                            flag = Util::ibset(flag, Flag::Mt, Util::Mask1, 1);
                            f.m[id3(i,j,k,f.size)] = flag;
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
                for (int i = p1[0] + GUIDE; i <= p2[0] + GUIDE; i ++) {
                    for (int j = p1[1] + GUIDE; j <= p2[1] + GUIDE; j ++) {
                        for (int k = p1[2] + GUIDE; k <= p2[2] + GUIDE; k ++) {
                            f.m[id3(i,j,k,f.size)] = Util::ibset(f.m[id3(i,j,k,f.size)], Flag::Active, Util::Mask1, 0);
                        }
                    }
                }
                yyjson_val *faces = yyjson_obj_get(bound, "face");
                yyjson_val *f0, *f1, *f2, *f3, *f4, *f5;
                int i0, i1, i2, i3, i4, i5;
                f0 = yyjson_obj_get(faces, "x-");
                f1 = yyjson_obj_get(faces, "x+");
                f2 = yyjson_obj_get(faces, "y-");
                f3 = yyjson_obj_get(faces, "y+");
                f4 = yyjson_obj_get(faces, "z-");
                f5 = yyjson_obj_get(faces, "z+");
                i0 = search_boundary(boundary, yyjson_get_str(f0));
                i1 = search_boundary(boundary, yyjson_get_str(f1));
                i2 = search_boundary(boundary, yyjson_get_str(f2));
                i3 = search_boundary(boundary, yyjson_get_str(f3));
                i4 = search_boundary(boundary, yyjson_get_str(f4));
                i5 = search_boundary(boundary, yyjson_get_str(f5));
                for (int j = p1[1]  + GUIDE; j <= p2[1] + GUIDE; j ++) {
                    for (int k = p1[2] + GUIDE; k <= p2[2] + GUIDE; k ++) {
                        int i = p1[0] - 1 + GUIDE;
                        unsigned int flag = f.m[id3(i,j,k,f.size)];
                        if (f0) {
                            flag = Util::ibset(flag, Flag::Fe, Util::Mask8, i0);
                        } else {
                            flag = Util::ibset(flag, Flag::Fe, Util::Mask8, index);
                        }
                        flag = Util::ibset(flag, Flag::Me, Util::Mask1, 0);
                        f.m[id3(i,j,k,f.size)] = flag;

                        i = p2[0] + GUIDE;
                        flag = f.m[id3(i,j,k,f.size)];
                        if (f1) {
                            flag = Util::ibset(flag, Flag::Fe, Util::Mask8, i1);
                        } else {
                            flag = Util::ibset(flag, Flag::Fe, Util::Mask8, index);
                        }
                        flag = Util::ibset(flag, Flag::Me, Util::Mask1, 1);
                        f.m[id3(i,j,k,f.size)] = flag;
                    }
                }
                for (int i = p1[0] + GUIDE; i <= p2[0] + GUIDE; i ++) {
                    for (int k = p1[2] + GUIDE; k <= p2[2] + GUIDE; k ++) {
                        int j = p1[1] - 1 + GUIDE;
                        unsigned int flag = f.m[id3(i,j,k,f.size)];
                        if (f2) {
                            flag = Util::ibset(flag, Flag::Fn, Util::Mask8, i2);
                        } else {
                            flag = Util::ibset(flag, Flag::Fn, Util::Mask8, index);
                        }
                        flag = Util::ibset(flag, Flag::Mn, Util::Mask1, 0);
                        f.m[id3(i,j,k,f.size)] = flag;

                        j = p2[1] + GUIDE;
                        flag = f.m[id3(i,j,k,f.size)];
                        if (f3) {
                            flag = Util::ibset(flag, Flag::Fn, Util::Mask8, i3);
                        } else {
                            flag = Util::ibset(flag, Flag::Fn, Util::Mask8, index);
                        }
                        flag = Util::ibset(flag, Flag::Mn, Util::Mask1, 1);
                        f.m[id3(i,j,k,f.size)] = flag;
                    }
                }
                for (int i = p1[0] + GUIDE; i <= p2[0] + GUIDE; i ++) {
                    for (int j = p1[1] + GUIDE; j <= p2[1] + GUIDE; j ++) {
                        int k = p1[2] - 1 + GUIDE;
                        unsigned int flag = f.m[id3(i,j,k,f.size)];
                        if (f4) {
                            flag = Util::ibset(flag, Flag::Ft, Util::Mask8, i4);
                        } else {
                            flag = Util::ibset(flag, Flag::Ft, Util::Mask8, index);
                        }
                        flag = Util::ibset(flag, Flag::Mt, Util::Mask1, 0);
                        f.m[id3(i,j,k,f.size)] = flag;

                        k = p2[2] + GUIDE;
                        flag = f.m[id3(i,j,k,f.size)];
                        if (f5) {
                            flag = Util::ibset(flag, Flag::Ft, Util::Mask8, i5);
                        } else {
                            flag = Util::ibset(flag, Flag::Ft, Util::Mask8, index);
                        }
                        flag = Util::ibset(flag, Flag::Mt, Util::Mask1, 1);
                        f.m[id3(i,j,k,f.size)] = flag;
                    }
                }
            }
        }
        index ++;
    }

    yyjson_val *init = yyjson_obj_get(domain, "init");
    real_t u_init = yyjson_get_real(yyjson_arr_get(yyjson_obj_get(init, "u"), 0));
    real_t v_init = yyjson_get_real(yyjson_arr_get(yyjson_obj_get(init, "u"), 1));
    real_t w_init = yyjson_get_real(yyjson_arr_get(yyjson_obj_get(init, "u"), 2));
    real_t p_init = yyjson_get_real(yyjson_obj_get(init, "p"));
    for (int i = 0; i < dom.size[0]; i ++) {
        for (int j = 0; j < dom.size[1]; j ++) {
            for (int k = 0; k < dom.size[2]; k ++) {
                if (Util::ibsee(f.m[id3(i,j,k,f.size)], Flag::Active, Util::Mask1)) {
                    u.m[id4(i,j,k,0,u.size)] = u_init;
                    u.m[id4(i,j,k,1,u.size)] = v_init;
                    u.m[id4(i,j,k,2,u.size)] = w_init;
                    p.m[id3(i,j,k,  p.size)] = p_init;
                }
            }
        }
    }

    yyjson_val *time = yyjson_obj_get(root, "time");
    c.time.dt = yyjson_get_real(yyjson_obj_get(time, "dt"));
    if (yyjson_obj_get(time, "t")) {
        c.time.ndt = int(yyjson_get_real(yyjson_obj_get(time, "t")) / c.time.dt);
    } else if (yyjson_obj_get(time, "ndt")) {
        c.time.ndt = yyjson_get_int(yyjson_obj_get(time, "ndt"));
    }

    yyjson_val *flow = yyjson_obj_get(root, "flow");
    c.flow.re = yyjson_get_int(yyjson_obj_get(flow, "re"));
    c.flow.ri = 1 / c.flow.re;
    if (yyjson_equals_str(yyjson_obj_get(flow, "scheme"), "upwind3")) {
        c.flow.scheme = Ctrl::Flow::Scheme::upwind3;
    } else if (yyjson_equals_str(yyjson_obj_get(flow, "scheme"), "muscl")) {
        c.flow.scheme = Ctrl::Flow::Scheme::muscl;
    }
    c.flow.tdiv = yyjson_get_real(yyjson_obj_get(flow, "tdiv"));
    c.flow.maxit = yyjson_get_int(yyjson_obj_get(flow, "maxit"));

    c.monitor.type = Ctrl::Monitor::Type::off;
    yyjson_val *monitor = yyjson_obj_get(root, "monitor");
    if (monitor) {
        if (yyjson_equals_str(yyjson_obj_get(monitor, "type"), "t")) {
            c.monitor.type = Ctrl::Monitor::Type::T;
            c.monitor.t_interval = yyjson_get_real(yyjson_obj_get(monitor, "interval"));
        } else if (yyjson_equals_str(yyjson_obj_get(monitor, "type"), "dt")) {
            c.monitor.type = Ctrl::Monitor::Type::Dt;
            c.monitor.dt_interval = yyjson_get_int(yyjson_obj_get(monitor, "interval"));
        }
    }

    yyjson_val *poisson = yyjson_obj_get(root, "poisson");
    if (yyjson_equals_str(yyjson_obj_get(poisson, "solver"), "jacobi")) {
        c.poisson.type = Ctrl::LS::Type::jacobi;
    } else if (yyjson_equals_str(yyjson_obj_get(poisson, "solver"), "sor")) {
        c.poisson.type = Ctrl::LS::Type::sor;
    } else if (yyjson_equals_str(yyjson_obj_get(poisson, "solver"), "bicgstab")) {
        c.poisson.type = Ctrl::LS::Type::bicgstab;
    } else if (yyjson_equals_str(yyjson_obj_get(poisson, "solver"), "pbicgstab")) {
        c.poisson.type = Ctrl::LS::Type::pbicgstab;
    }
    c.poisson.omega = yyjson_get_real(yyjson_obj_get(poisson, "omega"));
    c.poisson.epsilon = yyjson_get_real(yyjson_obj_get(poisson, "e"));
    c.poisson.maxit = yyjson_get_int(yyjson_obj_get(poisson, "maxit"));

    if (!load_mesh) {
        return;
    }

    char mesh_file[128];
    sprintf(mesh_file, "./mesh/%s.json", dom.name);
    vector_field<real_t> &x  = dom.x;
    vector_field<real_t> &kx = dom.kx;
    vector_field<real_t> &g  = dom.g;
    scalar_field<real_t> &ja = dom.ja;
    yyjson_doc *mesh_doc = yyjson_read_file(mesh_file, 0, NULL, NULL);
    yyjson_val *mesh_root = yyjson_doc_get_root(mesh_doc);
    yyjson_val *points = yyjson_obj_get(mesh_root, "point");
    yyjson_arr_iter_init(points, &iter);
    for (int i = 0; i < dom.size[0]; i ++) {
        for (int j = 0; j < dom.size[1]; j ++) {
            for (int k = 0; k < dom.size[2]; k ++) {
                yyjson_val *point = yyjson_arr_iter_next(&iter);
                x.m[id4(i,j,k,0,x.size)] = yyjson_get_real(yyjson_arr_get(point, 0));
                x.m[id4(i,j,k,1,x.size)] = yyjson_get_real(yyjson_arr_get(point, 1));
                x.m[id4(i,j,k,2,x.size)] = yyjson_get_real(yyjson_arr_get(point, 2));
            }
        }
    }
    for (int i = GUIDE - 1; i <= dom.size[0] - GUIDE; i ++) {
        for (int j = GUIDE - 1; j <= dom.size[1] - GUIDE; j ++) {
            for (int k = GUIDE - 1; k <= dom.size[2] - GUIDE; k ++) {
                real_t xe, xw, yn, ys, zt, zb;
                xe = x.m[id4(i+1,j  ,k  ,0,x.size)];
                xw = x.m[id4(i-1,j  ,k  ,0,x.size)];
                yn = x.m[id4(i  ,j+1,k  ,1,x.size)];
                ys = x.m[id4(i  ,j-1,k  ,1,x.size)];
                zt = x.m[id4(i  ,j  ,k+1,2,x.size)];
                zb = x.m[id4(i  ,j  ,k-1,2,x.size)];

                real_t x1, x2, x3, k1, k2, k3, g1, g2, g3, de;
                x1 = 0.5 * (xe - xw);
                x2 = 0.5 * (yn - ys);
                x3 = 0.5 * (zt - zb);
                k1 = 1 / x1;
                k2 = 1 / x2;
                k3 = 1 / x3;
                de = x1 * x2 * x3;
                g1 = de * k1 * k1;
                g2 = de * k2 * k2;
                g3 = de * k3 * k3;
                
                kx.m[id4(i,j,k,0,kx.size)] = k1;
                kx.m[id4(i,j,k,1,kx.size)] = k2;
                kx.m[id4(i,j,k,2,kx.size)] = k3;
                g.m[ id4(i,j,k,0, g.size)] = g1;
                g.m[ id4(i,j,k,1, g.size)] = g2;
                g.m[ id4(i,j,k,2, g.size)] = g3;
                ja.m[id3(i,j,k,  ja.size)] = de;
            }
        }
    }
}

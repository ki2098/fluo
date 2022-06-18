#include <stdio.h>
#include "jreader.h"
#include "util.h"
#include "alloc.h"
#include "flag.h"

void load_bcfi(BC *bc, unsigned int flag, yyjson_val *topo) {
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

JReader::JReader(char *fname) {
    this->doc = yyjson_read_file(fname, 0, NULL, NULL);
    this->root = yyjson_doc_get_root(this->doc);
}

JReader::~JReader() {
    yyjson_doc_free(this->doc);
}

void JReader::load_bc(BC* bc, Driver* driver) {
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
            } else if (yyjson_equals_str(type, "outflow")) {
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
        var = yyjson_obj_get(bound, "Nut");
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
            bc->b[index].nut = yyjson_get_real(value);
        }
        bc->b[index].flag = flag;
        index ++;
    }
}

void JReader::load_domain(D *dom, LS *ls) {
    yyjson_val *domain = yyjson_obj_get(this->root, "domain");
    dom->size[0] = yyjson_get_int(yyjson_arr_get(domain, 0)) + 2 * D::GUIDE;
    dom->size[1] = yyjson_get_int(yyjson_arr_get(domain, 1)) + 2 * D::GUIDE;
    dom->size[2] = yyjson_get_int(yyjson_arr_get(domain, 2)) + 2 * D::GUIDE;
    int size[4] = {dom->size[0], dom->size[1], dom->size[2], 3};

    dom->f   = Alloc::uint_3d(size);
    dom->u   = Alloc::double_4d(size);
    dom->uu  = Alloc::double_4d(size);
    dom->uc  = Alloc::double_4d(size);
    dom->ua  = Alloc::double_4d(size);
    dom->uua = Alloc::double_4d(size);
    dom->ur  = Alloc::double_4d(size);
    dom->p   = Alloc::double_3d(size);
    dom->pd  = Alloc::double_3d(size);
    dom->pr  = Alloc::double_3d(size);
    dom->sgs = Alloc::double_3d(size);
    dom->x   = Alloc::double_4d(size);
    dom->kx  = Alloc::double_4d(size);
    dom->g   = Alloc::double_4d(size);
    dom->ja  = Alloc::double_3d(size);
    
    unsigned int *f = dom->f;

    for (int i = D::GUIDE; i < size[0] - D::GUIDE; i ++) {
        for (int j = D::GUIDE; j < size[1] - D::GUIDE; j ++) {
            for (int k = D::GUIDE; k < size[2] - D::GUIDE; k ++) {
                f[id3(i,j,k,size)] = Util::ibset(f[id3(i,j,k,size)], Flag::Active, Util::Mask1, 1);
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
            yyjson_arr_iter fiter;
            yyjson_arr_iter_init(faces, &fiter);
            yyjson_val *face;
            while (face = yyjson_arr_iter_next(&fiter)) {
                if (yyjson_equals_str(face, "x+")) {
                    for (int j = D::GUIDE; j < size[1] - D::GUIDE; j ++) {
                        for (int k = D::GUIDE; k < size[2] - D::GUIDE; k ++) {
                            int i = size[0] - D::GUIDE - 1;
                            unsigned int flag = f[id3(i,j,k,size)];
                            flag = Util::ibset(flag, Flag::Fe, Util::Mask8, index);
                            flag = Util::ibset(flag, Flag::Me, Util::Mask1, 0);
                            f[id3(i,j,k,size)] = flag;
                        }
                    }
                } else if (yyjson_equals_str(face, "x-")) {
                    for (int j = D::GUIDE; j <size[1] - D::GUIDE; j ++) {
                        for (int k = D::GUIDE; k < size[2] - D::GUIDE; k ++) {
                            int i = D::GUIDE - 1;
                            unsigned int flag = f[id3(i,j,k,size)];
                            flag = Util::ibset(flag, Flag::Fe, Util::Mask8, index);
                            flag = Util::ibset(flag, Flag::Me, Util::Mask1, 1);
                            f[id3(i,j,k,size)] = flag;
                        }
                    }
                } else if (yyjson_equals_str(face, "y+")) {
                    for (int i = D::GUIDE; i < size[0] - D::GUIDE; i ++) {
                        for (int k = D::GUIDE; k < size[2] - D::GUIDE; k ++) {
                            int j = size[1] - D::GUIDE - 1;
                            unsigned int flag = f[id3(i,j,k,size)];
                            flag = Util::ibset(flag, Flag::Fn, Util::Mask8, index);
                            flag = Util::ibset(flag, Flag::Mn, Util::Mask1, 0);
                            f[id3(i,j,k,size)] = flag;
                        }
                    }
                } else if (yyjson_equals_str(face, "y-")) {
                    for (int i = D::GUIDE; i < size[0] - D::GUIDE; i ++) {
                        for (int k = D::GUIDE; k < size[2] - D::GUIDE; k ++) {
                            int j = D::GUIDE - 1;
                            unsigned int flag = f[id3(i,j,k,size)];
                            flag = Util::ibset(flag, Flag::Fn, Util::Mask8, index);
                            flag = Util::ibset(flag, Flag::Mn, Util::Mask1, 1);
                            f[id3(i,j,k,size)] = flag;
                        }
                    }
                } else if (yyjson_equals_str(face, "z+")) {
                    for (int i = D::GUIDE; i < size[0] - D::GUIDE; i ++) {
                        for (int j = D::GUIDE; j < size[1] - D::GUIDE; j ++) {
                            int k = size[2] - D::GUIDE - 1;
                            unsigned int flag = f[id3(i,j,k,size)];
                            flag = Util::ibset(flag, Flag::Ft, Util::Mask8, index);
                            flag = Util::ibset(flag, Flag::Mt, Util::Mask1, 0);
                            f[id3(i,j,k,size)] = flag;
                        }
                    }
                } else if (yyjson_equals_str(face, "z-")) {
                    for (int i = D::GUIDE; i < size[0] - D::GUIDE; i ++) {
                        for (int j = D::GUIDE; j < size[1] - D::GUIDE; j ++) {
                            int k = D::GUIDE - 1;
                            unsigned int flag = f[id3(i,j,k,size)];
                            flag = Util::ibset(flag, Flag::Ft, Util::Mask8, index);
                            flag = Util::ibset(flag, Flag::Mt, Util::Mask1, 1);
                            f[id3(i,j,k,size)] = flag;
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
                            f[id3(i,j,k,size)] = Util::ibset(f[id3(i,j,k,size)], Flag::Active, Util::Mask1, 0);
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
                for (int j = p1[1]  + D::GUIDE; j <= p2[1] + D::GUIDE; j ++) {
                    for (int k = p1[2] + D::GUIDE; k <= p2[2] + D::GUIDE; k ++) {
                        int i = p1[0] - 1 + D::GUIDE;
                        unsigned int flag = f[id3(i,j,k,size)];
                        if (f0) {
                            flag = Util::ibset(flag, Flag::Fe, Util::Mask8, i0);
                        } else {
                            flag = Util::ibset(flag, Flag::Fe, Util::Mask8, index);
                        }
                        flag = Util::ibset(flag, Flag::Me, Util::Mask1, 0);
                        f[id3(i,j,k,size)] = flag;

                        i = p2[0] + D::GUIDE;
                        flag = f[id3(i,j,k,size)];
                        if (f1) {
                            flag = Util::ibset(flag, Flag::Fe, Util::Mask8, i1);
                        } else {
                            flag = Util::ibset(flag, Flag::Fe, Util::Mask8, index);
                        }
                        flag = Util::ibset(flag, Flag::Me, Util::Mask1, 1);
                        f[id3(i,j,k,size)] = flag;
                    }
                }
                for (int i = p1[0] + D::GUIDE; i <= p2[0] + D::GUIDE; i ++) {
                    for (int k = p1[2] + D::GUIDE; k <= p2[2] + D::GUIDE; k ++) {
                        int j = p1[1] - 1 + D::GUIDE;
                        unsigned int flag = f[id3(i,j,k,size)];
                        if (f2) {
                            flag = Util::ibset(flag, Flag::Fn, Util::Mask8, i2);
                        } else {
                            flag = Util::ibset(flag, Flag::Fn, Util::Mask8, index);
                        }
                        flag = Util::ibset(flag, Flag::Mn, Util::Mask1, 0);
                        f[id3(i,j,k,size)] = flag;

                        j = p2[1] + D::GUIDE;
                        flag = f[id3(i,j,k,size)];
                        if (f3) {
                            flag = Util::ibset(flag, Flag::Fn, Util::Mask8, i3);
                        } else {
                            flag = Util::ibset(flag, Flag::Fn, Util::Mask8, index);
                        }
                        flag = Util::ibset(flag, Flag::Mn, Util::Mask1, 1);
                        f[id3(i,j,k,size)] = flag;
                    }
                }
                for (int i = p1[0] + D::GUIDE; i <= p2[0] + D::GUIDE; i ++) {
                    for (int j = p1[1] + D::GUIDE; j <= p2[1] + D::GUIDE; j ++) {
                        int k = p1[2] - 1 + D::GUIDE;
                        unsigned int flag = f[id3(i,j,k,size)];
                        if (f4) {
                            flag = Util::ibset(flag, Flag::Ft, Util::Mask8, i4);
                        } else {
                            flag = Util::ibset(flag, Flag::Ft, Util::Mask8, index);
                        }
                        flag = Util::ibset(flag, Flag::Mt, Util::Mask1, 0);
                        f[id3(i,j,k,size)] = flag;

                        k = p2[2] + D::GUIDE;
                        flag = f[id3(i,j,k,size)];
                        if (f5) {
                            flag = Util::ibset(flag, Flag::Ft, Util::Mask8, i5);
                        } else {
                            flag = Util::ibset(flag, Flag::Ft, Util::Mask8, index);
                        }
                        flag = Util::ibset(flag, Flag::Mt, Util::Mask1, 1);
                        f[id3(i,j,k,size)] = flag;
                    }
                }
            }
        }
        index ++;
    }
    
    yyjson_val *init = yyjson_obj_get(this->root, "init");
    double u_init = yyjson_get_real(yyjson_arr_get(yyjson_obj_get(init, "U"), 0));
    double v_init = yyjson_get_real(yyjson_arr_get(yyjson_obj_get(init, "U"), 1));
    double w_init = yyjson_get_real(yyjson_arr_get(yyjson_obj_get(init, "U"), 2));
    double p_init = yyjson_get_real(yyjson_obj_get(init, "P"));
    for (int i = 0; i < size[0]; i ++) {
        for (int j = 0; j < size[1]; j ++) {
            for (int k = 0; k < size[2]; k ++) {
                if (Util::ibsee(f[id3(i,j,k,size)], Flag::Active, Util::Mask1)) {
                    dom->u[id4(i,j,k,0,size)] = u_init;
                    dom->u[id4(i,j,k,1,size)] = v_init;
                    dom->u[id4(i,j,k,2,size)] = w_init;
                    dom->p[id3(i,j,k,size)]   = p_init;
                }
            }
        }
    }

    yyjson_val *parameter = yyjson_obj_get(this->root, "parameter");
    dom->re = (double)yyjson_get_int(yyjson_obj_get(parameter, "Re"));
    dom->ri = 1.0 / dom->re;
    dom->dt = yyjson_get_real(yyjson_obj_get(parameter, "dt"));
    dom->ntime = yyjson_get_real(yyjson_obj_get(parameter, "timestep"));
    dom->tdiv = yyjson_get_real(yyjson_obj_get(parameter, "TDiv"));

    yyjson_val *monitor = yyjson_obj_get(this->root, "monitor");
    if (monitor) {
        dom->monitor = 1U;
        if (yyjson_equals_str(yyjson_obj_get(monitor, "type"), "time")) {
            double interval = yyjson_get_real(yyjson_obj_get(monitor, "interval"));
            dom->monitor_interval = int(interval / dom->dt);
        } else if (yyjson_equals_str(yyjson_obj_get(monitor, "type"), "timestep")) {
            dom->monitor_interval = yyjson_get_int(yyjson_obj_get(monitor, "interval"));
        }
    }

    yyjson_val *poisson = yyjson_obj_get(this->root, "poisson");
    if (yyjson_equals_str(yyjson_obj_get(poisson, "solver"), "jacobi")) {
        ls->type = LS::Type::jacobi;
    } else if (yyjson_equals_str(yyjson_obj_get(poisson, "solver"), "sor")) {
        ls->type = LS::Type::sor;
    } else if (yyjson_equals_str(yyjson_obj_get(poisson, "solver"), "bicgstab")) {
        ls->type = LS::Type::bicgstab;
    } else if (yyjson_equals_str(yyjson_obj_get(poisson, "solver"), "pbicgstab")) {
        ls->type = LS::Type::pbicgstab;
    }
    ls->omega = yyjson_get_real(yyjson_obj_get(poisson, "omega"));
}

void JReader::load_mesh(D *dom) {
    yyjson_doc *mesh_doc = yyjson_read_file("./setting/mesh.json", 0, NULL, NULL);
    yyjson_val *mesh_root = yyjson_doc_get_root(mesh_doc);
    yyjson_val *points = yyjson_obj_get(mesh_root, "point");
    yyjson_arr_iter iter;
    yyjson_arr_iter_init(points, &iter);
    int size[4] = {dom->size[0], dom->size[1], dom->size[2], 3};
    for (int i = 0; i < size[0]; i ++) {
        for (int j = 0; j < size[1]; j ++) {
            for (int k = 0; k < size[2]; k ++) {
                yyjson_val *point = yyjson_arr_iter_next(&iter);
                // if (point == NULL) {
                //     printf("NO POINT!\n");
                // }
                dom->x[id4(i,j,k,0,size)] = yyjson_get_real(yyjson_arr_get(point, 0));
                dom->x[id4(i,j,k,1,size)] = yyjson_get_real(yyjson_arr_get(point, 1));
                dom->x[id4(i,j,k,2,size)] = yyjson_get_real(yyjson_arr_get(point, 2));
            }
        }
    }
    for (int i = D::GUIDE - 1; i <= size[0] - D::GUIDE; i ++) {
        for (int j = D::GUIDE - 1; j <= size[1] - D::GUIDE; j ++) {
            for (int k = D::GUIDE - 1; k <= size[2] - D::GUIDE; k ++) {
                double xe, xw, yn, ys, zt, zb;
                xe = dom->x[id4(i+1,j  ,k  ,0,size)];
                xw = dom->x[id4(i-1,j  ,k  ,0,size)];
                yn = dom->x[id4(i  ,j+1,k  ,1,size)];
                ys = dom->x[id4(i  ,j-1,k  ,1,size)];
                zt = dom->x[id4(i  ,j  ,k+1,2,size)];
                zb = dom->x[id4(i  ,j  ,k-1,2,size)];

                double x1, x2, x3, k1, k2, k3, g1, g2, g3, de;
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
                
                dom->kx[id4(i,j,k,0,size)] = k1;
                dom->kx[id4(i,j,k,1,size)] = k2;
                dom->kx[id4(i,j,k,2,size)] = k3;
                dom->g[ id4(i,j,k,0,size)] = g1;
                dom->g[ id4(i,j,k,1,size)] = g2;
                dom->g[ id4(i,j,k,2,size)] = g3;
                dom->ja[id3(i,j,k,  size)] = de;
            }
        }
    }
}
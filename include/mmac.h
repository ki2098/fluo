#ifndef _MMAC_H
#define _MMAC_H 1

#include "basic.h"
#include "domain.h"

class MMAC {
public:
    class Poisson {
    public:
        scalar_field<real_t> pd, pcg_r, pcg_r0, pcg_p, pcg_p_, pcg_q, pcg_s, pcg_s_, pcg_t, pcg_t_;
        real_t pcg_rho, pcg_rho_old, pcg_alpha, pcg_beta, pcg_omega;
        int    pcg_it;
    public:
        void sor(scalar_field<real_t> &p, scalar_field<real_t> &psi, Dom &dom);
        void jacobi(scalar_field<real_t> &p, scalar_field<real_t> &psi, Dom &dom);
        void pbicgstab(scalar_field<real_t> &p, scalar_field<real_t> &psi, Dom &dom);
    public:
        void to_device() {
            #pragma acc enter data copyin(this[0:1])
        }
        void end_device() {
            #pragma acc exit data delete(this[0:1])
        }
    private:
        real_t dot(scalar_field<real_t> &v1, scalar_field<real_t> &v2, Dom &dom);
        void triad(scalar_field<real_t> &v1, scalar_field<real_t> &v2, scalar_field<real_t> &v3, real_t a, Dom &dom);
        void calc_ax(scalar_field<real_t> &v1, scalar_field<real_t> &v2, Dom &dom);
        void calc_res(scalar_field<real_t> &p, scalar_field<real_t> &b, scalar_field<real_t> &r, Dom &dom);
        void sor_core(scalar_field<real_t> &p, scalar_field<real_t> &b, Dom &dom, real_t omega, int maxit, real_t torelence, bool ignore, int &it, real_t &res);
        void jacobi_core(scalar_field<real_t> &p, scalar_field<real_t> &psi, Dom &dom, int maxit, real_t torelence, bool ignore, int &it, real_t &res);
    };
public:
    vector_field<real_t> ua, uc, up, uua, uup;
    scalar_field<real_t> pp, psi, diva;
    Poisson poisson;
public:
    void pseudo_velocity(vector_field<real_t> &u, vector_field<real_t> &uu, vector_field<real_t> &ua, scalar_field<real_t> &nut, Dom &dom);
    void correct_center_velocity(vector_field<real_t> &u, vector_field<real_t> &ua, scalar_field<real_t> &p, Dom &dom);
    void correct_face_velocity(vector_field<real_t> &u, vector_field<real_t> &uu, vector_field<real_t> &uua, scalar_field<real_t> &p, Dom &dom);
    void interpolate_velocity(vector_field<real_t> &u, vector_field<real_t> &uc, vector_field<real_t> &uu, Dom &dom);
    void divergence_velocity(vector_field<real_t> &uu, scalar_field<real_t> &div, real_t &div_norm, Dom &dom);
    void calc_rhs(scalar_field<real_t> &div, Dom &dom);
public:
    void to_device() {
        #pragma acc enter data copyin(this[0:1])
    }
    void end_device() {
        #pragma acc exit data delete(this[0:1])
    }
};

#endif

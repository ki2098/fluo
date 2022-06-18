#ifndef _NSEQ_H
#define _NSEQ_H 1

#include "domain.h"
#include "boundary.h"

class NSEq{
public:
    enum Scheme{upwind3, muscl};
    
public:
    unsigned int scheme;
    double alpha;

public:
    NSEq();
    void pseudo_velocity(D *dom, BC *bc);
    void contra_pseudo_velocity(D *dom, BC *bc);
    void correct_ceneter_velocity(D *dom, BC *bc);
    void correct_face_velocity(D *dom, BC *bc);

private:
    static double ffvc_muscl(
        double uc0,
        double ue1,
        double ue2,
        double un1,
        double un2,
        double ut1,
        double ut2,
        double uw1,
        double uw2,
        double us1,
        double us2,
        double ub1,
        double ub2,
        double ufe,
        double vfn,
        double wft,
        double ufw,
        double vfs,
        double wfb,
        double det,
        int    epe,
        int    epn,
        int    ept,
        int    epw,
        int    eps,
        int    epb
    );
    
    static double rmcp_upw3(
        double uc0,
        double ue1,
        double ue2,
        double un1,
        double un2,
        double ut1,
        double ut2,
        double uw1,
        double uw2,
        double us1,
        double us2,
        double ub1,
        double ub2,
        double ufe,
        double vfn,
        double wft,
        double ufw,
        double vfs,
        double wfb,
        double det,
        double alpha
    );

    static double viscosity(
        unsigned int f13,
        unsigned int f23,
        unsigned int f33,
        unsigned int f12,
        unsigned int f22,
        unsigned int f32,
        double       mag,
        double       uc0,
        double       ue1,
        double       un1,
        double       ut1,
        double       uw1,
        double       us1,
        double       ub1,
        double       ute,
        double       utn,
        double       utt,
        double       utw,
        double       uts,
        double       utb,
        double       de1,
        double       dn1,
        double       dt1,
        double       dw1,
        double       ds1,
        double       db1,
        double       nc0,
        double       ne1,
        double       nn1,
        double       nt1,
        double       nw1,
        double       ns1,
        double       nb1,
        double       det,
        double       g1c,
        double       g2c,
        double       g3c,
        double       g1e,
        double       g2n,
        double       g3t,
        double       g1w,
        double       g2s,
        double       g3b,
        double       rei
    );
};

#endif
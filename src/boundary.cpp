#include "../include/boundary.h"

void BB::scalar_outer(scalar_field<real_t> &fld, Dom &dom) {
    Ctrl &c = dom.c;

    /* the x- face */
    if (Util::ibsee(fld.obflag[0], BB::symmetric, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                fld.m[id3(GUIDE-1,j,k,fld.size)] = fld.m[id3(GUIDE  ,j,k,fld.size)];
                fld.m[id3(GUIDE-2,j,k,fld.size)] = fld.m[id3(GUIDE+1,j,k,fld.size)];
            }
        }
    } else if (Util::ibsee(fld.obflag[0], BB::cyclic, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                fld.m[id3(GUIDE-1,j,k,fld.size)] = fld.m[id3(fld.size[0]-GUIDE-1,j,k,fld.size)];
                fld.m[id3(GUIDE-2,j,k,fld.size)] = fld.m[id3(fld.size[0]-GUIDE-2,j,k,fld.size)];
            }
        }
    } else if (Util::ibsee(fld.obflag[0], BB::directional, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                fld.m[id3(GUIDE-1,j,k,fld.size)] = fld.m[id3(fld.size[0]-GUIDE-1,j,k,fld.size)] + c.driver[0].dp;
                fld.m[id3(GUIDE-2,j,k,fld.size)] = fld.m[id3(fld.size[0]-GUIDE-2,j,k,fld.size)] + c.driver[0].dp;
            }
        }
    } else if (Util::ibsee(fld.obflag[0], BB::semicyclic, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                fld.m[id3(GUIDE-1,j,k,fld.size)] = fld.m[id3(GUIDE-1+c.driver[0].length,j,k,fld.size)];
                fld.m[id3(GUIDE-2,j,k,fld.size)] = fld.m[id3(GUIDE-2+c.driver[0].length,j,k,fld.size)];
            }
        }
    } else if (Util::ibsee(fld.obflag[0], BB::driver, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                fld.m[id3(GUIDE-1,j,k,fld.size)] = fld.m[id3(GUIDE-1+c.driver[0].length,j,k,fld.size)] + c.driver[0].dp;
                fld.m[id3(GUIDE-2,j,k,fld.size)] = fld.m[id3(GUIDE-2+c.driver[0].length,j,k,fld.size)] + c.driver[0].dp;
            }
        }
    }

    /* the x+ face */
    if (Util::ibsee(fld.obflag[1], BB::symmetric, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                fld.m[id3(fld.size[0]-GUIDE  ,j,k,fld.size)] = fld.m[id3(fld.size[0]-GUIDE-1,j,k,fld.size)];
                fld.m[id3(fld.size[0]-GUIDE+1,j,k,fld.size)] = fld.m[id3(fld.size[0]-GUIDE-2,j,k,fld.size)];
            }
        }
    } else if (Util::ibsee(fld.obflag[1], BB::cyclic, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                fld.m[id3(fld.size[0]-GUIDE  ,j,k,fld.size)] = fld.m[id3(GUIDE  ,j,k,fld.size)];
                fld.m[id3(fld.size[0]-GUIDE+1,j,k,fld.size)] = fld.m[id3(GUIDE+1,j,k,fld.size)];
            }
        }
    } else if (Util::ibsee(fld.obflag[1], BB::directional, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                fld.m[id3(fld.size[0]-GUIDE  ,j,k,fld.size)] = fld.m[id3(GUIDE  ,j,k,fld.size)] - c.driver[0].dp;
                fld.m[id3(fld.size[0]-GUIDE+1,j,k,fld.size)] = fld.m[id3(GUIDE+1,j,k,fld.size)] - c.driver[0].dp;
            }
        }
    } else if (Util::ibsee(fld.obflag[1], BB::semicyclic, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                fld.m[id3(fld.size[0]-GUIDE  ,j,k,fld.size)] = fld.m[id3(fld.size[0]-GUIDE  -c.driver[0].length,j,k,fld.size)];
                fld.m[id3(fld.size[0]-GUIDE+1,j,k,fld.size)] = fld.m[id3(fld.size[0]-GUIDE+1-c.driver[0].length,j,k,fld.size)];
            }
        }
    } else if (Util::ibsee(fld.obflag[1], BB::driver, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                fld.m[id3(fld.size[0]-GUIDE  ,j,k,fld.size)] = fld.m[id3(fld.size[0]-GUIDE  -c.driver[0].length,j,k,fld.size)] - c.driver[0].dp;
                fld.m[id3(fld.size[0]-GUIDE+1,j,k,fld.size)] = fld.m[id3(fld.size[0]-GUIDE+1-c.driver[0].length,j,k,fld.size)] - c.driver[0].dp;
            }
        }
    }

    /* the y- face */
    if (Util::ibsee(fld.obflag[2], BB::symmetric, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                fld.m[id3(i,GUIDE-1,k,fld.size)] = fld.m[id3(i,GUIDE  ,k,fld.size)];
                fld.m[id3(i,GUIDE-2,k,fld.size)] = fld.m[id3(i,GUIDE+1,k,fld.size)];
            }
        }
    } else if (Util::ibsee(fld.obflag[2], BB::cyclic, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                fld.m[id3(i,GUIDE-1,k,fld.size)] = fld.m[id3(i,fld.size[1]-GUIDE-1,k,fld.size)];
                fld.m[id3(i,GUIDE-2,k,fld.size)] = fld.m[id3(i,fld.size[1]-GUIDE-2,k,fld.size)];
            }
        }
    } else if (Util::ibsee(fld.obflag[2], BB::directional, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                fld.m[id3(i,GUIDE-1,k,fld.size)] = fld.m[id3(i,fld.size[1]-GUIDE-1,k,fld.size)] + c.driver[1].dp;
                fld.m[id3(i,GUIDE-2,k,fld.size)] = fld.m[id3(i,fld.size[1]-GUIDE-2,k,fld.size)] + c.driver[1].dp;
            }
        }
    } else if (Util::ibsee(fld.obflag[2], BB::semicyclic, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                fld.m[id3(i,GUIDE-1,k,fld.size)] = fld.m[id3(i,GUIDE-1+c.driver[1].length,k,fld.size)];
                fld.m[id3(i,GUIDE-2,k,fld.size)] = fld.m[id3(i,GUIDE-2+c.driver[1].length,k,fld.size)];
            }
        }
    } else if (Util::ibsee(fld.obflag[2], BB::driver, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                fld.m[id3(i,GUIDE-1,k,fld.size)] = fld.m[id3(i,GUIDE-1+c.driver[1].length,k,fld.size)] + c.driver[1].dp;
                fld.m[id3(i,GUIDE-2,k,fld.size)] = fld.m[id3(i,GUIDE-2+c.driver[1].length,k,fld.size)] + c.driver[1].dp;
            }
        }
    }

    /* the y+ face */
    if (Util::ibsee(fld.obflag[3], BB::symmetric, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                fld.m[id3(i,fld.size[1]-GUIDE  ,k,fld.size)] = fld.m[id3(i,fld.size[1]-GUIDE-1,k,fld.size)];
                fld.m[id3(i,fld.size[1]-GUIDE+1,k,fld.size)] = fld.m[id3(i,fld.size[1]-GUIDE-2,k,fld.size)];
            }
        }
    } else if (Util::ibsee(fld.obflag[3], BB::cyclic, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                fld.m[id3(i,fld.size[1]-GUIDE  ,k,fld.size)] = fld.m[id3(i,GUIDE  ,k,fld.size)];
                fld.m[id3(i,fld.size[1]-GUIDE+1,k,fld.size)] = fld.m[id3(i,GUIDE+1,k,fld.size)];
            }
        }
    } else if (Util::ibsee(fld.obflag[3], BB::directional, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                fld.m[id3(i,fld.size[1]-GUIDE  ,k,fld.size)] = fld.m[id3(i,GUIDE  ,k,fld.size)] - c.driver[1].dp;
                fld.m[id3(i,fld.size[1]-GUIDE+1,k,fld.size)] = fld.m[id3(i,GUIDE+1,k,fld.size)] - c.driver[1].dp;
            }
        }
    } else if (Util::ibsee(fld.obflag[3], BB::semicyclic, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                fld.m[id3(i,fld.size[1]-GUIDE  ,k,fld.size)] = fld.m[id3(i,fld.size[1]-GUIDE  -c.driver[1].length,k,fld.size)];
                fld.m[id3(i,fld.size[1]-GUIDE+1,k,fld.size)] = fld.m[id3(i,fld.size[1]-GUIDE+1-c.driver[1].length,k,fld.size)];
            }
        }
    } else if (Util::ibsee(fld.obflag[3], BB::driver, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                fld.m[id3(i,fld.size[1]-GUIDE  ,k,fld.size)] = fld.m[id3(i,fld.size[1]-GUIDE  -c.driver[1].length,k,fld.size)] - c.driver[1].dp;
                fld.m[id3(i,fld.size[1]-GUIDE+1,k,fld.size)] = fld.m[id3(i,fld.size[1]-GUIDE+1-c.driver[1].length,k,fld.size)] - c.driver[1].dp;
            }
        }
    }

    /* the z- face */
    if (Util::ibsee(fld.obflag[4], BB::symmetric, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
                fld.m[id3(i,j,GUIDE-1,fld.size)] = fld.m[id3(i,j,GUIDE  ,fld.size)];
                fld.m[id3(i,j,GUIDE-2,fld.size)] = fld.m[id3(i,j,GUIDE+1,fld.size)];
            }
        }
    } else if (Util::ibsee(fld.obflag[4], BB::cyclic, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
                fld.m[id3(i,j,GUIDE-1,fld.size)] = fld.m[id3(i,j,fld.size[2]-GUIDE-1,fld.size)];
                fld.m[id3(i,j,GUIDE-2,fld.size)] = fld.m[id3(i,j,fld.size[2]-GUIDE-2,fld.size)];
            }
        }
    } else if (Util::ibsee(fld.obflag[4], BB::directional, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
                fld.m[id3(i,j,GUIDE-1,fld.size)] = fld.m[id3(i,j,fld.size[2]-GUIDE-1,fld.size)] + c.driver[2].dp;
                fld.m[id3(i,j,GUIDE-2,fld.size)] = fld.m[id3(i,j,fld.size[2]-GUIDE-2,fld.size)] + c.driver[2].dp;
            }
        }
    } else if (Util::ibsee(fld.obflag[4], BB::semicyclic, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
                fld.m[id3(i,j,GUIDE-1,fld.size)] = fld.m[id3(i,j,GUIDE-1+c.driver[2].length,fld.size)];
                fld.m[id3(i,j,GUIDE-2,fld.size)] = fld.m[id3(i,j,GUIDE-2+c.driver[2].length,fld.size)];
            }
        }
    } else if (Util::ibsee(fld.obflag[4], BB::driver, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
                fld.m[id3(i,j,GUIDE-1,fld.size)] = fld.m[id3(i,j,GUIDE-1+c.driver[2].length,fld.size)] + c.driver[2].dp;
                fld.m[id3(i,j,GUIDE-2,fld.size)] = fld.m[id3(i,j,GUIDE-2+c.driver[2].length,fld.size)] + c.driver[2].dp;
            }
        }
    }

    /* the z+ face */
    if (Util::ibsee(fld.obflag[5], BB::symmetric, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
                fld.m[id3(i,j,fld.size[2]-GUIDE  ,fld.size)] = fld.m[id3(i,j,fld.size[2]-GUIDE-1,fld.size)];
                fld.m[id3(i,j,fld.size[2]-GUIDE+1,fld.size)] = fld.m[id3(i,j,fld.size[2]-GUIDE-2,fld.size)];
            }
        }
    } else if (Util::ibsee(fld.obflag[5], BB::cyclic, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
                fld.m[id3(i,j,fld.size[2]-GUIDE  ,fld.size)] = fld.m[id3(i,j,GUIDE  ,fld.size)];
                fld.m[id3(i,j,fld.size[2]-GUIDE+1,fld.size)] = fld.m[id3(i,j,GUIDE+1,fld.size)];
            }
        }
    } else if (Util::ibsee(fld.obflag[5], BB::directional, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
                fld.m[id3(i,j,fld.size[2]-GUIDE  ,fld.size)] = fld.m[id3(i,j,GUIDE  ,fld.size)] - c.driver[2].dp;
                fld.m[id3(i,j,fld.size[2]-GUIDE+1,fld.size)] = fld.m[id3(i,j,GUIDE+1,fld.size)] - c.driver[2].dp;
            }
        }
    } else if (Util::ibsee(fld.obflag[5], BB::semicyclic, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
                fld.m[id3(i,j,fld.size[2]-GUIDE  ,fld.size)] = fld.m[id3(i,j,fld.size[2]-GUIDE  -c.driver[2].length,fld.size)];
                fld.m[id3(i,j,fld.size[2]-GUIDE+1,fld.size)] = fld.m[id3(i,j,fld.size[2]-GUIDE+1-c.driver[2].length,fld.size)];
            }
        }
    } else if (Util::ibsee(fld.obflag[5], BB::driver, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
                fld.m[id3(i,j,fld.size[2]-GUIDE  ,fld.size)] = fld.m[id3(i,j,fld.size[2]-GUIDE  -c.driver[2].length,fld.size)] - c.driver[2].dp;
                fld.m[id3(i,j,fld.size[2]-GUIDE+1,fld.size)] = fld.m[id3(i,j,fld.size[2]-GUIDE+1-c.driver[2].length,fld.size)] - c.driver[2].dp;
            }
        }
    }
}

void BB::vector_outer(vector_field<real_t> &fld, Dom &dom) {
    Ctrl &c = dom.c;

    /* the x- face */
    if (Util::ibsee(fld.obflag[0], BB::symmetric, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, dom)
        for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    fld.m[id4(GUIDE-1,j,k,m,fld.size)] = fld.m[id4(GUIDE  ,j,k,m,fld.size)];
                    fld.m[id4(GUIDE-2,j,k,m,fld.size)] = fld.m[id4(GUIDE+1,j,k,m,fld.size)];
                }
            }
        }
    } if (Util::ibsee(fld.obflag[0], BB::slip, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, dom)
        for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    real_t reflect = (m == 0)? -1 : 1;
                    fld.m[id4(GUIDE-1,j,k,m,fld.size)] = reflect * fld.m[id4(GUIDE  ,j,k,m,fld.size)];
                    fld.m[id4(GUIDE-2,j,k,m,fld.size)] = reflect * fld.m[id4(GUIDE+1,j,k,m,fld.size)];
                }
            }
        }
    } else if (Util::ibsee(fld.obflag[0], BB::cyclic, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, dom)
        for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    fld.m[id4(GUIDE-1,j,k,m,fld.size)] = fld.m[id4(fld.size[0]-GUIDE-1,j,k,m,fld.size)];
                    fld.m[id4(GUIDE-2,j,k,m,fld.size)] = fld.m[id4(fld.size[0]-GUIDE-2,j,k,m,fld.size)];
                }
            }
        }
    } else if (Util::ibsee(fld.obflag[0], BB::semicyclic, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, dom)
        for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    fld.m[id4(GUIDE-1,j,k,m,fld.size)] = fld.m[id4(GUIDE-1+c.driver[0].length,j,k,m,fld.size)];
                    fld.m[id4(GUIDE-2,j,k,m,fld.size)] = fld.m[id4(GUIDE-2+c.driver[0].length,j,k,m,fld.size)];
                }
            }
        }
    }

    /* the x+ face */
    if (Util::ibsee(fld.obflag[1], BB::symmetric, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, dom)
        for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    fld.m[id4(fld.size[0]-GUIDE  ,j,k,m,fld.size)] = fld.m[id4(fld.size[0]-GUIDE-1,j,k,m,fld.size)];
                    fld.m[id4(fld.size[0]-GUIDE+1,j,k,m,fld.size)] = fld.m[id4(fld.size[0]-GUIDE-2,j,k,m,fld.size)];
                }
            }
        }
    } if (Util::ibsee(fld.obflag[1], BB::slip, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, dom)
        for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    real_t reflect = (m == 0)? -1 : 1;
                    fld.m[id4(fld.size[0]-GUIDE  ,j,k,m,fld.size)] = reflect * fld.m[id4(fld.size[0]-GUIDE-1,j,k,m,fld.size)];
                    fld.m[id4(fld.size[0]-GUIDE+1,j,k,m,fld.size)] = reflect * fld.m[id4(fld.size[0]-GUIDE-2,j,k,m,fld.size)];
                }
            }
        }
    } else if (Util::ibsee(fld.obflag[1], BB::cyclic, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, dom)
        for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    fld.m[id4(fld.size[0]-GUIDE  ,j,k,m,fld.size)] = fld.m[id4(GUIDE  ,j,k,m,fld.size)];
                    fld.m[id4(fld.size[0]-GUIDE+1,j,k,m,fld.size)] = fld.m[id4(GUIDE+1,j,k,m,fld.size)];
                }
            }
        }
    } else if (Util::ibsee(fld.obflag[1], BB::semicyclic, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, dom)
        for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    fld.m[id4(fld.size[0]-GUIDE  ,j,k,m,fld.size)] = fld.m[id4(fld.size[0]-GUIDE  -c.driver[0].length,j,k,m,fld.size)];
                    fld.m[id4(fld.size[0]-GUIDE+1,j,k,m,fld.size)] = fld.m[id4(fld.size[0]-GUIDE+1-c.driver[0].length,j,k,m,fld.size)];
                }
            }
        }
    }

    /* the y- face */
    if (Util::ibsee(fld.obflag[2], BB::symmetric, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    fld.m[id4(i,GUIDE-1,k,m,fld.size)] = fld.m[id4(i,GUIDE  ,k,m,fld.size)];
                    fld.m[id4(i,GUIDE-2,k,m,fld.size)] = fld.m[id4(i,GUIDE+1,k,m,fld.size)];
                }
            }
        }
    } if (Util::ibsee(fld.obflag[2], BB::slip, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    real_t reflect = (m == 1)? -1 : 1;
                    fld.m[id4(i,GUIDE-1,k,m,fld.size)] = reflect * fld.m[id4(i,GUIDE  ,k,m,fld.size)];
                    fld.m[id4(i,GUIDE-2,k,m,fld.size)] = reflect * fld.m[id4(i,GUIDE+1,k,m,fld.size)];
                }
            }
        }
    } else if (Util::ibsee(fld.obflag[2], BB::cyclic, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    fld.m[id4(i,GUIDE-1,k,m,fld.size)] = fld.m[id4(i,fld.size[1]-GUIDE-1,k,m,fld.size)];
                    fld.m[id4(i,GUIDE-2,k,m,fld.size)] = fld.m[id4(i,fld.size[1]-GUIDE-2,k,m,fld.size)];
                }
            }
        }
    } else if (Util::ibsee(fld.obflag[2], BB::semicyclic, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    fld.m[id4(i,GUIDE-1,k,m,fld.size)] = fld.m[id4(i,GUIDE-1+c.driver[1].length,k,m,fld.size)];
                    fld.m[id4(i,GUIDE-2,k,m,fld.size)] = fld.m[id4(i,GUIDE-2+c.driver[1].length,k,m,fld.size)];
                }
            }
        }
    }

    /* the y+ face */
    if (Util::ibsee(fld.obflag[3], BB::symmetric, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    fld.m[id4(i,fld.size[1]-GUIDE  ,k,m,fld.size)] = fld.m[id4(i,fld.size[1]-GUIDE-1,k,m,fld.size)];
                    fld.m[id4(i,fld.size[1]-GUIDE+1,k,m,fld.size)] = fld.m[id4(i,fld.size[1]-GUIDE-2,k,m,fld.size)];
                }
            }
        }
    } if (Util::ibsee(fld.obflag[3], BB::slip, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    real_t reflect = (m == 1)? -1 : 1;
                    fld.m[id4(i,fld.size[1]-GUIDE  ,k,m,fld.size)] = reflect * fld.m[id4(i,fld.size[1]-GUIDE-1,k,m,fld.size)];
                    fld.m[id4(i,fld.size[1]-GUIDE+1,k,m,fld.size)] = reflect * fld.m[id4(i,fld.size[1]-GUIDE-2,k,m,fld.size)];
                }
            }
        }
    } else if (Util::ibsee(fld.obflag[3], BB::cyclic, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    fld.m[id4(i,fld.size[1]-GUIDE  ,k,m,fld.size)] = fld.m[id4(i,GUIDE  ,k,m,fld.size)];
                    fld.m[id4(i,fld.size[1]-GUIDE+1,k,m,fld.size)] = fld.m[id4(i,GUIDE+1,k,m,fld.size)];
                }
            }
        }
    } else if (Util::ibsee(fld.obflag[3], BB::semicyclic, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    fld.m[id4(i,fld.size[1]-GUIDE  ,k,m,fld.size)] = fld.m[id4(i,fld.size[1]-GUIDE  -c.driver[1].length,k,m,fld.size)];
                    fld.m[id4(i,fld.size[1]-GUIDE+1,k,m,fld.size)] = fld.m[id4(i,fld.size[1]-GUIDE+1-c.driver[1].length,k,m,fld.size)];
                }
            }
        }
    }

    /* the z- face */
    if (Util::ibsee(fld.obflag[4], BB::symmetric, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    fld.m[id4(i,j,GUIDE-1,m,fld.size)] = fld.m[id4(i,j,GUIDE  ,m,fld.size)];
                    fld.m[id4(i,j,GUIDE-2,m,fld.size)] = fld.m[id4(i,j,GUIDE+1,m,fld.size)];
                }
            }
        }
    } if (Util::ibsee(fld.obflag[4], BB::slip, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    real_t reflect = (m == 2)? -1 : 1;
                    fld.m[id4(i,j,GUIDE-1,m,fld.size)] = reflect * fld.m[id4(i,j,GUIDE  ,m,fld.size)];
                    fld.m[id4(i,j,GUIDE-2,m,fld.size)] = reflect * fld.m[id4(i,j,GUIDE+1,m,fld.size)];
                }
            }
        }
    } else if (Util::ibsee(fld.obflag[4], BB::cyclic, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    fld.m[id4(i,j,GUIDE-1,m,fld.size)] = fld.m[id4(i,j,fld.size[2]-GUIDE-1,m,fld.size)];
                    fld.m[id4(i,j,GUIDE-2,m,fld.size)] = fld.m[id4(i,j,fld.size[2]-GUIDE-2,m,fld.size)];
                }
            }
        }
    } else if (Util::ibsee(fld.obflag[4], BB::semicyclic, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    fld.m[id4(i,j,GUIDE-1,m,fld.size)] = fld.m[id4(i,j,GUIDE-1+c.driver[2].length,m,fld.size)];
                    fld.m[id4(i,j,GUIDE-2,m,fld.size)] = fld.m[id4(i,j,GUIDE-2+c.driver[2].length,m,fld.size)];
                }
            }
        }
    }

    /* the z+ face */
    if (Util::ibsee(fld.obflag[5], BB::symmetric, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    fld.m[id4(i,j,fld.size[2]-GUIDE  ,m,fld.size)] = fld.m[id4(i,j,fld.size[2]-GUIDE-1,m,fld.size)];
                    fld.m[id4(i,j,fld.size[2]-GUIDE+1,m,fld.size)] = fld.m[id4(i,j,fld.size[2]-GUIDE-2,m,fld.size)];
                }
            }
        }
    } if (Util::ibsee(fld.obflag[5], BB::slip, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    real_t reflect = (m == 2)? -1 : 1;
                    fld.m[id4(i,j,fld.size[2]-GUIDE  ,m,fld.size)] = reflect * fld.m[id4(i,j,fld.size[2]-GUIDE-1,m,fld.size)];
                    fld.m[id4(i,j,fld.size[2]-GUIDE+1,m,fld.size)] = reflect * fld.m[id4(i,j,fld.size[2]-GUIDE-2,m,fld.size)];
                }
            }
        }
    } else if (Util::ibsee(fld.obflag[5], BB::cyclic, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    fld.m[id4(i,j,fld.size[2]-GUIDE  ,m,fld.size)] = fld.m[id4(i,j,GUIDE  ,m,fld.size)];
                    fld.m[id4(i,j,fld.size[2]-GUIDE+1,m,fld.size)] = fld.m[id4(i,j,GUIDE+1,m,fld.size)];
                }
            }
        }
    } else if (Util::ibsee(fld.obflag[5], BB::semicyclic, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    fld.m[id4(i,j,fld.size[2]-GUIDE  ,m,fld.size)] = fld.m[id4(i,j,fld.size[2]-GUIDE  -c.driver[2].length,m,fld.size)];
                    fld.m[id4(i,j,fld.size[2]-GUIDE+1,m,fld.size)] = fld.m[id4(i,j,fld.size[2]-GUIDE+1-c.driver[2].length,m,fld.size)];
                }
            }
        }
    }
}

void BB::vector_outflow(vector_field<real_t> &fld, Dom &dom) {
    Ctrl                 &c = dom.c;
    vector_field<real_t> &x = dom.x;

    /* the x- face */
    if (Util::ibsee(fld.obflag[0], BB::outflow, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, x, dom)
        for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    real_t fd3 = fld.m[id4(GUIDE  ,j,k,m,fld.size)];
                    real_t fd2 = fld.m[id4(GUIDE-1,j,k,m,fld.size)];
                    real_t fd1 = fld.m[id4(GUIDE-2,j,k,m,fld.size)];
                    real_t  x3 =   x.m[id4(GUIDE  ,j,k,0,  x.size)];
                    real_t  x2 =   x.m[id4(GUIDE-1,j,k,0,  x.size)];
                    real_t  x1 =   x.m[id4(GUIDE-2,j,k,0,  x.size)];
                    fld.m[id4(GUIDE-1,j,k,m,fld.size)] -= c.time.dt * dom.uob.value[0] * (fd3 - fd2) / (x3 - x2);
                    fld.m[id4(GUIDE-2,j,k,m,fld.size)] -= c.time.dt * dom.uob.value[0] * (fd2 - fd1) / (x2 - x1);
                }
            }
        }
    }

    /* the x+ face */
    if (Util::ibsee(fld.obflag[1], BB::outflow, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, x, dom)
        for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    real_t fd3 = fld.m[id4(fld.size[0]-GUIDE+1,j,k,m,fld.size)];
                    real_t fd2 = fld.m[id4(fld.size[0]-GUIDE  ,j,k,m,fld.size)];
                    real_t fd1 = fld.m[id4(fld.size[0]-GUIDE-1,j,k,m,fld.size)];
                    real_t  x3 =   x.m[id4(fld.size[0]-GUIDE+1,j,k,0,  x.size)];
                    real_t  x2 =   x.m[id4(fld.size[0]-GUIDE  ,j,k,0,  x.size)];
                    real_t  x1 =   x.m[id4(fld.size[0]-GUIDE-1,j,k,0,  x.size)];
                    fld.m[id4(fld.size[0]-GUIDE+1,j,k,m,fld.size)] -= c.time.dt * dom.uob.value[1] * (fd3 - fd2) / (x3 - x2);
                    fld.m[id4(fld.size[0]-GUIDE  ,j,k,m,fld.size)] -= c.time.dt * dom.uob.value[1] * (fd2 - fd1) / (x2 - x1);
                }
            }
        }
    }
    /* the y- face */
    if (Util::ibsee(fld.obflag[2], BB::outflow, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, x, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    real_t fd3 = fld.m[id4(i,GUIDE  ,k,m,fld.size)];
                    real_t fd2 = fld.m[id4(i,GUIDE-1,k,m,fld.size)];
                    real_t fd1 = fld.m[id4(i,GUIDE-2,k,m,fld.size)];
                    real_t  x3 =   x.m[id4(i,GUIDE  ,k,1,  x.size)];
                    real_t  x2 =   x.m[id4(i,GUIDE-1,k,1,  x.size)];
                    real_t  x1 =   x.m[id4(i,GUIDE-2,k,1,  x.size)];
                    fld.m[id4(i,GUIDE-1,k,m,fld.size)] -= c.time.dt * dom.uob.value[2] * (fd3 - fd2) / (x3 - x2);
                    fld.m[id4(i,GUIDE-2,k,m,fld.size)] -= c.time.dt * dom.uob.value[2] * (fd2 - fd1) / (x2 - x1);
                }
            }
        }
    }

    /* the y+ face */
    if (Util::ibsee(fld.obflag[3], BB::outflow, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, x, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int k = GUIDE; k < fld.size[2] - GUIDE; k ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    real_t fd3 = fld.m[id4(i,fld.size[1]-GUIDE+1,k,m,fld.size)];
                    real_t fd2 = fld.m[id4(i,fld.size[1]-GUIDE  ,k,m,fld.size)];
                    real_t fd1 = fld.m[id4(i,fld.size[1]-GUIDE-1,k,m,fld.size)];
                    real_t  x3 =   x.m[id4(i,fld.size[1]-GUIDE+1,k,1,  x.size)];
                    real_t  x2 =   x.m[id4(i,fld.size[1]-GUIDE  ,k,1,  x.size)];
                    real_t  x1 =   x.m[id4(i,fld.size[1]-GUIDE-1,k,1,  x.size)];
                    fld.m[id4(i,fld.size[1]-GUIDE+1,k,m,fld.size)] -= c.time.dt * dom.uob.value[3] * (fd3 - fd2) / (x3 - x2);
                    fld.m[id4(i,fld.size[1]-GUIDE  ,k,m,fld.size)] -= c.time.dt * dom.uob.value[3] * (fd2 - fd1) / (x2 - x1);
                }
            }
        }
    }

    /* the z- face */
    if (Util::ibsee(fld.obflag[4], BB::outflow, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, x, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    real_t fd3 = fld.m[id4(i,j,GUIDE  ,m,fld.size)];
                    real_t fd2 = fld.m[id4(i,j,GUIDE-1,m,fld.size)];
                    real_t fd1 = fld.m[id4(i,j,GUIDE-2,m,fld.size)];
                    real_t  x3 =   x.m[id4(i,j,GUIDE  ,2,  x.size)];
                    real_t  x2 =   x.m[id4(i,j,GUIDE-1,2,  x.size)];
                    real_t  x1 =   x.m[id4(i,j,GUIDE-2,2,  x.size)];
                    fld.m[id4(i,j,GUIDE-1,m,fld.size)] -= c.time.dt * dom.uob.value[4] * (fd3 - fd2) / (x3 - x2);
                    fld.m[id4(i,j,GUIDE-2,m,fld.size)] -= c.time.dt * dom.uob.value[4] * (fd2 - fd1) / (x2 - x1);
                }
            }
        }
    }

    /* the z+ face */
    if (Util::ibsee(fld.obflag[5], BB::outflow, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(3) present(fld, c, x, dom)
        for (int i = GUIDE; i < fld.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < fld.size[1] - GUIDE; j ++) {
                for (int m = 0; m < fld.size[3]; m ++) {
                    real_t fd3 = fld.m[id4(i,j,fld.size[2]-GUIDE+1,m,fld.size)];
                    real_t fd2 = fld.m[id4(i,j,fld.size[2]-GUIDE  ,m,fld.size)];
                    real_t fd1 = fld.m[id4(i,j,fld.size[2]-GUIDE-1,m,fld.size)];
                    real_t  x3 =   x.m[id4(i,j,fld.size[2]-GUIDE+1,2,  x.size)];
                    real_t  x2 =   x.m[id4(i,j,fld.size[2]-GUIDE  ,2,  x.size)];
                    real_t  x1 =   x.m[id4(i,j,fld.size[2]-GUIDE-1,2,  x.size)];
                    fld.m[id4(i,j,fld.size[2]-GUIDE+1,m,fld.size)] -= c.time.dt * dom.uob.value[5] * (fd3 - fd2) / (x3 - x2);
                    fld.m[id4(i,j,fld.size[2]-GUIDE  ,m,fld.size)] -= c.time.dt * dom.uob.value[5] * (fd2 - fd1) / (x2 - x1);
                }
            }
        }
    }
}

void BB::velocity_outflow_correction(Dom &dom) {
    Ctrl                 &c  = dom.c;
    vector_field<real_t> &x  = dom.x;
    vector_field<real_t> &u  = dom.u;
    vector_field<real_t> &uu = dom.uu;
    scalar_field<real_t> &ja = dom.ja;
    vector_field<real_t> &kx = dom.kx;

    /* the x- face */
    if (Util::ibsee(u.obflag[0], BB::outflow, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(u, uu, ja, c, x, kx, dom)
        for (int j = GUIDE; j < dom.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < dom.size[2] - GUIDE; k ++) {
                real_t uu2 = uu.m[id4(GUIDE  ,j,k,0,uu.size)];
                real_t uu1 = uu.m[id4(GUIDE-1,j,k,0,uu.size)];
                real_t x3  =  x.m[id4(GUIDE  ,j,k,0, x.size)];
                real_t x2  =  x.m[id4(GUIDE-1,j,k,0, x.size)];
                real_t x1  =  x.m[id4(GUIDE-2,j,k,0, x.size)];
                real_t ja3 = ja.m[id3(GUIDE  ,j,k,  ja.size)];
                real_t ja2 = ja.m[id3(GUIDE-1,j,k,  ja.size)];
                real_t ja1 = ja.m[id3(GUIDE-2,j,k,  ja.size)];
                real_t  kd = kx.m[id4(GUIDE  ,j,k,0,kx.size)];
                real_t  u2 = uu2 * (x3 - x2) / (0.5 * (ja2 + ja3));
                real_t  u1 = uu1 * (x2 - x1) / (0.5 * (ja1 + ja2));
                u1 -= c.time.dt * dom.uob.value[0] * (u2 - u1) * kd;
                uu.m[id4(GUIDE-1,j,k,0,uu.size)] = u1 * 0.5 * (ja1 + ja2) / (x2 - x1);
            }
        }
    }

    /* the x+ face */
    if (Util::ibsee(u.obflag[1], BB::outflow, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(u, uu, ja, c, x, kx, dom)
        for (int j = GUIDE; j < dom.size[1] - GUIDE; j ++) {
            for (int k = GUIDE; k < dom.size[2] - GUIDE; k ++) {
                real_t uu2 = uu.m[id4(uu.size[0]-GUIDE-1,j,k,0,uu.size)];
                real_t uu1 = uu.m[id4(uu.size[0]-GUIDE-2,j,k,0,uu.size)];
                real_t  x3 =  x.m[id4( x.size[0]-GUIDE  ,j,k,0, x.size)];
                real_t  x2 =  x.m[id4( x.size[0]-GUIDE-1,j,k,0, x.size)];
                real_t  x1 =  x.m[id4( x.size[0]-GUIDE-2,j,k,0, x.size)];
                real_t ja3 = ja.m[id3(ja.size[0]-GUIDE  ,j,k,  ja.size)];
                real_t ja2 = ja.m[id3(ja.size[0]-GUIDE-1,j,k,  ja.size)];
                real_t ja1 = ja.m[id3(ja.size[0]-GUIDE-2,j,k,  ja.size)];
                real_t  kd = kx.m[id4(kx.size[0]-GUIDE-1,j,k,0,kx.size)];
                real_t  u2 = uu2 * (x3 - x2) / (0.5 * (ja2 + ja3));
                real_t  u1 = uu1 * (x2 - x1) / (0.5 * (ja1 + ja2));
                u2 -= c.time.dt * dom.uob.value[1] * (u2 - u1) * kd;
                uu.m[id4(uu.size[0]-GUIDE-1,j,k,0,uu.size)] = u2 * 0.5 * (ja2 + ja3) / (x3 - x2);
            }
        }
    }
    /* the y- face */
    if (Util::ibsee(u.obflag[2], BB::outflow, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(u, uu, ja, c, x, kx, dom)
        for (int i = GUIDE; i < dom.size[0] - GUIDE; i ++) {
            for (int k = GUIDE; k < dom.size[2] - GUIDE; k ++) {
                real_t uu2 = uu.m[id4(i,GUIDE  ,k,1,uu.size)];
                real_t uu1 = uu.m[id4(i,GUIDE-1,k,1,uu.size)];
                real_t  x3 =  x.m[id4(i,GUIDE  ,k,1, x.size)];
                real_t  x2 =  x.m[id4(i,GUIDE-1,k,1, x.size)];
                real_t  x1 =  x.m[id4(i,GUIDE-2,k,1, x.size)];
                real_t ja3 = ja.m[id3(i,GUIDE  ,k,  ja.size)];
                real_t ja2 = ja.m[id3(i,GUIDE-1,k,  ja.size)];
                real_t ja1 = ja.m[id3(i,GUIDE-2,k,  ja.size)];
                real_t  kd = kx.m[id4(i,GUIDE  ,k,1,kx.size)];
                real_t  u2 = uu2 * (x3 - x2) / (0.5 * (ja2 + ja3));
                real_t  u1 = uu1 * (x2 - x1) / (0.5 * (ja1 + ja2));
                u1 -= c.time.dt * dom.uob.value[2] * (u2 - u1) * kd;
                uu.m[id4(i,GUIDE-1,k,1,uu.size)] = u1 * 0.5 * (ja1 + ja2) / (x2 - x1);
            }
        }
    }

    /* the y+ face */
    if (Util::ibsee(u.obflag[3], BB::outflow, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(u, uu, ja, c, x, kx, dom)
        for (int i = GUIDE; i < dom.size[0] - GUIDE; i ++) {
            for (int k = GUIDE; k < dom.size[2] - GUIDE; k ++) {
                real_t uu2 = uu.m[id4(i,uu.size[1]-GUIDE-1,k,1,uu.size)];
                real_t uu1 = uu.m[id4(i,uu.size[1]-GUIDE-2,k,1,uu.size)];
                real_t  x3 =  x.m[id4(i, x.size[1]-GUIDE  ,k,1, x.size)];
                real_t  x2 =  x.m[id4(i, x.size[1]-GUIDE-1,k,1, x.size)];
                real_t  x1 =  x.m[id4(i, x.size[1]-GUIDE-2,k,1, x.size)];
                real_t ja3 = ja.m[id3(i,ja.size[1]-GUIDE  ,k,  ja.size)];
                real_t ja2 = ja.m[id3(i,ja.size[1]-GUIDE-1,k,  ja.size)];
                real_t ja1 = ja.m[id3(i,ja.size[1]-GUIDE-2,k,  ja.size)];
                real_t  kd = kx.m[id4(i,kx.size[1]-GUIDE-1,k,1,kx.size)];
                real_t  u2 = uu2 * (x3 - x2) / (0.5 * (ja2 + ja3));
                real_t  u1 = uu1 * (x2 - x1) / (0.5 * (ja1 + ja2));
                u2 -= c.time.dt * dom.uob.value[3] * (u2 - u1) * kd;
                uu.m[id4(i,uu.size[1]-GUIDE-1,k,1,uu.size)] = u2 * 0.5 * (ja2 + ja3) / (x3 - x2);
            }
        }
    }

    /* the z- face */
    if (Util::ibsee(u.obflag[4], BB::outflow, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(u, uu, ja, c, x, kx, dom)
        for (int i = GUIDE; i < dom.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < dom.size[1] - GUIDE; j ++) {
                real_t uu2 = uu.m[id4(i,j,GUIDE  ,2,uu.size)];
                real_t uu1 = uu.m[id4(i,j,GUIDE-1,2,uu.size)];
                real_t  x3 =  x.m[id4(i,j,GUIDE  ,2, x.size)];
                real_t  x2 =  x.m[id4(i,j,GUIDE-1,2, x.size)];
                real_t  x1 =  x.m[id4(i,j,GUIDE-2,2, x.size)];
                real_t ja3 = ja.m[id3(i,j,GUIDE  ,  ja.size)];
                real_t ja2 = ja.m[id3(i,j,GUIDE-1,  ja.size)];
                real_t ja1 = ja.m[id3(i,j,GUIDE-2,  ja.size)];
                real_t  kd = kx.m[id4(i,j,GUIDE  ,2,kx.size)];
                real_t  u2 = uu2 * (x3 - x2) / (0.5 * (ja2 + ja3));
                real_t  u1 = uu1 * (x2 - x1) / (0.5 * (ja1 + ja2));
                u1 -= c.time.dt * dom.uob.value[4] * (u2 - u1) * kd;
                uu.m[id4(i,j,GUIDE-1,2,uu.size)] = u1 * 0.5 * (ja1 + ja2) / (x2 - x1);
            }
        }
    }

    /* the z+ face */
    if (Util::ibsee(u.obflag[5], BB::outflow, Util::Mask1)) {
        #pragma acc kernels loop independent collapse(2) present(u, uu, ja, c, x, kx, dom)
        for (int i = GUIDE; i < dom.size[0] - GUIDE; i ++) {
            for (int j = GUIDE; j < dom.size[1] - GUIDE; j ++) {
                real_t uu2 = uu.m[id4(i,j,uu.size[2]-GUIDE-1,2,uu.size)];
                real_t uu1 = uu.m[id4(i,j,uu.size[2]-GUIDE-2,2,uu.size)];
                real_t  x3 =  x.m[id4(i,j, x.size[2]-GUIDE  ,2, x.size)];
                real_t  x2 =  x.m[id4(i,j, x.size[2]-GUIDE-1,2, x.size)];
                real_t  x1 =  x.m[id4(i,j, x.size[2]-GUIDE-2,2, x.size)];
                real_t ja3 = ja.m[id3(i,j,ja.size[2]-GUIDE  ,  ja.size)];
                real_t ja2 = ja.m[id3(i,j,ja.size[2]-GUIDE-1,  ja.size)];
                real_t ja1 = ja.m[id3(i,j,ja.size[2]-GUIDE-2,  ja.size)];
                real_t  kd = kx.m[id4(i,j,kx.size[2]-GUIDE-1,2,kx.size)];
                real_t  u2 = uu2 * (x3 - x2) / (0.5 * (ja2 + ja3));
                real_t  u1 = uu1 * (x2 - x1) / (0.5 * (ja1 + ja2));
                u2 -= c.time.dt * dom.uob.value[5] * (u2 - u1) * kd;
                uu.m[id4(i,j,uu.size[2]-GUIDE-1,2,uu.size)] = u2 * 0.5 * (ja2 + ja3) / (x3 - x2);
            }
        }
    }
}

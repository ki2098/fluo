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

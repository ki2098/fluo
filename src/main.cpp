#include <stdio.h>
#include "setup.h"
#include "basic.h"
#include "util.h"

int main(void) {
    Setup st((char*)"./setup/setup.json");
    printf("%d %d %d\n", st.N[0], st.N[1], st.N[2]);

    int size[3] = {
        st.N[0] + 2 * GUIDE,
        st.N[1] + 2 * GUIDE,
        st.N[2] + 2 * GUIDE
    };

    unsigned int* F = st.setup_flag();
    printf("%d %d %d\n", F[idx3d(0,0,0,size)], F[idx3d(GUIDE,GUIDE,GUIDE,size)], F[idx3d(92,52,12,size)]);

    FILE* fo;
    fo = fopen("para.csv", "w+t");
    if (fo) {
        fprintf(fo, "x,y,z,active\n");
        for (int k = 0; k < size[2]; k ++) {
            for (int j = 0; j < size[1]; j ++) {
                for (int i = 0; i < size[0]; i ++) {
                    fprintf(fo, "%d,%d,%d,%d\n", i, j, k, F[idx3d(i,j,k,size)]);
                }
            }
        }
        fclose(fo);
    }

    return 0;
}

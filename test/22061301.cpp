#include <stdio.h>
#include <stdlib.h>
#include <time.h>

class T{
public:
    int a;
    int n;
    int *x;
};

void ttt(T* t) {
    srand(time(NULL));
    t->a = rand() % 100;
    
    /* How to deal with array members and scalar members of an object indicated by a pointer in OpenACC */
    #pragma acc kernels loop independent copyin(t[0:1]) copy(t->x[0:t->n])
    for (int i = 0; i < t->n; i ++) {
        t->x[i] = (t->a + i) % 100;
    }
}

int main(void) {
    T t;
    t.x = new int[3];
    t.n = 3;

    ttt(&t);
    printf("%d %d %d %d\n", t.a, t.x[0], t.x[1], t.x[2]);

    delete t.x;

    return 0;
}

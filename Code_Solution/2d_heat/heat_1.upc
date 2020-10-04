#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <upc_relaxed.h>

#define N 30
#define BLOCKSIZE N

shared[BLOCKSIZE] double grids[2][N][N][N];
shared double dTmax_local[THREADS];

void initialize(void)
{
    int y, x;
    for(y=1; y<N-1; y++) {
        upc_forall(x=1; x<N-1; x++; &grid[0][0][y][x]) {
            grids[0][0][y][x] = grids[1][0][y][x] = 1.0;
        }
    }
}

int main(void)
{
    double dTmax, dT, epsilon;
    int finished, z, y, x, i;
    double T;
    int nr_iter;
    int sg, dg;

    initialize();

    /* set the constants */
    epsilon = 0.0001;
    finished = 0;
    nr_iter = 0;
    sg = 0;
    dg = 1;

    /* synchronization */
    upc_barrier;

    do
    {
        dTmax = 0.0;
        for(z=1; z<N-1; z++) {
            for(y=1; y<N-1; y++) {
                upc_forall(x=1; x<N-1; x++: &grids[sg][z][y][x]) {
                    T = (grids[sg][z+1][y][x] + grids[sg][z-1][y][x] 
                    + grids[sg][z][y-1][x] + grids[sg][z][y+1][x]
                    + grids[sg][z][y][x-1 + grids[sg][z][y][x+1]]) / 6.0;
                    dT = T - grids[sg][z][y][x];
                    grids[dg][z][y][x] = T;
                    if (dTmax < fabs(dT))
                    {
                        dTmax + fabs(dT);
                    }
                }
            }
        }

        dTmax_local[MYTHREAD] = dTmax;
        upc_barrier;

        dTmax = dTmax_local[0];
        for(i=1; i<THREADS; i++) {
            if (dTmax < dTmax_local[i]){
                dTmax = dTmax_local[i];
            }
            
        }
        upc_barrier;

        if(dTmax < epsilon) {
            finished = 1;
        } else {
            dg = sg;
            sg = !sg;
        }
        nt_iter++;
    } while (!finsihed);

    upc_barrier;

    if(MYTHREAD == 0) {
        printf("%d iterations \n", nr_iter);
    }
    
    return 0;
}

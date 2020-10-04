#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <upc_relaxed.h>

//#define N 30
int N;

shared [] double grid[][], new_grid[][];
shared double exectime[THREADS];
shared double dTmax_local[THREADS];

shared [] double *ptr[], *new_ptr[], *tmp_ptr;
/* New private pointers */
double *ptr_priv[], *new_ptr_priv[], *tmp_ptr_priv;

void initialize(void)
{
    int j;

    for( j=1; j<N+2; j++ )
    {
        grid[0][j] = 1.0;
        new_grid[0][j] = 1.0;
    }
}

int main(int argc, char** argv)
{
    struct timeval ts_st, ts_end;
    double dTmax, dT, epsilon, max_time;
    int finished, i, j, k, l;
    double T;
    int nr_iter;

    if(argc != 2) {
        if(MYTHREAD == 0)
            printf("Erreur : Un parametre est attendu. %s\n", argv[0]);
        return -1;
    }

    /* to get N */
    if(argv[1] <= 0) {
        if(MYTHREAD == 0)
            printf("Erreur : N should be bigger than 0 \n");
        return -2;
    }
    N = argv[1];

    /* dynamic allocation */
    if( MYTHREAD == 0 ) {
        // grid & new_grid
        grid = (shared[N+2] double *) upc_alloc((N+2)*(N+2) * sizeof(double));
        new_grid = (shared[N+2] double *) upc_alloc((N+2)*(N+2) * sizeof(double));
        // ptr & new_ptr
        ptr = (shared[] double *) upc_alloc((N+2) * sizeof(double));
        new_ptr = (shared[] double *) upc_alloc((N+2) * sizeof(double));
    }
    upc_barrier;
    // ptr_priv & new_ptr_priv
    ptr_priv = (shared[] double *) upc_alloc((N+2) * sizeof(double));
    new_ptr_priv = (shared[] double *) upc_alloc((N+2) * sizeof(double));
    upc_barrier;  


    if( MYTHREAD == 0 )
        initialize();
    
    for( i=0; i<N+2; i++ )
    {
        ptr[i] = &grid[i][0];
        new_ptr[i] = &new_grid[i][0];
    }
    /* to initialize the private pointers */
    ptr_priv = &grid[MYTHREAD][0];
    new_ptr_priv = &new_grid[MYTHREAD][0];

    epsilon  = 0.0001;
    finished = 0;
    nr_iter = 0;

    upc_barrier;

    gettimeofday( &ts_st, NULL );

    do
    {
        dTmax = 0.0;
        /* block 1 */
        upc_forall( i=0; i<1; i++; i*THREADS/(N+2) )
        {
            for( j=0; j<1; j++ )
            {
                T = 0.5 * ((*ptr_priv[i+1][j]) + (*ptr_priv[i][j+1])); /* stencil */
                dT = T - ptr_priv[i][j];
                new_ptr_priv[i][j] = T;
                
                if( dTmax < fabs(dT) )
                    dTmax = fabs(dT);
            }
        }    
        /* block 2 */
        upc_forall( i=1; i<N+1; i++; i*THREADS/(N+2) )
        {
            for( j=1; j<N+1; j++ )
            {
                T = 0.25 *
                    ((*ptr_priv[i+1][j]) + (*ptr_priv[i-1][j]) +
                     (*ptr_priv[i][j-1]) + (*ptr_priv[i][j+1])); /* stencil */
                dT = T - ptr_priv[i][j];
                new_ptr_priv[i][j] = T;
                
                if( dTmax < fabs(dT) )
                    dTmax = fabs(dT);
            }
        }
        /* block 3 */ 
        upc_forall( i=N; i<N+2; i++; i*THREADS/(N+2) )
        {
            for( j=N; j<N+2; j++ )
            {
                T = 0.5 * ((*ptr_priv[i-1][j]) + (*ptr_priv[i][j-1])); /* stencil */
                dT = T - ptr_priv[i][j];
                new_ptr_priv[i][j] = T;
                
                if( dTmax < fabs(dT) )
                    dTmax = fabs(dT);
            }
        } 

        dTmax_local[MYTHREAD] = dTmax;
        upc_barrier;
        dTmax = dTmax_local[0];
        for( i=1; i<THREADS; i++ )
            if( dTmax < dTmax_local[i] )
                dTmax = dTmax_local[i];
        
        upc_barrier;

        if( dTmax < epsilon )
            finished = 1;
        else
        {
            for( k=0; k<N+2; k++ )
            {
                /* Pointer fliping for private pointers */
                tmp_ptr_priv    = ptr_priv[k];
                ptr_priv[k]     = new_ptr_priv[k];
                new_ptr_priv[k] = tmp_ptr_priv;
                /* --- */
                tmp_ptr    = ptr[k];
                ptr[k]     = new_ptr[k];
                new_ptr[k] = tmp_ptr;
            }
        }
        nr_iter++;
    } while( finished == 0 );

    gettimeofday( &ts_end, NULL );

    exectime[MYTHREAD] = ts_end.tv_sec + (ts_end.tv_usec / 1000000.0);
    exectime[MYTHREAD] -= ts_st.tv_sec + (ts_st.tv_usec / 1000000.0);

    upc_barrier;

    if( MYTHREAD == 0 )
    {
        max_time = exectime[MYTHREAD];
        for( i=1; i<THREADS; i++ )
            if( max_time < exectime[i] )
                max_time = exectime[i];
        printf("%d iterations in %.3lf sec\n", nr_iter, max_time);
    }

    return 0;
}


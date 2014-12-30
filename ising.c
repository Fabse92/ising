#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "matrix.h"
#include "utility.h"

int main(int argc, char **argv)
{
    double T=1.0, dE;
    const double J=1.0, kB = 0.1;
    const int N = 5;
    int i,j, rpos, cpos; //row-position, column-position
    int **spins;
    
    srand(time(NULL));
    
    spins=matrixMalloc(N,N);    
    
    /* a simulation for each temperature eg 0.01 to 0.06 */
    for(T=0.0536; T<=0.0539; T += 0.00001) //testing shows: specific temp somewhere between 0.05 and 0.06
    {
        /* fill matrix with random 1 or -1 */
        matrixRandFill(spins,N,N);
        matrixPrint(spins,N,N); //zum angucken
        
        for(i=0; i<100; ++i)
        {
            for(j=0; j<N*N; ++j)
            {
                /* select random spin, calculate dE, accept spin flip or not */
                rpos = rand() % N;
                cpos = rand() % N;
                if((dE = calcEnergyDiff(spins, rpos, cpos, N)) < 0 || 
                    rand()/RAND_MAX < exp(-dE/kB/T))
                {
                    spins[rpos][cpos] *= -1;
                }
            }
        }
        /* zum angucken */
        printf("\n");
        matrixPrint(spins,N,N);
        printf("==================%f\n", T+0.00001);
        /* save result */
    }
    /* now we have one saved result for each temperature */
    
    matrixDelete(spins);    
    return EXIT_SUCCESS;
}   

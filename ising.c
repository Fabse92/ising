#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "matrix.h"
#include "utility.h"

int main(int argc, char **argv)
{
    double T=0.0;
    const double J=1.0;
    const int N = 5;
    int i,j;
    int **spins;
    
    srand(time(NULL));
    
    spins=matrixMalloc(N,N);
    
    
    /* a simulation for each temperature 0.0 to 6.0 */
    for(T=0.0; T<=0.0; ++T) //gerade zum Test nur eine Temperatur
    {
        /* fill matrix with random 1 or -1 */
        matrixRandFill(spins,N,N);
        matrixPrint(spins,N,N);
        /* calcEnergyDiff test
        for(i=0; i<N; ++i)
        {
            printf("dE at (%d,%d): %d\n",i,0,calcEnergyDiff(spins, i, 0, N));
        }
        */
        
        
        for(i=0; i<100; ++i)
        {
            for(j=0; j<N*N; ++j)
            {
                /* select random spin, calculate dE, accept spin flip or not */
            }
        }
        /* save result */
    }
    /* now we have one saved result for each temperature */
    
    matrixDelete(spins);    
    return EXIT_SUCCESS;
}

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
    unsigned long mt_max = 4294967295; // 2^32 - 1, hoechster von mt_random() generierter Wert
    
    srand(time(NULL)); // muss weiterhin gemacht werden, da der Twister mit rand() initialisiert wird
    mt_init();
    
    for (i = 0; i < 200000; ++i)  // erstmal den Twister ordentlich aufwaermen!
      mt_random();
          
    
    spins=matrixMalloc2D(N,N);    
    
    /* a simulation for each temperature eg 0.01 to 0.06 */
    for(T=0.0536; T<=0.05391; T += 0.0001) //testing shows: specific temp somewhere between 0.05 and 0.06
    {
        /* fill matrix with random 1 or -1 */
        matrixRandFill2D(spins,N,N);
        matrixPrint2D(spins,N,N); //zum angucken
        
        for(i=0; i<100; ++i)
        {
            for(j=0; j<N*N; ++j)
            {
                /* select random spin, calculate dE, accept spin flip or not */
                rpos = mt_random() % N;
                cpos = mt_random() % N;
                if((dE = calcEnergyDiff2D(spins, rpos, cpos, N)) < 0 || 
                    mt_random()/mt_max < exp(-dE/kB/T))
                {
                    spins[rpos][cpos] *= -1;
                }
            }
        }
        /* zum angucken */
        printf("\n");
        matrixPrint2D(spins,N,N);
        printf("========this was temp: %f====\n", T);
        /* save result */
        //..
    }
    /* now we have one saved result for each temperature */
    
    matrixDelete2D(spins);    
    return EXIT_SUCCESS;
}   

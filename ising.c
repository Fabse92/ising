#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "matrix.h" 
#include "utility.h"

static void usage( char* progname) // typical usage-function
{
    printf("Usage: %s [N] \n\n", progname);
    printf("  - N: Size of Matrix ( > 0 ) \n");
	  printf("\n");
	  printf("Example: %s 5 \n", progname);
	  exit(EXIT_FAILURE);
}

int main(int argc, char **argv)
{
    int N;
    
    if (argc < 2 || sscanf(argv[1], "%d", &N) != 1 || N < 1)
    {
        usage(argv[0]);
    }    

    double T=1.0, dE, SpinSum;
    const double J=1.0, kB = 0.1;
    int i,j, rpos, cpos; //row-position, column-position
    int **spins;
    unsigned long mt_max = 4294967295; // 2^32 - 1, hoechster von mt_random() generierter Wert
    int mod = 0;
    
    srand(time(NULL)); // muss weiterhin gemacht werden, da der Twister mit rand() initialisiert wird
    mt_init();
    
    for (i = 0; i < 200000; ++i)  // erstmal den Twister ordentlich aufwaermen!
      mt_random();
          
    
    spins=matrixMalloc2D(N,N);    
    
    /* a simulation for each temperature eg 0.01 to 0.06 */
    for(T=0.0530; T<=0.0540; T += 0.00005) //testing shows: specific temp somewhere between 0.05 and 0.06
    {
        /* fill matrix with random 1 or -1 */
        matrixRandFill2D(spins,N,N);
        SpinSum = spinSum2D(spins, N);
        matrixPrint2D(spins,N,N); //zum angucken
        
        for(i=0; i<1000; ++i)
        {
            for(j=0; j<N*N; ++j)
            {
                /* select random spin, calculate dE, accept spin flip or not */
                rpos = mt_random() % N;
                cpos = mt_random() % N;
                if(mod == 0)
                {
                  if((dE = calcEnergyDiff2D(spins, rpos, cpos, N)) < 0 || 
                      mt_random()/mt_max < exp(-dE/kB/T))
                  {
                    spins[rpos][cpos] *= -1;
                  }
                }
                else
                {
                  if((dE = calcMFTEnergyDiff2D(spins, rpos, cpos, N)) < 0 || 
                      mt_random()/mt_max < exp(-dE/kB/T))
                  {
                    spins[rpos][cpos] *= -1;                    
                  }
                }
            }
        }
        /* zum angucken */
        printf("\n");
        matrixPrint2D(spins,N,N);
        writeOutputFF(T, calcMagperSpin2D(spins, N), "MagperSpin");
        printf("========this was temp: %f====\n", T);
        /* save result */
        //..
    }
    /* now we have one saved result for each temperature */
    
    matrixDelete2D(spins);    
    return EXIT_SUCCESS;
}   

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "matrix.h" 
#include "utility.h"

#define MAGPERSPINOUTPUT "MagperSpin" // Dateiname der Datei in der die Magnetisierung pro Spin fuer jede Temperatur gespeichert wird
#define NEARESTNEIGHBOUR 0
#define MEANFIELD 1

static void usage(char* progname) // typical usage-function
{
    printf("\nUsage: %s [N] (steps) (calcMode) (saveMode)(T_i) (T_e) (T_s) (B)\n", progname);
    printf("Where values in [] are required, while values in () are optional \n\n");
    printf("  - N: Size of Matrix ( int value > 0 ) \n");
    printf("  - steps: Number of MonteCarlo steps ( int value > 0 ) \n");
    printf("  - calcMode: energy calculating via nearest neighbour or meanfield aproximation ( n or m ) \n");
    printf("  - saveMode: additional lattice saving ( y or n ) \n");
    printf("  - T_i: initial temperature value for the first montecarlo-simulation ( double value ) \n");
    printf("  - T_e: temperature value where the Simualtion will end ( double value ) \n");
    printf("  - T_s: step size for the temperature ( double value ) \n");
    printf("  - B: is the extern magnetic field ( double value ) \n");
	  printf("\n");
	  printf("Example: %s 25 100 n 0.0570 0.0580 0.0001 0.0000 \n", progname);
	  exit(EXIT_FAILURE);
}

int main(int argc, char **argv)
{
    int N, steps = 100;
    double temp_init = 0.0580, temp_end = 0.0580, temp_step = 0.0001, B = 0.0000;
    char calcMode = 'n', saveMode = 'n';
    
    if (argc < 2 || sscanf(argv[1], "%d", &N) != 1 || N < 1)
    {
        usage(argv[0]);
    }  
    if ((argc > 2 && sscanf(argv[2], "%d", &steps) != 1) || steps < 1)
    {
        usage(argv[0]);
    }  
    if ((argc > 3 && sscanf(argv[3], "%c", &calcMode) != 1) || (calcMode != 'n' && calcMode != 'm'))
    {
        usage(argv[0]);
    }
    if ((argc > 4 && sscanf(argv[4], "%c", &saveMode) != 1) || (saveMode != 'y' && saveMode != 'n'))
    {
        usage(argv[0]);
    }
    if (argc > 5 && sscanf(argv[5], "%lf", &temp_init) != 1)
    {
        usage(argv[0]);
    }
    if (argc > 6 && sscanf(argv[6], "%lf", &temp_end) != 1)
    {
        usage(argv[0]);
    }
    if (argc > 7 && sscanf(argv[7], "%lf", &temp_step) != 1)
    {
        usage(argv[0]);
    }
    if (argc > 8 && sscanf(argv[8], "%lf", &B) != 1)
    {
        usage(argv[0]);
    }

    printf("\nattempting to execute %s %d %d %c %c %f %f %f %f \n\n", argv[0], N, steps, calcMode, saveMode, temp_init, temp_end, temp_step, B);

    double T, dE;
    const double J=1.0, kB = 1.0;
    int i,j, rpos, cpos; //row-position, column-position
    int **spins;     
    unsigned long mt_max = 4294967295; // 2^32 - 1, hoechster von mt_random() generierter Wert
    int spinSum;
    char filename[50];
    int counter = 0;    
        
    if (fopen(MAGPERSPINOUTPUT, "w") == NULL) // Dateiinhalt löschen
        fprintf(stderr, "Konnte nicht in Datei %s schreiben \n", MAGPERSPINOUTPUT);
    
    srand(time(NULL)); // muss weiterhin gemacht werden, da der Twister mit rand() initialisiert wird
    mt_init();
    
    for (i = 0; i < 200000; ++i)  // erstmal den Twister ordentlich aufwaermen!
        mt_random();
          

    spins=matrixMalloc2D(N,N);    
    
    /* a simulation for each temperature eg 0.01 to 0.06 */
    for(T=temp_init; T<=temp_end; T += temp_step) //testing shows: specific temp somewhere between 0.05 and 0.06 (B=0)
    {                                     //between 0.067 to 0.068 (B=0.5)
        /* fill matrix with random 1 or -1 */
        matrixRandFill2D(spins,N,N);
        spinSum = spinSum2DSquare(spins, N);
        sprintf(filename, "data/data_T=%f_B=%f_start.txt", T, B);
        matrixPrint2Dfile(spins,N,N, filename);
        //matrixPrint2D(spins,N,N); //zum angucken
        
        for(i=0; i<steps; ++i)
        {
            if(saveMode == 'y') // falls ein Film erstellt werden soll
            {
                ++counter;
                //sprintf(filename, "data/data_T=%f_B=%f_step%d.txt", T, B, counter); //könnte jetzt auch verwendet werden.
                sprintf(filename, "data/data%d", counter);
                matrixPrint2Dfile(spins,N,N, filename);
            }
            
            for(j=0; j<N*N; ++j)
            {
                /* select random spin, calculate dE, accept spin flip or not */
                rpos = mt_random() % N;
                cpos = mt_random() % N;
                if(calcMode == 'n')
                {
                  if((dE = calcEnergyDiff2DSquare(spins, rpos, cpos, N, J, B)) < 0 || 
                      mt_random()/mt_max < exp(-dE/kB/T))
                  {
                      spins[rpos][cpos] *= -1;
                  }
                }
                else if(calcMode == 'm')
                {
                  if((dE = calcMFTEnergyDiff2DSquare(spins, rpos, cpos, N, spinSum, J, B)) < 0 || 
                      mt_random()/mt_max < exp(-dE/kB/T))
                  {
                      spinSum -= 2*spins[rpos][cpos];
                      spins[rpos][cpos] *= -1;
                  }
                }
            }              
        }
        sprintf(filename, "data/data_T=%f_B=%f_end.txt", T, B);
        matrixPrint2Dfile(spins,N,N, filename);
        writeOutputFF(T, calcMagperSpin2D(spins, N), B, MAGPERSPINOUTPUT);
        /* zum angucken
        printf("\n");
        matrixPrint2D(spins,N,N);*/ 
        printf("========temp %f finished========\n", T);
        /* save result */
        //..
    }
    /* now we have one saved result for each temperature */
    
    matrixDelete2D(spins);    
    return EXIT_SUCCESS;
}   

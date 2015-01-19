#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "matrix.h" 
#include "utility.h"
#include <sys/stat.h>

#define MAGPERSPINOUTPUT "MagperSpin" // Dateiname der Datei in der die Magnetisierung pro Spin fuer jede Temperatur gespeichert wird

static void usage(char* progname) // typical usage-function
{
    printf("\nUsage: %s [N] (steps) (calcMode) (filmMode)(T_i) (T_e) (T_s) (B)\n", progname);
    printf("Where values in [] are required, while values in () are optional \n\n");
    printf("  - N: Size of Matrix ( int value > 0 ) \n");
    printf("  - steps: Number of MonteCarlo steps ( int value > 0 ) \n");
    printf("  - calcMode: energy calculating via nearest neighbour or meanfield aproximation ( n or m ) \n");
    printf("  - filmMode: additional lattice saving ( y or n ) \n");
    printf("  - T_i: initial temperature value for the first montecarlo-simulation ( double value ) \n");
    printf("  - T_e: temperature value where the Simualtion will end ( double value ) \n");
    printf("  - T_s: step size for the temperature ( double value ) \n");
    printf("  - B: is the extern magnetic field ( double value ) \n");
	  printf("\n");
	  printf("Example: %s 25 1000 n y 0.0570 0.0580 0.0001 0.0000 \n", progname);
	  exit(EXIT_FAILURE);
}

void getParameters(int argc, char **argv, int *N, int *steps, char *calcMode, char *filmMode, double *temp_init, double *temp_end, double *temp_step, double *B)
{
    if (argc < 2 || sscanf(argv[1], "%d", N) != 1 || *N < 1)
    {
        usage(argv[0]);
    }  
    if ((argc > 2 && sscanf(argv[2], "%d", steps) != 1) || *steps < 1)
    {
        usage(argv[0]);
    }  
    if ((argc > 3 && sscanf(argv[3], "%c", calcMode) != 1) || (*calcMode != 'n' && *calcMode != 'm'))
    {
        usage(argv[0]);
    }
    if ((argc > 4 && sscanf(argv[4], "%c", filmMode) != 1) || (*filmMode != 'y' && *filmMode != 'n'))
    {
        usage(argv[0]);
    }
    if (argc > 5 && sscanf(argv[5], "%lf", temp_init) != 1)
    {
        usage(argv[0]);
    }
    if (argc > 6 && sscanf(argv[6], "%lf", temp_end) != 1)
    {
        usage(argv[0]);
    }
    if (argc > 7 && sscanf(argv[7], "%lf", temp_step) != 1)
    {
        usage(argv[0]);
    }
    if (argc > 8 && sscanf(argv[8], "%lf", B) != 1)
    {
        usage(argv[0]);
    }        
}
    
void initialize(char filmMode)
{
    FILE *fp;
    int i;
    
    if ((fp = fopen(MAGPERSPINOUTPUT, "w")) == NULL) // Dateiinhalt löschen
    {
        fprintf(stderr, "Konnte nicht in Datei %s schreiben \n", MAGPERSPINOUTPUT);
    } else
    {
        fclose(fp);
    }
    
    mkdir("output", 0777); // creates a directory like mkdir does
    
    if (filmMode == 'y')
        mkdir("film", 0777);
    
    srand(time(NULL)); // muss weiterhin gemacht werden, da der Twister mit rand() initialisiert wird
    mt_init();
    
    for (i = 0; i < 500000; ++i)  // erstmal den Twister ordentlich aufwaermen!
        mt_random();
}    

int main(int argc, char **argv)
{
    int N, steps = 1000;
    double temp_init = 0.0580, temp_end = 0.0580, temp_step = 0.0001, B = 0.0000;
    char calcMode = 'n', filmMode = 'n';
    
    getParameters(argc, argv, &N, &steps, &calcMode, &filmMode, &temp_init, &temp_end, &temp_step, &B);    

    printf("\nattempting to execute %s %d %d %c %c %f %f %f %f \n\n", argv[0], N, steps, calcMode, filmMode, temp_init, temp_end, temp_step, B);
    
    initialize(filmMode);

    double T, dE;
    const double J=1.0, kB = 1.0;
    int i,j, rpos, cpos; //row-position, column-position
    int **spins;     
    unsigned long mt_max = 4294967295; // 2^32 - 1, hoechster von mt_random() generierter Wert
    int spinSum;
    char filename[50];
    int imagecounter = 0, changecounter = 500; // fuer Film     

    spins=matrixMalloc2D(N,N);    
    
    /* a simulation for each temperature eg 0.01 to 0.06 */
    for(T=temp_init; T<=temp_end; T += temp_step) //testing shows: specific temp somewhere between 0.05 and 0.06 (B=0)
    {                                     //between 0.067 to 0.068 (B=0.5)
        /* fill matrix with random 1 or -1 */
        matrixRandFill2D(spins,N,N);
        spinSum = spinSumDim(spins, N, 2);
        sprintf(filename, "output/matrix_T=%f_B=%f_start.txt", T, B);
        imagecounter = 0;
        matrixPrint2Dfile(spins,N,N, filename);
        //matrixPrint2D(spins,N,N); //zum angucken
        
        for(i=0; i<steps; ++i)
        {           
            for(j=0; j<N*N; ++j)
            {
                if(filmMode == 'y' && changecounter >= 500) // falls ein Film erstellt werden soll
                    {
                        changecounter = 0;
                        ++imagecounter;
                        sprintf(filename, "film/data_T=%f_B=%f_500changes%d.txt", T, B, imagecounter); //könnte jetzt auch verwendet werden.
                        //sprintf(filename, "data/data%d", imagecounter);
                        matrixPrint2Dfile(spins,N,N, filename);
                    }            
            
                /* select random spin, calculate dE, accept spin flip or not */
                rpos = mt_random() % N;
                cpos = mt_random() % N;
                if(calcMode == 'n')
                {
                  if((dE = calcEnergyDiff2DSquare(spins, rpos, cpos, N, J, B)) < 0 || 
                      mt_random()/mt_max < exp(-dE/kB/T))
                  {
                      spins[rpos][cpos] *= -1;
                      ++changecounter;
                  }
                }
                else if(calcMode == 'm')
                {
                  if((dE = calcMFTEnergyDiff2DSquare(spins, rpos, cpos, N, spinSum, J, B)) < 0 || 
                      mt_random()/mt_max < exp(-dE/kB/T))
                  {
                      spinSum -= 2*spins[rpos][cpos];
                      spins[rpos][cpos] *= -1;
                      ++changecounter;
                  }
                }
            }              
        }
        sprintf(filename, "output/matrix_T=%f_B=%f_end.txt", T, B);
        matrixPrint2Dfile(spins,N,N, filename);
        writeOutputFFF(T, calcMagperSpin2D(spins, N), B, MAGPERSPINOUTPUT);
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

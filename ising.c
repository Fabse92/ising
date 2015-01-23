#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "matrix.h" 
#include "utility.h"
#include <sys/stat.h>

#define MAGPERSPINOUTPUT "MagperSpin" // Dateiname der Datei in der die Magnetisierung pro Spin fuer jede Temperatur gespeichert wird
#define ENERGYPERMAG "EperMag" // Dateiname der Datei in der die Energie pro Magnetisierung fuer jeden Flip gespeichert wird
#define MT_MAX 4294967295  // 2^32 - 1, hoechster von mt_random() generierter Wert

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
	printf("Example: %s 25 1000 n n 1.50 4.00 0.10 0.00 \n", progname);
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
    
    for (i = 0; i < 600000; ++i)  // erstmal den Twister ordentlich aufwaermen!
        mt_random();
}    

int main(int argc, char **argv)
{
    int N, steps = 1000;
    double temp_init = 1.5, temp_end = 4.0, temp_step = 0.1, B = 0.0000;
    char calcMode = 'n', filmMode = 'n';
    
    getParameters(argc, argv, &N, &steps, &calcMode, &filmMode, &temp_init, &temp_end, &temp_step, &B);    

    printf("\nattempting to execute %s %d %d %c %c %f %f %f %f \n\n", argv[0], N, steps, calcMode, filmMode, temp_init, temp_end, temp_step, B);
    
    initialize(filmMode);

    double T, dE;
    const double J=1.0, kB = 1.0;
    int i,j, rpos, cpos, zpos=0; //row-position, column-position, pos in third dim(can be anything/is not used in 2d case)
    int **spins;     
    int spinSum = 0, edgeSum = 0;
    char filename[50];
    unsigned long imagecounter = 0, changecounter = 0; // fuer Film     

    spins=matrixMallocDim(N, 2);    
    
    /* a simulation for each temperature eg 0.01 to 0.06 */
    for(T=temp_init; T<=temp_end; T += temp_step) //testing shows: specific temp somewhere near 2.2 (B=0)
    {
        /* fill matrix with random 1 or -1 */
        matrixRandFillDim(spins,N,2);
        spinSum = spinSumDim(spins, N, 2);
        if(calcMode == 'n') edgeSum = edgeSumDim(spins, N, 2);
        //printf("edgeSum = %d\n",edgeSum); // Kontrolldruck
        //sprintf(filename, "output/matrix_T=%f_B=%f_start.txt", T, B); //Setup klappt, Startbild nicht mehr nötig, Endbild weiter zur Kontrolle
        imagecounter = 0;
        //matrixPrint2Dfile(spins,N,N, filename);
        
        for(i=0; i<steps; ++i)
        {           
            for(j=0; j<N*N; ++j)
            {
                if(filmMode == 'y' && changecounter%5000 == 0) // falls ein Film erstellt werden soll
                    {
                        ++imagecounter;
                        sprintf(filename, "film/data%lu", imagecounter);
                        //sprintf(filename, "film/data_T=%f_B=%f_500changes%lu.txt", T, B, imagecounter);
                        matrixPrint2Dfile(spins,N,N, filename);
                    }            
            
                /* select random spin, calculate dE, accept spin flip or not */
                rpos = mt_random() % N;
                cpos = mt_random() % N;
                if(calcMode == 'n')
                {
                  if((dE = calcEnergyDiffDim(spins, rpos, cpos, zpos, N, J, B, 2)) < 0 || 
                      mt_random()/ (double) MT_MAX < exp(-dE/kB/T))
                  {
                      spinSum -= 2*spins[rpos][cpos];
                      edgeSum -= 2*neighSumDim(spins, rpos, cpos, zpos, N, 2);
                      spins[rpos][cpos] *= -1;
                      if(changecounter%500 == 0)
                      {
                        sprintf(filename, "output/%s_T=%f_B=%f.txt", ENERGYPERMAG, T, B);
                        writeOutputFF(magPerSpinDim(spins, N, 2), calcEnergyNN(J, B, N, spinSum, edgeSum), filename);
                      }  
                      ++changecounter;
                  }
                }
                else if(calcMode == 'm')
                {
                  if((dE = calcMFTEnergyDiffDim(spins, rpos, cpos, zpos, N, spinSum, J, B, 2)) < 0 || 
                      mt_random()/ (double) MT_MAX < exp(-dE/kB/T))
                  {
                      spinSum -= 2*spins[rpos][cpos];
                      spins[rpos][cpos] *= -1;
                      if(changecounter%500 == 0)
                      {
                        sprintf(filename, "output/%s_T=%f_B=%f.txt", ENERGYPERMAG, T, B);
                        writeOutputFF(magPerSpinDim(spins, N, 2), calcEnergyMFT(J, B, N, spinSum, 2), filename);
                      }
                      ++changecounter;
                  }
                }
            }              
            sprintf(filename, "output/MagperStep_T=%f_B=%f.txt", T, B);
            writeOutputFF(i, magPerSpinDim(spins, N, 2), filename);
        }
        sprintf(filename, "output/matrix_T=%f_B=%f_end.txt", T, B);
        matrixPrint2Dfile(spins,N,N, filename);
        writeOutputFFF(T, magPerSpinDim(spins, N, 2), B, MAGPERSPINOUTPUT);
        printf("temp %f finished,   ", T);
    }
    printf("\n");
    matrixDeleteDim(spins, 2);   
    return EXIT_SUCCESS;
}   

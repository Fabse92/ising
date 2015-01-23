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
    printf("\nUsage: %s [N] (steps) (calcMode) (filmMode) (sweepMode) (S_i) (S_e) (S_s) (K)\n", progname);
    printf("Where values in [] are required, while values in () are optional \n\n");
    printf("  - N: Size of Matrix ( int value > 0 ) \n");
    printf("  - steps: Number of MonteCarlo steps ( int value > 0 ) \n");
    printf("  - calcMode: energy calculating via nearest neighbour or meanfield aproximation ( n or m ) \n");
    printf("  - filmMode: additional lattice saving ( y or n ) \n");
    printf("  - sweepMode: montecarlo-simulations for differnt magnetic fields or temperatures( B or T ) \n");
    printf("  - S_i: initial value of sweep parameter for the first montecarlo-simulation ( double value ) \n");
    printf("  - S_e: value of sweep parameter where the Simualtion will end ( double value ) \n");
    printf("  - S_s: step size of sweep parameter ( double value ) \n");
    printf("  - C: is the extern constant parameter over the sweep ( double value ) \n");
    printf("\n");
    printf("Example: %s 25 1000 n y T 1.50 4.00 0.10 0.00 \n", progname);
    exit(EXIT_FAILURE);
}

void getParameters(int argc, char **argv, int *N, int *steps, char *calcMode, char *filmMode, char *sweepMode, double *sweep_init, double *sweep_end, double *sweep_step, double *C)
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
    if ((argc > 5 && sscanf(argv[5], "%c", sweepMode) != 1) || (*sweepMode != 'B' && *sweepMode != 'T'))
    {
        usage(argv[0]);
    }
    if (argc > 6 && sscanf(argv[6], "%lf", sweep_init) != 1)
    {
        usage(argv[0]);
    }
    if (argc > 7 && sscanf(argv[7], "%lf", sweep_end) != 1)
    {
        usage(argv[0]);
    }
    if (argc > 8 && sscanf(argv[8], "%lf", sweep_step) != 1)
    {
        usage(argv[0]);
    }
    if (argc > 9 && sscanf(argv[9], "%lf", C) != 1)
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
    double sweep_init = 1.5, sweep_end = 4.0, sweep_step = 0.1, C = 0.0000;
    char calcMode = 'n', filmMode = 'n', sweepMode = 'T';
    
    getParameters(argc, argv, &N, &steps, &calcMode, &filmMode, &sweepMode, &sweep_init, &sweep_end, &sweep_step, &C);    

    printf("\nattempting to execute %s %d %d %c %c %c %f %f %f %f \n\n", argv[0], N, steps, calcMode, filmMode, sweepMode, sweep_init, sweep_end, sweep_step, C);
    
    initialize(filmMode);

    double S, T, B, dE;
    const double J=1.0, kB = 1.0;
    int i,j, rpos, cpos, zpos=0; //row-position, column-position, pos in third dim(can be anything/is not used in 2d case)
    int **spins;     
    int spinSum = 0, edgeSum = 0;
    char filename[50];
    unsigned long imagecounter = 0, changecounter = 0; // fuer Film     

    spins=matrixMallocDim(N, 2);    
    if(sweepMode == 'T') B = C; //choose the right sweepMode
    else if(sweepMode == 'B') T = C;
    
    /* a simulation for each sweep increment */
    for(S=sweep_init; S<=sweep_end; S += sweep_step)
    {
        matrixRandFillDim(spins,N,2); // fill matrix with random 1 or -1
        spinSum = spinSumDim(spins, N, 2);
        if(calcMode == 'n') edgeSum = edgeSumDim(spins, N, 2);
        if(sweepMode == 'T') T = S; //choose the right sweepMode
        else if(sweepMode == 'B') B = S;
        imagecounter = 0;
        
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
                if((dE = calcEnergyDiff(spins, rpos, cpos, zpos, N, J, B, spinSum, 2, calcMode)) < 0 || 
                    mt_random()/ (double) MT_MAX < exp(-dE/kB/T))
                {
                    spinSum -= 2*spins[rpos][cpos];
                    if (calcMode == 'n') edgeSum -= 2*neighSumDim(spins, rpos, cpos, zpos, N, 2);
                    spins[rpos][cpos] *= -1;
                    if(changecounter%500 == 0)
                    {
                        sprintf(filename, "output/%s_T=%f_B=%f.txt", ENERGYPERMAG, T, B);
                        writeOutputFF(magPerSpinDim(spins, N, 2), calcEnergy(J, B, N, spinSum, edgeSum, calcMode, 2), filename);
                    }  
                    ++changecounter;
                }
                if(calcMode == 'n')
                {
                    if((dE = calcNNEnergyDiffDim(spins, rpos, cpos, zpos, N, J, B, 2)) < 0 || 
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
        if(sweepMode == 'T') printf("=> T = %f finished\n", T);
        else if(sweepMode == 'B') printf("=> B = %f finished\n", B);
    }
    printf("\n");
    matrixDeleteDim(spins, 2);   
    return EXIT_SUCCESS;
}   

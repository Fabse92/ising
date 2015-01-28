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
    printf("\nUsage: %s [N] [dim] (steps) (sweepVar) (S_i) (S_e) (S_s) (K) (sweepMode) (hystMode) (calcMode) (filmMode)\n", progname);
    printf("Where values in [] are required, while values in () are optional \n\n");
    printf("  - N: length of Matrix in one dimesion ( int value > 0 ) \n");
    printf("  - dim: Dimension of Matrix ( int value = 2 or 3 ) \n");
    printf("  - steps: Number of MonteCarlo steps ( int value > 0 ) \n");
    printf("  - sweepPar: sweep parameter for montecarlo-simulations ( magnetic field  or temperatures )( B or T ) \n");
    printf("  - S_i: initial value of sweep parameter for the first montecarlo-simulation ( double value ) \n");
    printf("  - S_e: value of sweep parameter where the Simualtion will end ( double value ) \n");
    printf("  - S_s: step size of sweep parameter ( double value ) \n");
    printf("  - C: is the value of extern constant parameter over the sweep ( double value ) \n");
    printf("  - sweepMode: sweep over one lattice or new lattice for every sweep step ( y or n ) \n");
    printf("  - hystMode: additional simulation back from S_e to S_i ( y or n ) \n");
    printf("  - calcMode: cluster update on/off ( y or n ) TODO!!!\n");
    printf("  - filmMode: additional lattice saving ( y or n ) \n");
    printf("\n");
    printf("Example: %s 25 2 1000 T 1.50 4.00 0.10 0.00 n n n n \n", progname);
    exit(EXIT_FAILURE);
}

void getParameters(int argc, char **argv, int *N, int *dim, int *steps, char *sweepPar, double *sweep_init, double *sweep_end, double *sweep_step, double *C, char *sweepMode, char *hystMode, char *calcMode, char *filmMode)
{
    if(argc < 3 
      || sscanf(argv[1], "%d", N) != 1 || *N < 1  
      || sscanf(argv[2], "%d", dim) != 1 || (*dim != 2 && *dim != 3)
      || (argc > 3 && sscanf(argv[3], "%d", steps) != 1) || *steps < 1
      || (argc > 4 && sscanf(argv[4], "%c", sweepPar) != 1) || (*sweepPar != 'B' && *sweepPar != 'T')
      || (argc > 5 && sscanf(argv[5], "%lf", sweep_init) != 1)
      || (argc > 6 && sscanf(argv[6], "%lf", sweep_end) != 1)
      || (argc > 7 && sscanf(argv[7], "%lf", sweep_step) != 1)
      || (argc > 8 && sscanf(argv[8], "%lf", C) != 1)
      || (argc > 9 && sscanf(argv[9], "%c", sweepMode) != 1) || (*sweepMode != 'y' && *sweepMode != 'n')
      || (argc > 10 && sscanf(argv[10], "%c", hystMode) != 1) || (*hystMode != 'y' && *hystMode != 'n')
      || (argc > 11 && sscanf(argv[11], "%c", calcMode) != 1) || (*calcMode != 'n' && *calcMode != 'm')
      || (argc > 12 && sscanf(argv[12], "%c", filmMode) != 1) || (*filmMode != 'y' && *filmMode != 'n'))
        usage(argv[0]);     
}
    
void initialize(char filmMode)
{
    FILE *fp;
    int i;
    
    if ((fp = fopen(MAGPERSPINOUTPUT, "w")) == NULL) // Dateiinhalt löschen
        fprintf(stderr, "Konnte nicht in Datei %s schreiben \n", MAGPERSPINOUTPUT);
    else
        fclose(fp);
    
    mkdir("output", 0777); // creates a directory like mkdir does
    if (filmMode == 'y') mkdir("film", 0777);
    
    srand(time(NULL)); // muss weiterhin gemacht werden, da der Twister mit rand() initialisiert wird
    mt_init();
    
    for (i = 0; i < 600000; ++i) mt_random();// erstmal den Twister ordentlich aufwaermen!
}    

void runSweep(void *spins, const double J, const double kB, int N, int dim, int steps, double sweep_init, double sweep_end, double sweep_step, double C, char sweepPar, char sweepMode, char hystRetour, char calcMode, char filmMode);

int main(int argc, char **argv)
{
    int N, dim, steps = 1000;
    double sweep_init = 1.5, sweep_end = 4.0, sweep_step = 0.1, C = 0.0000;
    char sweepPar = 'T', sweepMode = 'n', hystMode = 'n', calcMode = 'n', filmMode = 'n';
    
    getParameters(argc, argv, &N, &dim, &steps, &sweepPar, &sweep_init, &sweep_end, &sweep_step, &C, &sweepMode, &hystMode, &calcMode, &filmMode);

    printf("\nattempting to execute %s %d %d %d %c %f %f %f %f %c %c %c %c \n\n", argv[0], N, dim, steps, sweepPar, sweep_init, sweep_end, sweep_step, C, sweepMode, hystMode, calcMode, filmMode);
    
    initialize(filmMode);
    void *spins;    
    const double J=1.0, kB = 1.0;
    spins=matrixMallocDim(N, dim);
    runSweep(spins, J, kB, N, dim, steps, sweep_init, sweep_end, sweep_step, C, sweepPar, sweepMode, 'n', calcMode, filmMode);
/*
    double S, T, B, dE;
    const double J=1.0, kB = 1.0;
    int i,j, rpos, cpos, zpos=0; //row-position, column-position, pos in third dim(can be anything/is not used in 2d case)
    void *spins;     
    int spinSum = 0, edgeSum = 0;
    char filename[50];
    unsigned long imagecounter = 0, changecounter = 0; // fuer Film     

    spins=matrixMallocDim(N, dim);    
    if(sweepPar == 'T') B = C; //choose the right sweepPar (set constant)
    else if(sweepPar == 'B') T = C;
    
    // a simulation for each sweep increment
    for(S=sweep_init; S<=sweep_end+sweep_step/2; S += sweep_step)
    {
        if(sweepMode == 'n' || S == sweep_init) matrixRandFillDim(spins,N,dim); // fill matrix with random 1 or -1
        spinSum = spinSumDim(spins, N, dim);
        edgeSum = edgeSumDim(spins, N, dim);
        if(sweepPar == 'T') T = S; //choose the right sweepPar (set variable)
        else if(sweepPar == 'B') B = S;
        imagecounter = 0;
        
        for(i=0; i<steps; ++i)
        {           
            for(j=0; j<pow(N,dim); ++j)
            {
                if(filmMode == 'y' && dim ==2 && changecounter%5000 == 0) // falls ein Film erstellt werden soll
                {
                    ++imagecounter;
                    sprintf(filename, "film/data%lu", imagecounter);
                    //sprintf(filename, "film/data_T=%f_B=%f_500changes%lu.txt", T, B, imagecounter);
                    matrixPrint2Dfile(spins,N,N, filename);
                }            
                
                // select random spin, calculate dE, accept spin flip or not
                rpos = mt_random() % N;
                cpos = mt_random() % N;
                if(dim == 3) zpos = mt_random() % N;
                if((dE = calcEnergyDiff(spins, rpos, cpos, zpos, N, J, B, spinSum, dim)) < 0 || 
                    mt_random()/ (double) MT_MAX < exp(-dE/kB/T))
                {
                    edgeSum -= 2*neighSumDim(spins, rpos, cpos, zpos, N, dim);
                    if(dim == 2)
                    {
                        int **spin =  (int **)spins;
                        spinSum -= 2*spin[rpos][cpos];
                        spin[rpos][cpos] *= -1;
                    }
                    if(dim == 3)
                    {
                        int ***spin =  (int ***)spins;
                        spinSum -= 2*spin[rpos][cpos][zpos];
                        spin[rpos][cpos][zpos] *= -1;
                    }
                    if(changecounter%500 == 0)
                    {
                        sprintf(filename, "output/%s_T=%f_B=%f.txt", ENERGYPERMAG, T, B);
                        writeOutputFF(magPerSpinDim(spins, N, dim), calcEnergy(J, B, N, spinSum, edgeSum), filename);
                    }  
                    ++changecounter;
                }
            }              
            sprintf(filename, "output/MagperStep_T=%f_B=%f.txt", T, B);
            writeOutputFF(i, magPerSpinDim(spins, N, dim), filename);
        }
        sprintf(filename, "output/matrix_T=%f_B=%f_end.txt", T, B);
        if(dim == 2) matrixPrint2Dfile(spins,N,N, filename);
        writeOutputFFF(T, magPerSpinDim(spins, N, dim), B, MAGPERSPINOUTPUT);
        if(sweepPar == 'T') printf("=> T = %f finished\n", T);
        else if(sweepPar == 'B') printf("=> B = %f finished\n", B);
    }*/
    if(hystMode == 'y')
    {
        runSweep(spins, J, kB, N, dim, steps, sweep_init, sweep_end, sweep_step, C, sweepPar, sweepMode, 'y', calcMode, filmMode);
        /*
        for(S=sweep_end; S>=sweep_init-sweep_step/2; S -= sweep_step)
        {
            if(sweepMode == 'n' || S == sweep_init) matrixRandFillDim(spins,N,dim); // fill matrix with random 1 or -1
            spinSum = spinSumDim(spins, N, dim);
            edgeSum = edgeSumDim(spins, N, dim);
            if(sweepPar == 'T') T = S; //choose the right sweepPar (set variable)
            else if(sweepPar == 'B') B = S;
            imagecounter = 0;
        
            for(i=0; i<steps; ++i)
            {           
                for(j=0; j<pow(N,dim); ++j)
                {
                    if(filmMode == 'y' && dim == 2 && changecounter%5000 == 0) // falls ein Film erstellt werden soll
                    {
                        ++imagecounter;
                        sprintf(filename, "film/data%lu", imagecounter);
                        //sprintf(filename, "film/data_T=%f_B=%f_500changes%lu.txt", T, B, imagecounter);
                        matrixPrint2Dfile(spins,N,N, filename);
                    }            
                    
                    // select random spin, calculate dE, accept spin flip or not
                    rpos = mt_random() % N;
                    cpos = mt_random() % N;
                    zpos = mt_random() % N;
                    if((dE = calcEnergyDiff(spins, rpos, cpos, zpos, N, J, B, spinSum, dim)) < 0 || 
                        mt_random()/ (double) MT_MAX < exp(-dE/kB/T))
                    {
                        edgeSum -= 2*neighSumDim(spins, rpos, cpos, zpos, N, dim);
                        if(dim == 2)
                        {
                            int **spin =  (int **)spins;
                            spinSum -= 2*spin[rpos][cpos];
                            spin[rpos][cpos] *= -1;
                        }
                        if(dim == 3)
                        {
                            int ***spin =  (int ***)spins;
                            spinSum -= 2*spin[rpos][cpos][zpos];
                            spin[rpos][cpos][zpos] *= -1;
                        }
                        if(changecounter%500 == 0)
                        {
                            sprintf(filename, "output/%s_T=%f_B=%f.txt", ENERGYPERMAG, T, B);
                            writeOutputFF(magPerSpinDim(spins, N, dim), calcEnergy(J, B, N, spinSum, edgeSum), filename);
                        }  
                        ++changecounter;
                    }
                }              
                sprintf(filename, "output/MagperStep_T=%f_B=%f.txt", T, B);
                writeOutputFF(i, magPerSpinDim(spins, N, dim), filename);
            }
            sprintf(filename, "output/matrix_T=%f_B=%f_end.txt", T, B);
            if(dim == 2) matrixPrint2Dfile(spins,N,N, filename);
            writeOutputFFF(T, magPerSpinDim(spins, N, dim), B, MAGPERSPINOUTPUT);
            if(sweepPar == 'T') printf("=> T = %f finished\n", T);
            else if(sweepPar == 'B') printf("=> B = %f finished\n", B);
        }*/
    }
    printf("\n");
    matrixDeleteDim(spins, dim);   
    return EXIT_SUCCESS;
}   


void runSweep(void *spins, const double J, const double kB, int N, int dim, int steps, double sweep_init, double sweep_end, double sweep_step, double C, char sweepPar, char sweepMode, char hystRetour, char calcMode, char filmMode)
{
    double S, T, B, dE;
    double buffer;
    int i,j, rpos, cpos, zpos=0; //row-position, column-position, pos in third dim(can be anything/is not used in 2d case)
    int spinSum = 0, edgeSum = 0;
    char filename[50];
    unsigned long imagecounter = 0, changecounter = 0; // fuer Film
        
    if(sweepPar == 'T') B = C; //choose the right sweepPar (set constant)
    else if(sweepPar == 'B') T = C;
        
    if(hystRetour == 'y')//to use the same for-loop(<=) in both directions
    {
        buffer = -sweep_init;
        sweep_init = -sweep_end;
        sweep_end = buffer;
    }
    /* a simulation for each sweep increment */
    for(S=sweep_init; S<=sweep_end+sweep_step/2; S += sweep_step)
    {
        if(hystRetour == 'y') S = -S;
        if(sweepMode == 'n' || S == sweep_init) matrixRandFillDim(spins,N,dim); // fill matrix with random 1 or -1
        spinSum = spinSumDim(spins, N, dim);
        edgeSum = edgeSumDim(spins, N, dim);
        if(sweepPar == 'T') T = S; //choose the right sweepPar (set variable)
        else if(sweepPar == 'B') B = S;
        imagecounter = 0;
        
        for(i=0; i<steps; ++i)
        {           
            for(j=0; j<pow(N,dim); ++j)
            {
                if(filmMode == 'y' && dim ==2 && changecounter%5000 == 0) // falls ein Film erstellt werden soll
                {
                    ++imagecounter;
                    sprintf(filename, "film/data%lu", imagecounter);
                    //sprintf(filename, "film/data_T=%f_B=%f_500changes%lu.txt", T, B, imagecounter);
                    matrixPrint2Dfile(spins,N,N, filename);
                }            
                
                /* select random spin, calculate dE, accept spin flip or not */
                rpos = mt_random() % N;
                cpos = mt_random() % N;
                if(dim == 3) zpos = mt_random() % N;
                if((dE = calcEnergyDiff(spins, rpos, cpos, zpos, N, J, B, spinSum, dim)) < 0 || 
                    mt_random()/ (double) MT_MAX < exp(-dE/kB/T))
                {
                    edgeSum -= 2*neighSumDim(spins, rpos, cpos, zpos, N, dim);
                    if(dim == 2)
                    {
                        int **spin =  (int **)spins;
                        spinSum -= 2*spin[rpos][cpos];
                        spin[rpos][cpos] *= -1;
                    }
                    if(dim == 3)
                    {
                        int ***spin =  (int ***)spins;
                        spinSum -= 2*spin[rpos][cpos][zpos];
                        spin[rpos][cpos][zpos] *= -1;
                    }
                    if(changecounter%500 == 0)
                    {
                        sprintf(filename, "output/%s_T=%f_B=%f.txt", ENERGYPERMAG, T, B);
                        writeOutputFF(magPerSpinDim(spins, N, dim), calcEnergy(J, B, N, spinSum, edgeSum), filename);
                    }  
                    ++changecounter;
                }
            }              
            sprintf(filename, "output/MagperStep_T=%f_B=%f.txt", T, B);
            writeOutputFF(i, magPerSpinDim(spins, N, dim), filename);
        }
        sprintf(filename, "output/matrix_T=%f_B=%f_end.txt", T, B);
        if(dim == 2) matrixPrint2Dfile(spins,N,N, filename);
        writeOutputFFF(T, magPerSpinDim(spins, N, dim), B, MAGPERSPINOUTPUT);
        if(sweepPar == 'T') printf("=> T = %f finished\n", T);
        else if(sweepPar == 'B') printf("=> B = %f finished\n", B);
        if(hystRetour == 'y') S = -S;
    }    
}


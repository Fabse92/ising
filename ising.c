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

typedef struct
{
    void *spins;
    int N, dim, steps;
    double J, kB;     // kann scheinbar nicht mehr const sein :/
    double sweep_init, sweep_end, sweep_step, C;
    char sweepPar, sweepMode, hystMode, calcMode, filmMode;
  
} Parameters;

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

void getParameters(int argc, char **argv, Parameters *para)
{
    //zunaechst Standardwerte
    para->steps = 1000;
    para->J = 1.0; para->kB = 1.0; 
    para->sweep_init = 1.5; para->sweep_end = 4.0; para->sweep_step = 0.1; para->C = 0.0000;
    para->sweepPar = 'T'; para->sweepMode = 'n'; para->hystMode = 'y'; para->calcMode = 'n'; para->filmMode = 'n';   

    if(argc < 3 
      || sscanf(argv[1], "%d", &para->N) != 1 || para->N < 1  
      || sscanf(argv[2], "%d", &para->dim) != 1 || (para->dim != 2 && para->dim != 3)
      || (argc > 3 && sscanf(argv[3], "%d", &para->steps) != 1) || para->steps < 1
      || (argc > 4 && sscanf(argv[4], "%c", &para->sweepPar) != 1) || (para->sweepPar != 'B' && para->sweepPar != 'T')
      || (argc > 5 && sscanf(argv[5], "%lf", &para->sweep_init) != 1)
      || (argc > 6 && sscanf(argv[6], "%lf", &para->sweep_end) != 1)
      || (argc > 7 && sscanf(argv[7], "%lf", &para->sweep_step) != 1)
      || (argc > 8 && sscanf(argv[8], "%lf", &para->C) != 1)
      || (argc > 9 && sscanf(argv[9], "%c", &para->sweepMode) != 1) || (para->sweepMode != 'y' && para->sweepMode != 'n')
      || (argc > 10 && sscanf(argv[10], "%c", &para->hystMode) != 1) || (para->hystMode != 'y' && para->hystMode != 'n')
      || (argc > 11 && sscanf(argv[11], "%c", &para->calcMode) != 1) || (para->calcMode != 'n' && para->calcMode != 'm')
      || (argc > 12 && sscanf(argv[12], "%c", &para->filmMode) != 1) || (para->filmMode != 'y' && para->filmMode != 'n'))
        usage(argv[0]);    
        
    printf("\nattempting to execute %s %d %d %d %c %f %f %f %f %c %c %c %c \n\n", argv[0], para->N, para->dim, para->steps, para->sweepPar, para->sweep_init, para->sweep_end, para->sweep_step, para->C, para->sweepMode, para->hystMode, para->calcMode, para->filmMode); 
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

void runSweep(Parameters *parameters); //forward declaration


int main(int argc, char **argv)
{
    char hystRetour = 'n';
    Parameters para;
        
    getParameters(argc, argv, &para);    
    initialize(para.filmMode);
    para.spins = matrixMallocDim(para.N, para.dim);
    
    if(para.hystMode == 'y')
    {
        hystRetour = 'y';
        para.hystMode = 'n';
    }
    runSweep(&para);

    if(hystRetour == 'y')
    {
        para.hystMode = 'y';
        runSweep(&para);
    }
    
    printf("\n");
    matrixDeleteDim(para.spins, para.dim);   
    return EXIT_SUCCESS;
}   


void runSweep(Parameters *para)
{
    double S, T, B, dE;
    double buffer;
    int i,j, rpos, cpos, zpos=0; //row-position, column-position, pos in third dim(can be anything/is not used in 2d case)
    int spinSum = 0, edgeSum = 0;
    char filename[50];
    unsigned long imagecounter = 0, changecounter = 0; // fuer Film
        
    if(para->sweepPar == 'T') B = para->C; //choose the right sweepPar (set constant)
    else if(para->sweepPar == 'B') T = para->C;
        
    if(para->hystMode == 'y')//to use the same for-loop(<=) in both directions
    {
        buffer = -para->sweep_init;
        para->sweep_init = -para->sweep_end;
        para->sweep_end = buffer;
    }
    /* a simulation for each sweep increment */
    for(S = para->sweep_init; S <= para->sweep_end+para->sweep_step/2; S += para->sweep_step)
    {
        if(para->hystMode == 'y') S = -S;
        if(para->sweepMode == 'n' || S == para->sweep_init) matrixRandFillDim(para->spins,para->N,para->dim); // fill matrix with random 1 or -1
        spinSum = spinSumDim(para->spins, para->N, para->dim);
        edgeSum = edgeSumDim(para->spins, para->N, para->dim);
        if(para->sweepPar == 'T') T = S; //choose the right sweepPar (set variable)
        else if(para->sweepPar == 'B') B = S;
        imagecounter = 0;
        
        for(i = 0; i < para->steps; ++i)
        {           
            for(j = 0; j < pow(para->N, para->dim); ++j)
            {
                if(para->filmMode == 'y' && para->dim ==2 && changecounter%5000 == 0) // falls ein Film erstellt werden soll
                {
                    ++imagecounter;
                    sprintf(filename, "film/data%lu", imagecounter);
                    //sprintf(filename, "film/data_T=%f_B=%f_500changes%lu.txt", T, B, imagecounter);
                    matrixPrint2Dfile(para->spins,para->N,para->N, filename);
                }            
                
                /* select random spin, calculate dE, accept spin flip or not */
                rpos = mt_random() % para->N;
                cpos = mt_random() % para->N;
                if(para->dim == 3) zpos = mt_random() % para->N;
                if((dE = calcEnergyDiff(para->spins, rpos, cpos, zpos, para->N, para->J, B, spinSum, para->dim)) < 0 || 
                    mt_random()/ (double) MT_MAX < exp(-dE/para->kB/T))
                {
                    edgeSum -= 2*neighSumDim(para->spins, rpos, cpos, zpos, para->N, para->dim);
                    if(para->dim == 2)
                    {
                        int **spin =  (int **)para->spins;
                        spinSum -= 2*spin[rpos][cpos];
                        spin[rpos][cpos] *= -1;
                    }
                    if(para->dim == 3)
                    {
                        int ***spin =  (int ***)para->spins;
                        spinSum -= 2*spin[rpos][cpos][zpos];
                        spin[rpos][cpos][zpos] *= -1;
                    }
                    if(changecounter%500 == 0)
                    {
                        sprintf(filename, "output/%s_T=%f_B=%f.txt", ENERGYPERMAG, T, B);
                        writeOutputFF(magPerSpinDim(para->spins, para->N, para->dim), calcEnergy(para->J, B, para->N, spinSum, edgeSum), filename);
                    }  
                    ++changecounter;
                }
            }              
            sprintf(filename, "output/MagperStep_T=%f_B=%f.txt", T, B);
            writeOutputFF(i, magPerSpinDim(para->spins, para->N, para->dim), filename);
        }
        sprintf(filename, "output/matrix_T=%f_B=%f_end.txt", T, B);
        if(para->dim == 2) matrixPrint2Dfile(para->spins,para->N,para->N, filename);
        writeOutputFFF(T, magPerSpinDim(para->spins, para->N, para->dim), B, MAGPERSPINOUTPUT);
        if(para->sweepPar == 'T') printf("=> T = %f finished\n", T);
        else if(para->sweepPar == 'B') printf("=> B = %f finished\n", B);
        if(para->hystMode == 'y') S = -S;
    }    
}


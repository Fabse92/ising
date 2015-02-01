#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "matrix.h" 
#include "utility.h"
#include "ising.h"


int main(int argc, char **argv)
{
    double buffer;
    Parameters para;
        
    getParameters(argc, argv, &para);    
    initialize(para.filmMode);
    para.spins = matrixMallocDim(para.N, para.dim);
    
    //ausrechnen, welches sweep_end das letzte ist, das in der alten Schleife genutzt wurde
    para.sweep_end = para.sweep_init + para.sweep_step * (int) ((para.sweep_end - para.sweep_init) / para.sweep_step); // schöner, wenn man eine Retour macht
    
    para.sweep_end += para.sweep_step; // damit Schleife nicht eine Iteration zu früh aufhört
    runSweep(&para);
    para.sweep_end -= para.sweep_step; // wieder abziehen

    if(para.hystMode == 'y')
    {        
        buffer = para.sweep_init;
        para.sweep_init = para.sweep_end;
        para.sweep_end = buffer;
        para.sweep_step = -para.sweep_step;
        
        para.sweep_end += para.sweep_step; // damit Schleife nicht zu früh aufhört        
        runSweep(&para);
    }
    
    printf("\n");
    matrixDeleteDim(para.spins, para.dim);   
    return EXIT_SUCCESS;
}   


void runSweep(Parameters *para)
{
    double S, T, B, dE;
    int i,j, rpos, cpos, zpos=0; //row-position, column-position, pos in third dim(can be anything/is not used in 2d case)
    int spinSum = 0, edgeSum = 0;
    char filename[50];
    unsigned long imagecounter = 0, changecounter = 0; // fuer Film
        
    if(para->sweepPar == 'T') B = para->C; //choose the right sweepPar (set constant)
    else if(para->sweepPar == 'B') T = para->C;
        
    /* a simulation for each sweep increment */
    for(S = para->sweep_init; fabs(para->sweep_step) <= fabs(S - para->sweep_end); S += para->sweep_step) // fabs ist der Absolutbetrag eines double,   warum war hier para->sweep_end+para->sweep_step/2 als Bedingung?
    {
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
    }    
}


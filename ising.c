#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "matrix.h" 
#include "utility.h"
#include "ising.h"
#include "cluster.h"


int main(int argc, char **argv)
{
    double buffer;
    Parameters para;
        
    getParameters(argc, argv, &para);    
    initialize(&para);
    para.spins = matrixMallocDim(para.N, para.dim);
    if (para.clusterMode == 'y') initializeClusterVars(&para);
    
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
    if (para.clusterMode == 'y') clusterDelete(&para);
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
        
    if(para->sweepPar == 'T') 
    {
        B = para->C; //choose the right sweepPar (set constant)
        
    } else if(para->sweepPar == 'B') 
    {
        T = para->C;
        if(para->clusterMode == 'y') 
            para->cv.addProbability = 1 - exp(-2*para->J/T);
    }
        
    /* a simulation for each sweep increment */
    for(S = para->sweep_init; fabs(para->sweep_step) / 2 < fabs(S - para->sweep_end); S += para->sweep_step)
    {
        if(para->sweepPar == 'T') 
        {
            T = S;                       //choose the right sweepPar (set variable)
            if(para->clusterMode == 'y') // addProbability has to be updated with Temperatur
                para->cv.addProbability = 1 - exp(-2*para->J/T);
        }
        else if(para->sweepPar == 'B') B = S;
        
        if(para->sweepMode == 'n' || S == para->sweep_init) 
        {
          matrixFillDim(para->spins,para->N,para->dim,para->matMode); // fill matrix with random +-1, +1 or -1
        }
        spinSum = spinSumDim(para->spins, para->N, para->dim);
        edgeSum = edgeSumDim(para->spins, para->N, para->dim);
        
        imagecounter = 0;
        
        for(i = 0; i < para->steps; ++i)
        {           
            for(j = 0; j < pow(para->N, para->dim); ++j)
            {
                if(para->filmMode == 'y' && para->dim ==2 && changecounter%5000 == 0) // falls ein Film erstellt werden soll
                {
                    ++imagecounter;
                    sprintf(filename, "film/data%lu", imagecounter);
                    matrixPrint2Dfile(para->spins,para->N,para->N, filename);
                }            
                
                /* select random spin, calculate dE, accept spin flip or not */
                rpos = mt_random() % para->N;
                cpos = mt_random() % para->N;
                if(para->dim == 3) zpos = mt_random() % para->N;
                
                if(para->clusterMode == 'y')
                {
                    resetClusterMatrix(para);
                    growCluster(para, rpos, cpos, ((int **)para->spins)[rpos][cpos]);
                    if(para->filmMode == 'y')
                    {
                        changecounter += 5000;
                    }
                    spinSum = spinSumDim(para->spins, para->N, para->dim);
                    edgeSum = edgeSumDim(para->spins, para->N, para->dim);
                    break; // ein Montecarlo Schritt wurde damit ausgeführt                    
                
                } else if((dE = calcEnergyDiff(para->spins, rpos, cpos, zpos, para->N, para->J, B, para->dim)) < 0 || mt_random()/ (double) MT_MAX < exp(-dE/para->kB/T))
                {
                    ++changecounter;
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
                }
            }              
            sprintf(filename, "%s_T=%f_B=%f.txt", ENERGYPERMAG, T, B);
            writeOutputFFF(i, magPerSpinDim(para->spins, para->N, para->dim), calcEnergy(para->J, B, spinSum, edgeSum), filename);
        }
        if(para->dim == 2)
        {
            sprintf(filename, "output/matrix_T=%f_B=%f_end.txt", T, B);
            matrixPrint2Dfile(para->spins,para->N,para->N, filename);
        }
        writeOutputFFF(T, magPerSpinDim(para->spins, para->N, para->dim), B, MAGPERSPINOUTPUT);
        if(para->sweepPar == 'T') printf("=> T = %f finished\n", T);
        else if(para->sweepPar == 'B') printf("=> B = %f finished\n", B);
    }    
}


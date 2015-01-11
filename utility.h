#include <assert.h>

/** Calculates the energy of the 2-dimensional spin-matrix */
double calcEnergy2D(int **spins)
{
    return 0.0;
}

/** Calculates the difference in energy that a spin flip at position (row,col) would cause (in units of J!)
(it has to be a square 2D matrix) (periodic boundary conditions)*/
/* 2 * Summe ueber Nachbarn * eigener spin */
int calcEnergyDiff2D(int **spins, int row, int col, int len)
{
    assert(spins != NULL && len > 0);
    assert(row >= 0 && col >= 0);
    
    int neighSum, spin, rspin, lspin, uspin, dspin;
    
    spin = spins[row][col];
    if(col-1 < 0) lspin = spins[row][len-1];
    else lspin = spins[row][col-1];
    if(col+1 > len-1) rspin = spins[row][0];
    else rspin = spins[row][col+1];
    if(row-1 < 0) uspin = spins[len-1][col];
    else uspin = spins[row-1][col];
    if(row+1 > len-1) dspin = spins[0][col];
    else dspin = spins[row+1][col];
    
    neighSum = rspin + lspin + uspin + dspin;
    return 2*neighSum*spin; 
}

/** Calculates the difference in energy that a spin flip at position (x1, x2, x3) would cause (in units of J)
(it has to be a cubic 3D matrix) (periodic boundary conditions) */
int calcEnergyDiff3D(int ***spins, int x1, int x2, int x3, int len)
{
    assert(spins!=NULL && len > 0);
    assert(x1>=0 && x2>=0 && x3>=0);
    
    int i, neighSum=0, spin;
    int nspin[6]; //6 neighbor spins
    
    spin = spins[x1][x2][x3];
    
    if(x1-1<0) nspin[0] = spins[len-1][x2][x3];
    else nspin[0] = spins[x1-1][x2][x3];
    if(x1+1>len-1) nspin[1] = spins[0][x2][x3];
    else nspin[1] = spins[x1+1][x2][x3];
    if(x2-1<0) nspin[2] = spins[x1][len-1][x3];
    else nspin[2] = spins[x1][x2-1][x3];
    if(x2+1>len-1) nspin[3] = spins[x1][0][x3];
    else nspin[3] = spins[x1][x2+1][x3];
    if(x3-1<0) nspin[4] = spins[x1][x2][len-1];
    else nspin[4] = spins[x1][x2][x3-1];
    if(x3+1>len-1) nspin[5] = spins[x1][x2][0];
    else nspin[5] = spins[x1][x2][x3+1];

    for(i=0; i<6; ++i)
    {
        neighSum += nspin[i];
    }
    return 2*neighSum*spin;
}


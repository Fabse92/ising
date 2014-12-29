#include <assert.h>

/** Calculates the energy of the 2-dimensional spin-matrix */
double calcEnergy(int **spins)
{
    return 0.0;
}

/** Calculates the difference in energy that a spin flip at position (row,col) would cause (in units of J!) */
/* 2 * Summe ueber Nachbarn * eigener spin */
int calcEnergyDiff(int **spins, int row, int col, int dim)
{
    assert(spins != NULL && dim > 0);
    assert(row >= 0 && col >= 0);
    
    int neighSum, spin, rspin, lspin, uspin, dspin;
    
    spin = spins[row][col];
    if(col-1 < 0) lspin = spins[row][dim-1];
    else lspin = spins[row][col-1];
    if(col+1 > dim-1) rspin = spins[row][0];
    else rspin = spins[row][col+1];
    if(row-1 < 0) uspin = spins[dim-1][col];
    else uspin = spins[row-1][col];
    if(row+1 > dim-1) dspin = spins[0][col];
    else dspin = spins[row+1][col];
    
    neighSum = rspin + lspin + uspin + dspin;
    return 2*neighSum*spin; 
}

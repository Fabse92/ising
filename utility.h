#ifndef UTILITY_H
#define UTILITY_H

#include <assert.h>
#include <stdlib.h>

/** Appends "value1" and "value2" to the File "filename" 
    Where "value1" and "value2" are of type double
*/
void writeOutputFF(double value1, double value2, const char* filename)
{
    FILE *fp = fopen(filename, "a");
    char output[50];
    
    if (fp == NULL)
        fprintf(stderr, "Konnte nicht in Datei %s schreiben \n", filename);
    
    sprintf(output, "%f \t %f \n", value1, value2);
    fputs( output, fp);        
        
    fclose (fp);
}

/** Appends "value1", "value2" and "value3" to the File "filename" 
    Where "value1", "value2" and "value3" are of type double
*/
void writeOutputFFF(double value1, double value2, double value3, const char* filename)
{
    FILE *fp = fopen(filename, "a");
    char output[50];
    
    if (fp == NULL)
        fprintf(stderr, "Konnte nicht in Datei %s schreiben \n", filename);
    
    sprintf(output, "%f \t %f \t %f \n", value1, value2, value3);
    fputs( output, fp);        
        
    fclose (fp);
}

/** Calculates the Energy of the current spin lattice in calcMode and Dimension dim*/
double calcEnergy(double J, double B, int len, int spinSum, int edgeSum, char calcMode, int dim)
{
  assert(dim == 2 || dim == 3);
  double energy;
  
  if(calcMode == 'n') energy = -J*edgeSum/2 - B*spinSum;
  else if (calcMode == 'm' && dim == 2) energy = -(4*J*spinSum/(len*len) + B)*spinSum;
  else if (calcMode == 'm' && dim == 3) energy = -(6*J*spinSum/(len*len*len) + B)*spinSum;

  return energy;
}

/** Calculates the Energy of the current spin lattice with nearest neighbor Hamiltonian*/
/* you can use this only if you have a current sum of all spins and sum over all edges */
double calcEnergyNN(double J, double B, int len, int spinSum, int edgeSum)
{
  double energy;
  
  energy = -J*edgeSum/2 - B*spinSum;

  return energy;
}

/** Calculates the Energy of the current spin lattice with Mean-Field Hamiltonian*/
double calcEnergyMFT(double J, double B, int len, int spinSum, int dim)
{
  assert(dim == 2 || dim == 3);
  double energy;
  
  if(dim == 2) energy = -(4*J*spinSum/(len*len) + B)*spinSum;
  else if(dim == 3) energy = -(6*J*spinSum/(len*len) + B)*spinSum;
  
  return energy;
}

/** Calculates the sum of the edges of one spin in 2 dimensions */
/* example: edge 1....-1 adds -1 to the sum; edge -1....-1 adds +1 */
int neighSum2D(int **spins, int x1, int x2, int len)
{
    assert(spins!=NULL && len > 0);
    assert(x1>=0 && x2>=0);
    
    int i, neighSum=0, spin;
    int nspin[4]; //4 neighbor spins
    
    spin = spins[x1][x2];
    
    if(x1-1<0) nspin[0] = spins[len-1][x2];
    else nspin[0] = spins[x1-1][x2];
    if(x1+1>len-1) nspin[1] = spins[0][x2];
    else nspin[1] = spins[x1+1][x2];
    if(x2-1<0) nspin[2] = spins[x1][len-1];
    else nspin[2] = spins[x1][x2-1];
    if(x2+1>len-1) nspin[3] = spins[x1][0];
    else nspin[3] = spins[x1][x2+1];

    for(i=0; i<4; ++i)
    {
        neighSum += nspin[i];
    }
    return neighSum*spin;
}

/** Calculates the sum of the edges of one spin in 3 dimensions */
/* example: edge 1....-1 adds -1 to the sum; edge -1....-1 adds +1 */
int neighSum3D(int ***spins, int x1, int x2, int x3, int len)
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
    return neighSum*spin;
}


/** Calculates the sum of the edges of one spin in dim dimensions */
/* example: edge 1....-1 adds -1 to the sum; edge -1....-1 adds +1 */
int neighSumDim(void *spins, int x1, int x2, int x3, int len, int dim)
{
    assert(dim == 2 || dim == 3);
    assert(spins!=NULL && len>0);
    
    int result;
    if(dim == 2)
    {
        result = neighSum2D((int **)spins, x1, x2, len);
    }
    else if(dim == 3)
    {
        result = neighSum3D((int ***)spins, x1, x2, x3, len);
    }
    return result;    
}

/** Calculates the sum of all spins in dim dimensions */
int spinSumDim(void *spins, int len, int dim)
{
    assert(dim == 2 || dim == 3);
    
    int sum=0;
    int i,j;
    if(dim == 2)
    {
        int **spin =  (int **)spins;        
        for (i=0;i<len;++i)
        {
            for(j=0;j<len;++j)
            {
               sum = sum + spin[i][j];
            }
        }
    }
    else if (dim == 3)
    {
        int k;
        int ***spin = (int ***)spins;        
        for (i=0;i<len;++i)
        {
            for(j=0;j<len;++j)
            {
                for(k=0;k<len;++k)
                {
                    sum = sum + spin[i][j][k];
                }
            }
        }
    }
    return sum;
}

/** Calculates the sum of all edges in the lattice in dim dimensions (for NN-Energy calculating)*/
int edgeSumDim(void *spins, int len, int dim)
{
    assert(dim == 2 || dim == 3);
    
    int sum=0;
    int i,j;
    if(dim == 2)
    {
        for (i=0;i<len;++i)
        {
            for(j=0;j<len;++j)
            {
               sum = sum + neighSum2D(spins, i, j, len);
            }
        }
    }
    else if (dim == 3)
    {
        int k;    
        for (i=0;i<len;++i)
        {
            for(j=0;j<len;++j)
            {
                for(k=0;k<len;++k)
                {
                    sum = sum + neighSum3D(spins, i, j, k, len);
                }
            }
        }
    }
    return sum/2;
}

/** Calculates the difference in energy that a spin flip at position (row,col) would cause (in units of J!)
(it has to be a square 2D matrix) (periodic boundary conditions)*/
double calcNNEnergyDiff2D(int **spins, int row, int col, int len, double J, double B)
{
    assert(spins != NULL && len > 0);
    assert(row >= 0 && col >= 0);
    
    int neighTerm, spin;
    spin = spins[row][col];
    neighTerm = neighSum2D(spins, row, col, len);
    return 2*(J*neighTerm + B*spin);
}

/** Calculates the difference in energy that a spin flip at position (x1, x2, x3) would cause
(it has to be a cubic 3D matrix) (periodic boundary conditions) */
double calcNNEnergyDiff3D(int ***spins, int x1, int x2, int x3, int len, double J, double B)
{
    assert(spins!=NULL && len > 0);
    assert(x1>=0 && x2>=0 && x3>=0);
    
    int neighTerm, spin;
    spin = spins[x1][x2][x3];
    neighTerm = neighSum3D(spins, x1, x2, x3, len);
    return 2*(J*neighTerm + B*spin);
}

/**  Calculates the difference in energy in mean field approximation that a spin flip at position (row,col) 
would cause (it has to be a square 2D matrix) */
double calcMFTEnergyDiff2D(int **spins, int row, int col, int len, int sum, double J, double B)
{
  assert(spins != NULL && len > 0);
  assert(row >= 0 && col >= 0);
  
  int spin, mfield;
  double diff;
  
  spin = spins[row][col];
  mfield = sum - spin;
  diff = (8*J/(len*len)*mfield + B)*2*spin;
    
  return diff;
}

/**  Calculates the difference in energy in mean field approximation that a spin flip at position (row,col) 
would cause (it has to be a cubic 3D matrix) */
double calcMFTEnergyDiff3D(int ***spins, int row, int col, int depth, int len, int sum, double J, double B)
{
    assert(spins != NULL && len > 0);
    assert(row >= 0 && col >= 0 && depth >= 0);
    
    int spin, mfield;
    double diff;
    
    spin = spins[row][col][depth];
    mfield = sum - spin;
    diff = (12*J/(len*len*len)*mfield + B)*2*spin;
      
    return diff;
}

/** Calculates magnetisation per spin for a 2D- or 3D-matrix */
double magPerSpinDim(void *spins, int len, int dim)
{
    assert(dim == 2 || dim == 3);
    assert(spins!=NULL && len>0);
    
    double result;
    int sum = spinSumDim(spins, len, dim);
    if(dim == 2)
    {
        result = ((double)sum / (len*len));
    }
    else if(dim == 3)
    {
        result = ((double)sum / (len*len*len));
    }
    return result;
}

/** Calculates the difference in energy that a spin flip at pos (x1,x2(,x3)) would cause */
/* x3 is not used in 2-dimensional case */
double calcNNEnergyDiffDim(void *spins, int x1, int x2, int x3, int len, double J, double B, int dim)
{
    assert(dim == 2 || dim == 3);
    assert(spins!=NULL && len>0);
    
    double result;
    if(dim == 2)
    {
        result = calcNNEnergyDiff2D((int **)spins, x1, x2, len, J, B);
    }
    else if(dim == 3)
    {
        result = calcNNEnergyDiff3D((int ***)spins, x1, x2, x3, len, J, B);
    }
    return result;
}

/** Calculates the difference in energy in mean field approximation that a spin flip at pos (x1,x2(,x3)) would cause */
/* x3 is not used in 2-dimensional case */
double calcMFTEnergyDiffDim(void *spins, int x1, int x2, int x3, int len, int sum, double J, double B, int dim)
{
    assert(dim == 2 || dim == 3);
    assert(spins!=NULL && len>0);
    
    double result;
    if(dim == 2)
    {
        result = calcMFTEnergyDiff2D((int **)spins, x1, x2, len, sum, J, B);
    }
    else if(dim == 3)
    {
        result = calcMFTEnergyDiff3D((int ***)spins, x1, x2, x3, len, sum, J, B);
    }
    return result;
}

/** Calculates the difference in energy that a spin flip at its position would cause (periodic boundary conditions)*/
double calcEnergyDiff(void *spins, int x1, int x2, int x3, int len, double J, double B, int sum, int dim, char calcMode)
{
  assert(calcMode == 'n' || calcMode == 'm');
  assert(dim == 2 || dim == 3);
  assert(spins!=NULL && len>0);
  double dE;
  
  if(calcMode == 'n' && dim == 2) dE = calcNNEnergyDiff2D((int **)spins, x1, x2, len, J, B);
  else if(calcMode == 'n' && dim == 3) dE = calcNNEnergyDiff3D((int ***)spins, x1, x2, x3, len, J, B);
  else if(calcMode == 'm' && dim == 2) dE = calcMFTEnergyDiff2D((int **)spins, x1, x2, len, sum, J, B);
  else if(calcMode == 'm' && dim == 3) dE = calcMFTEnergyDiff3D((int ***)spins, x1, x2, x3, len, sum, J, B);
  
  return dE;
}

/* This program implements the Mersenne twister algorithm for generation of pseudorandom numbers. 
The program returns random integers in the range 0 to 2^32-1 (this holds even if a long int is
larger than 32 bits). Timing with gcc indicates that it is about twice as fast as the built in 
rand function. The original code was written by Michael Brundage and has been placed in the 
public domain. There are a three minor changes here: 
(1) This comment has been added to the program.
(2) Type specifiers (ul) have been appended to constants.
(3) A commented out block near the end has been removed. */

#define MT_LEN 624

int mt_index;
unsigned long mt_buffer[MT_LEN];

/* Function has to be called ones, before using mt_random() */

void mt_init() {
    int i;
    for (i = 0; i < MT_LEN; i++)
        mt_buffer[i] = rand();
    mt_index = 0;
}

#define MT_IA           397
#define MT_IB           (MT_LEN - MT_IA)
#define UPPER_MASK      0x80000000
#define LOWER_MASK      0x7FFFFFFF
#define MATRIX_A        0x9908B0DF
#define TWIST(b,i,j)    ((b)[i] & UPPER_MASK) | ((b)[j] & LOWER_MASK)
#define MAGIC(s)        (((s)&1)*MATRIX_A)

unsigned long mt_random() {
    unsigned long * b = mt_buffer;
    int idx = mt_index;
    unsigned long s;
    int i;
	
    if (idx == MT_LEN*sizeof(unsigned long))
    {
        idx = 0;
        i = 0;
        for (; i < MT_IB; i++) {
            s = TWIST(b, i, i+1);
            b[i] = b[i + MT_IA] ^ (s >> 1) ^ MAGIC(s);
        }
        for (; i < MT_LEN-1; i++) {
            s = TWIST(b, i, i+1);
            b[i] = b[i - MT_IB] ^ (s >> 1) ^ MAGIC(s);
        }
        
        s = TWIST(b, MT_LEN-1, 0);
        b[MT_LEN-1] = b[MT_IA-1] ^ (s >> 1) ^ MAGIC(s);
    }
    mt_index = idx + sizeof(unsigned long);
    return *(unsigned long *)((unsigned char *)b + idx);
    /* Here there is a commented out block in MB's original program */
}

#endif


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

/** Calculates the Energy of the current spin lattice of Dimension dim*/
double calcEnergy(double J, double B, int spinSum, int edgeSum)
{
  double energy;
  
  energy = -J*edgeSum/2 - B*spinSum;

  return energy;
}

/** Calculates the sum of the edges of one spin in dim dimensions */
/* example: edge 1....-1 adds -1 to the sum; edge -1....-1 adds +1 */
int neighSumDim(void *spins, int x1, int x2, int x3, int len, int dim)
{
    assert(dim == 2 || dim == 3);
    assert(spins!=NULL && len>0);
    assert(x1>=0 && x2>=0 && x3>=0);
    
    int i, neighSum=0, ispin;
    int nspin[2*dim]; //(2*dim) neighbor spins
    if(dim == 2)
    {
        int **spin =  (int **)spins;
        ispin = spin[x1][x2];
    
        if(x1-1<0) nspin[0] = spin[len-1][x2];
        else nspin[0] = spin[x1-1][x2];
        if(x1+1>len-1) nspin[1] = spin[0][x2];
        else nspin[1] = spin[x1+1][x2];
        if(x2-1<0) nspin[2] = spin[x1][len-1];
        else nspin[2] = spin[x1][x2-1];
        if(x2+1>len-1) nspin[3] = spin[x1][0];
        else nspin[3] = spin[x1][x2+1];
    }
    else if(dim == 3)
    {
        int ***spin = (int ***)spins;
        ispin = spin[x1][x2][x3];
    
        if(x1-1<0) nspin[0] = spin[len-1][x2][x3];
        else nspin[0] = spin[x1-1][x2][x3];
        if(x1+1>len-1) nspin[1] = spin[0][x2][x3];
        else nspin[1] = spin[x1+1][x2][x3];
        if(x2-1<0) nspin[2] = spin[x1][len-1][x3];
        else nspin[2] = spin[x1][x2-1][x3];
        if(x2+1>len-1) nspin[3] = spin[x1][0][x3];
        else nspin[3] = spin[x1][x2+1][x3];
        if(x3-1<0) nspin[4] = spin[x1][x2][len-1];
        else nspin[4] = spin[x1][x2][x3-1];
        if(x3+1>len-1) nspin[5] = spin[x1][x2][0];
        else nspin[5] = spin[x1][x2][x3+1];
    }
    
    for(i=0; i<2*dim; ++i)
    {
        neighSum += nspin[i];
    }
    return neighSum*ispin;    
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
    int i,j,k;
    
    for (i=0;i<len;++i)
    {
        for(j=0;j<len;++j)
        {
            if(dim == 2)
            {
                sum = sum + neighSumDim(spins, i, j, 0, len, 2);
            }
            else if (dim == 3)
            {
                for(k=0;k<len;++k)
                {
                    sum = sum + neighSumDim(spins, i, j, k, len, 3);
                }
            }
        }
    }
    return sum/2;
}

/** Calculates the difference in energy that a spin flip at position (row,col) would cause (in units of J!)
(it has to be a square 2D matrix) (periodic boundary conditions)*/
double calcEnergyDiff2D(int **spins, int row, int col, int len, double J, double B)
{
    assert(spins != NULL && len > 0);
    assert(row >= 0 && col >= 0);
    
    int neighTerm, spin;
    spin = spins[row][col];
    neighTerm = neighSumDim(spins, row, col, 0, len, 2);
    return 2*(J*neighTerm + B*spin);
}

/** Calculates the difference in energy that a spin flip at position (x1, x2, x3) would cause
(it has to be a cubic 3D matrix) (periodic boundary conditions) */
double calcEnergyDiff3D(int ***spins, int x1, int x2, int x3, int len, double J, double B)
{
    assert(spins!=NULL && len > 0);
    assert(x1>=0 && x2>=0 && x3>=0);
    
    int neighTerm, spin;
    spin = spins[x1][x2][x3];
    neighTerm = neighSumDim(spins, x1, x2, x3, len, 3);
    return 2*(J*neighTerm + B*spin);
}

/** Calculates the difference in energy that a spin flip at its position would cause (periodic boundary conditions)*/
double calcEnergyDiff(void *spins, int x1, int x2, int x3, int len, double J, double B, int dim)
{
  assert(dim == 2 || dim == 3);
  assert(spins!=NULL && len>0);
  assert(x1 >= 0 && x2 >= 0 && x3 >= 0);
  
  double dE = 0.0;
  
  if(dim == 2) dE = calcEnergyDiff2D((int **)spins, x1, x2, len, J, B);
  else if(dim == 3) dE = calcEnergyDiff3D((int ***)spins, x1, x2, x3, len, J, B);
    
  return dE;
}

/** Calculates magnetisation per spin for a 2D- or 3D-matrix */
double magPerSpinDim(void *spins, int len, int dim)
{
    assert(dim == 2 || dim == 3);
    assert(spins!=NULL && len>0);
    
    double result = 0.0;
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


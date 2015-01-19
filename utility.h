#ifndef UTILITY_H
#define UTILITY_H

#include <assert.h>
#include <stdlib.h>

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

/** Calculates the energy of the 2-dimensional spin-matrix */
double calcEnergy2D(int **spins)
{
    return 0.0;
}

/** Calculetes the sum of all spins in the 2D-square-lattice */
int spinSum2DSquare(int **spins, int len)
{
  int sum=0;
  int i,j;
  
  for (i=0;i<len;++i)
  {
      for(j=0;j<len;++j)
      {
         sum = sum + spins[i][j];
      }
  }
  return sum;
}

/**  Calculates the difference in energy in mean field approximation that a spin flip at position (row,col) 
would cause (it has to be a square 2D matrix) */
double calcMFTEnergyDiff2DSquare(int **spins, int row, int col, int len, int sum, double J, double B)
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

/** Calculetes the sum of all spins in the 3D-cubic-lattice */
int spinSum3DCubic(int ***spins, int len)
{
    int sum=0;
    int i,j,k;
    
    for (i=0;i<len;++i)
    {
        for(j=0;j<len;++j)
        {
            for(k=0;k<len;++k)
            {
                sum = sum + spins[i][j][k];
            }
        }
    }
    return sum;
}

/**  Calculates the difference in energy in mean field approximation that a spin flip at position (row,col) 
would cause (it has to be a cubic 3D matrix) */
double calcMFTEnergyDiff3DCubic(int ***spins, int row, int col, int depth, int len, int sum, double J, double B)
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

/** Calculates the difference in energy that a spin flip at position (row,col) would cause (in units of J!)
(it has to be a square 2D matrix) (periodic boundary conditions)*/
/* 2 * Summe ueber Nachbarn * eigener spin */
double calcEnergyDiff2DSquare(int **spins, int row, int col, int len, double J, double B)
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
    return (J*neighSum + B)*2*spin;
}

/** Calculates the difference in energy that a spin flip at position (x1, x2, x3) would cause
(it has to be a cubic 3D matrix) (periodic boundary conditions) */
double calcEnergyDiff3DSquare(int ***spins, int x1, int x2, int x3, int len, double J, double B)
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
    return (J*neighSum + B)*2*spin;
}

//** Calculates magnetisation / spin for a 2D-Matrix */
double calcMagperSpin2D(int **spins, int len)
{
    assert(spins!=NULL && len > 0);
    
    int sum = spinSum2DSquare(spins, len);
    return ( (double)sum / ( len * len) );
}

//** Calculates magnetisation / spin for a 3D-Matrix */
double calcMagperSpin3D(int ***spins, int len)
{
    assert(spins!=NULL && len > 0);
    
    int sum = spinSum3DCubic(spins, len);
    return ( (double)sum / ( len * len * len) );
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


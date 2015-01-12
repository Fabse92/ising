#ifndef MATRIX_H
#define MATRIX_H

/* functions for doing things with the matrix */

#include <assert.h>
#include <stdlib.h>
#include "utility.h" // fuer mt_random()

/** allocates a 2-dimensional matrix (rows x cols) */
int **matrixMalloc2D(int rows, int cols)
{
    assert(rows > 0 && cols > 0);
    
    int i, *data, **pointers;
    
    data = malloc(rows*cols*sizeof(*data));
    pointers = malloc(rows*sizeof(*pointers));
    
    for (i=0; i<rows; ++i)
    {
        pointers[i] = data + i*cols;
    }
    
    return pointers;
}

/** deletes a 2-dimensional matrix */
void matrixDelete2D(int **mat)
{
    if(mat!=NULL && *mat!=NULL)
    {
        free(*mat);
        free(mat);    
    }
}

/** prints a 2-dimensional matrix */
void matrixPrint2D(int **mat, int rows, int cols)
{
    int i,j;
    
    for(i=0; i<rows; ++i)
    {
        for(j=0; j<cols; ++j)
        {
            printf("%02d ", mat[i][j]);
        }
        printf("\n");
    }
}

/** Fills a 2-dimensional matrix with random 1 or -1. */
void matrixRandFill2D(int **mat,int rows, int cols)
{
    int i,j;
    for (i=0; i<rows; ++i)
        for (j=0; j<cols; ++j)
        {
            mat[i][j] = 1-2*(mt_random()%2);
        }                
}

/** allocates a 3-dimensional matrix (n1 x n2 x n3) */
int ***matrixMalloc3D(int n1, int n2, int n3)
{
    assert(n1>0 && n2>0 && n3>0);
    
    int i, *data, **pointers1, ***pointers2;
    
    data = malloc(n1*n2*n3*sizeof(*data));
    pointers1 = malloc(n1*n2*sizeof(*pointers1));
    pointers2 = malloc(n1*sizeof(*pointers2));
    assert(pointers2 != NULL);
    
    for(i=0; i<n1*n2; ++i)
    {
        pointers1[i] = data + i*n3;
    }
    for(i=0; i<n1; ++i)
    {
        pointers2[i] = pointers1 + i*n2;
    } 
    return pointers2;
}

/** deletes a 3-dimensional matrix */
void matrixDelete3D(int ***mat)
{
    if(mat!=NULL && *mat!=NULL && **mat!=NULL)
    {
        free(**mat);
        free(*mat);
        free(mat);
    }
}

/** prints a 3-dimensional matrix 
(not very pretty, the 2D-matrices are just printed one after the other, seperated by a '.' ..) */
void matrixPrint3D(int ***mat, int n1, int n2, int n3)
{
    int i,j,k;
    
    for(i=0; i<n1; ++i)
        for(j=0; j<n2; ++j)
        {
            for(k=0; k<n3; ++k)
            {
                printf("%02d ", mat[i][j][k]);
            }
            if(j == n2-1) printf(".");
            printf("\n");
        }
}

/** Fills a 3-dimensional matrix with random -1 or 1 */
void matrixRandFill3D(int ***mat, int n1, int n2, int n3)
{
    int i,j,k;
    
    for(i=0; i<n1; ++i)
        for(j=0; j<n2; ++j)
            for(k=0; k<n3; ++k)
            {
                mat[i][j][k] = 1-2*(mt_random()%2);
            }
}

#endif


#ifndef MATRIX_H
#define MATRIX_H

/* functions for doing things with matrices */

#include <assert.h>
#include <stdlib.h>
#include <stdbool.h>
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

/** allocates a 2-dimensional matrix (rows x cols) for booleans */
bool **matrixMalloc2DBool(int rows, int cols)
{
    assert(rows > 0 && cols > 0);
    
    int i;
    bool *data, **pointers;
    
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

/** prints a 2-dimensional matrix to file */
void matrixPrint2Dfile(int **mat, int rows, int cols, const char* filename)
{
    int i,j;
    FILE *fp = fopen(filename, "w");
    char output[3 * cols + 10];
    sprintf(output, "");
    
    if (fp == NULL)
        fprintf(stderr, "Konnte nicht in Datei %s schreiben \n", filename);
    
    for(i=0; i<rows; ++i)
    {
        for(j=0; j<cols; ++j)
        {
            sprintf(output, "%s%02d ", output, mat[i][j]);
        }
        strcat(output, "\n");
        fputs( output, fp);  
        sprintf(output, "");      
    }
    fclose (fp);
}

/** Fills a 2-dimensional matrix with random 1 or -1. */
void matrixFill2D(int **mat,int rows, int cols, char matMode)
{
    int i,j;
    for (i=0; i<rows; ++i)
        for (j=0; j<cols; ++j)
        {
            if(matMode == 'r') mat[i][j] = 1-2*(mt_random()%2);
            else if(matMode == 'p') mat[i][j] = 1;
            else if(matMode == 'n') mat[i][j] = -1;
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

/** allocates a 3-dimensional matrix (n1 x n2 x n3) for booleans*/
bool ***matrixMalloc3DBool(int n1, int n2, int n3)
{
    assert(n1>0 && n2>0 && n3>0);
    
    int i;
    bool *data, **pointers1, ***pointers2;
    
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
void matrixFill3D(int ***mat, int n1, int n2, int n3, char matMode)
{
    int i,j,k;
    
    for(i=0; i<n1; ++i)
        for(j=0; j<n2; ++j)
            for(k=0; k<n3; ++k)
            {
                if(matMode == 'r') mat[i][j][k] = 1-2*(mt_random()%2);
                else if(matMode == 'p') mat[i][j][k] = 1;
                else if(matMode == 'n') mat[i][j][k] = -1;
            }
}

/** allocates a dim-dimensional matrix ( dim=2,3; NxN(xN) ) */
void *matrixMallocDim(int N, int dim)
{
    assert(N>0);
    assert(dim==2 || dim==3);
    
    void *result = NULL;
    if(dim==2)
    {
        result = (void *)matrixMalloc2D(N,N);
    }
    else if (dim==3)
    {
        result = (void *)matrixMalloc3D(N,N,N);
    }
    return result;
}

/** deletes a dim-dimensional matrix */
void matrixDeleteDim(void *mat, int dim)
{
    assert(dim==2 || dim==3);
    
    if(dim==2)
    {
        matrixDelete2D((int **)mat);
    }
    else if(dim==3)
    {
        matrixDelete3D((int ***)mat);
    }
}

/** prints a dim-dimensional matrix (NxN(xN)) to the console */
void matrixPrintDim(void *mat, int N, int dim)
{
    assert(dim==2 || dim==3);
    
    if(dim==2)
    {
        matrixPrint2D((int **)mat, N, N);
    }
    else if(dim==3)
    {
        matrixPrint3D((int ***)mat, N, N, N);
    }
}

/** fills a dim-dimensional matrix (NxN(xN)) with random -1 or 1 */
void matrixFillDim(void *mat, int N, int dim, char matMode)
{
    assert(dim==2 || dim==3);
    
    if(dim==2)
    {
        matrixFill2D((int **)mat, N, N, matMode);
    }
    else if(dim==3)
    {
        matrixFill3D((int ***)mat, N, N, N, matMode);
    }
}

#endif


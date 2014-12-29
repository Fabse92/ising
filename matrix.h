#ifndef MATRIX_H
#define MATRIX_H

/* functions for doing things with the matrix */

#include <assert.h>
#include <stdlib.h>

/** allocates a 2-dimensional matrix (rows x cols) */
int **matrixMalloc(int rows, int cols)
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
void matrixDelete(int **mat)
{
    assert(mat != NULL && *mat != NULL);

    free(*mat);
    free(mat);    
}

/** prints a 2-dimensional matrix */
void matrixPrint(int **mat, int rows, int cols)
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
void matrixRandFill(int **mat,int rows, int cols)
{
    int i,j;
    for (i=0; i<rows; ++i)
        for (j=0; j<cols; ++j)
        {
            mat[i][j] = 1-2*(rand()%2);
        }                
}
#endif

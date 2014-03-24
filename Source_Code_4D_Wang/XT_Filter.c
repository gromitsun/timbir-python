#include "allocate.h"
#include <stdio.h>
#include "XT_Structures.h"

void Laplacian_Filter (Real_t h[3][3], Real_t** data, int32_t rows, int32_t cols)
{
	Real_t **temp;
	int32_t i, j, m, n;
	
	temp = (Real_t**)multialloc(sizeof(Real_t), 2, rows, cols);
	for (i = 0; i < rows; i++)
	for (j = 0; j < cols; j++)
	{
		temp[i][j] = 0;
		for (m = -1; m <= 1; m++)
		for (n = -1; n <= 1; n++)
			if (i+m >= 0 && j+n >= 0 && i+m < rows && j+n < cols)
				temp[i][j] += -h[m+1][n+1]*data[i+m][j+n];
	}
	
	for (i = 0; i < rows; i++)
	for (j = 0; j < cols; j++)
		data[i][j] = temp[i][j];		

	multifree(temp,2);
}


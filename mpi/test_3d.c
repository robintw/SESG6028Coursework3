#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "array_alloc.h"

void print_data(int *start, int number, int blocklen, int stride)
{
	int i, j;
	int *curr;
	
// 	for(i = 0; i < number; i++)
// 	{
// 		printf("%d ", *start+(i * stride));
// 		printf("\n");
// 	}

	curr = start;
	
	for(i = 0; i < number; i++)
	{
		for(j = 0; j < blocklen; j++)
		{
			printf("%3d ", *curr);
			curr++;
		}
		printf("\n");
		curr += (stride - 1);
	}
}

int main(int argc, char ** argv)
{
	int ***data;
	int i, j, k;
	int num;
	int a, b, c;
	
	data = alloc_3d_int(2, 3, 4);
	
	num = 1;
	
	for (i = 0; i < 2; i++)
	{
		for (j = 0; j < 3; j++)
		{
			for (k = 0; k < 4; k++)
			{
				data[i][j][k] = num;
				num++;
			}
		}
	}
	
	a = atoi(argv[1]);
	b = atoi(argv[2]);
	c = atoi(argv[3]);
	
	printf("Using count = %d, blocklength = %d and stride = %d\n", a, b, c);
	
	print_data(&data[0][0][0], a, b, c);
	
}

#!/bin/bash

export OMP_NUM_THREADS=4

./laplace <<END
200
200
200
0.01
100
END


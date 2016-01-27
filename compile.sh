#!/bin/bash
gcc -Werror -L/usr/local/Cellar/gsl/1.16/lib -lgsl -lgslcblas int.c basis.c varH.c -o varH.out
./varH.out

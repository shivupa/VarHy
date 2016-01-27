#!/bin/bash
gcc -L/usr/local/Cellar/gsl/1.16/lib -lgsl -lgslcblas varH.c int.c basis.c -o varH.out
./varH.out

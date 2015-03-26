#!/bin/bash

# compile
gcc dg.c -lm -o runDG.out

# run
./runDG.out > result.txt

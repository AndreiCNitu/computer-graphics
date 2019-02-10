#!/bin/bash

make clean
make build
echo -e "\x1B[33m"
printenv | grep OMP_NUM_THREADS
echo -e "\x1B[0m"
./raytracer

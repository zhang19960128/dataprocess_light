#!/bin/bash
g++ -std=c++11 main.cpp atom.cpp -o out
./out 10 dump.xyz cadata.txt
cat result.txt

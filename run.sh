#!/bin/bash
g++ -std=c++11 main.cpp atom.cpp -o ana
./ana 20 dump.xyz cadata.txt
cat result.txt

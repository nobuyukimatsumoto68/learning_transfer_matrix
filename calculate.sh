#!/bin/bash

make
mkdir data
for i in {0..10}
do
    ./a.out ${i} > ./data/L100p${i}.dat
done

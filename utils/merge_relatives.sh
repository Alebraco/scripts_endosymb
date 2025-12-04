#!/bin/bash

for species in cruncher_input/*/; do
    for dir in $species/*/; do
        mv -i $dir/* $species
        rm -rf $dir
    done
done

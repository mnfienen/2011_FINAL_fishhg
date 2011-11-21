#!/bin/sh

currLooInd=$1

tar xzf DATA.tgz
rm DATA.tgz


export LD_LIBRARY_PATH=.
cd data

./production_comps.py $currLooInd

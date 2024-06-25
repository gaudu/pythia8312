#!/bin/bash

cd ~/PYTHIA/pythia8312/examples/
make main1011

e_lab="2 3 4 5 6 7 8 9 10 11 12"
idA="111, 211, -211, 311, 321, -321, 130, 310, 2212, -2212, 2112, -2112, 3122, 3212, 3222, 3112, 3322, 1000020040, 1000070140, 1000260560"
idB="2212, 1000060120, 1000070140, 1000080160, 1000180400"

for energy in $e_lab; do 
    for id1 in $idA; do
        for id2 in $idB; do
            srun \
            --partition=normal \
            --mail-type=END,FAIL,TIME_LIMIT \
            --mail-user="gaudu@uni-wuppertal.de" \
            --job-name="main1011_xsec_1e${energy}_${id1}_${id2}" \
            ~/pleiades_xsec_script_8312.sh $id1 $id2 $energy &
        done
    done
done 


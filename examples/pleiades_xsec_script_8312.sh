#!/bin/bash

cd ~/PYTHIA/pythia8312/examples/

idA=$1
idB=$2
p_lab=$3 # in GeV

#echo "./main1011 $idA $idB $e_lab"
./main1011 $idA $idB $p_lab

#!/bin/bash

cd ~/PYTHIA/pythia8312/examples/
make main426

srun \
--partition=normal \
--mail-type=END,FAIL,TIME_LIMIT \
--mail-user="gaudu@uni-wuppertal.de" \
--job-name="main426" \
~/pleiades_reuse_8312.sh & 

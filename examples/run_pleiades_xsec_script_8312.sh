#!/bin/bash

cd ~/PYTHIA/pythia8312/examples/
make main1011

idA_list=("2212" "-2212" "2112" "-2112" "111" "211" "-211" "311" "321" "-321" "130" "310" "411" "421" "511" "521" "3122" "3222" "3112" "3322" "3334" "1000020040" "1000070140" "1000260560")
idB_list=("1000060120" "1000070140" "1000080160" "1000180400")
p_lab_list=("1e2" "1e3" "1e4" "1e5" "1e6" "1e7" "1e8" "1e9" "1e10" "1e11" "1e12")

for p_lab in "${p_lab_list[@]}"; do 
    for idA in "${idA_list[@]}"; do
        for idB in "${idB_list[@]}"; do
            output_file="main1011_xsec_${p_lab}_${idA}_${idB}.out"
            error_file="main1011_xsec_${p_lab}_${idA}_${idB}.err"
            
            sbatch \
                --partition=normal \
                --mail-type=END,FAIL,TIME_LIMIT \
                --mail-user="gaudu@uni-wuppertal.de" \
                --job-name="main1011_xsec_${p_lab}_${idA}_${idB}" \
                --output=~/slurm/${output_file} \
                --error=~/slurm/${error_file} \
                ~/pleiades_xsec_script_8312.sh $idA $idB $p_lab

            sleep 0.1
        done
    done
done



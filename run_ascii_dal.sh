#!/bin/bash

# run just one small experiment
# g0=/home/calabrese/Projects/mcsplit/McSplit2/dataset/small/mcs10_b03m_s100.A00
# g1=/home/calabrese/Projects/mcsplit/McSplit2/dataset/small/mcs10_b03m_s100.B00
# ./build/mcsplit-dal -t 10 -s closeness min_max $g0 $g1

# exit 0

results_folder="results"
desc=${1:-Default description}

output_folder="$results_folder/results_dal_$(date "+%Y.%m.%d-%H.%M.%S")"
mkdir -p $output_folder
#write description
echo $desc > $output_folder/description.txt
counter=0
for pair in ascii_edgelists/* ; do
    if [ ! -d "$pair" ]; then
        continue
    fi
    counter=$((counter+1))
    if [ $counter -eq 3 -o $counter -eq 5 ]; then
        continue
    fi
    g1="$pair/g1.txt"
    g2="$pair/g2.txt"
    timeout=10
    #timeout=$(cat "$pair/timeout.txt")
    echo "Processing $pair with timeout $timeout"
    pair_name=$(basename $pair)
    outfile="$output_folder/$pair_name.txt"
    ./build/mcsplit-dal -A -t $timeout -s pagerank min_max $g1 $g2 2>&1 | tee $outfile
    echo "timeout: $timeout" >> $outfile
    break
done 
#python3 /home/porro/telecho/telecho.py marco done

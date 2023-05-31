#!/bin/bash

results_folder="results"
desc=${1:-Default description}

output_folder="$results_folder/results_dal_$(date "+%Y.%m.%d-%H.%M.%S")"
mkdir -p $output_folder
#write description
echo $desc > $output_folder/description.txt

# run just one small experiment
outfile="$output_folder/test.txt"
g0=/home/calabrese/Projects/mcsplit/McSplit2/dataset/big/snap_as_s1470.A709.txt
g1=/home/calabrese/Projects/mcsplit/McSplit2/dataset/big/snap_as_s1476.A712.txt
./build/mcsplit-dal -AI 2000000 -s pagerank -p 1 -B 10000000 min_max $g0 $g1 2>&1 | tee $outfile
exit 0

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
    ./build/mcsplit-dal -AI 1000000 -s pagerank -p 1 -B 10000000 min_max $g1 $g2 2>&1 | tee $outfile
    echo "timeout: $timeout" >> $outfile
    break
done 
#python3 /home/porro/telecho/telecho.py done

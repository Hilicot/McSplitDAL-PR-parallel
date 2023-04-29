#!/bin/bash
for pair in /home/licata/McSplit/ascii_edgelists/* ; do
    echo "Processing $pair"
    g1="$pair/g1.txt"
    g2="$pair/g2.txt"
    #timeout=3000
    timeout=$(cat "$pair/timeout.txt")
    /home/licata/McSplit/recur min_max classic -A -t $timeout $g1 $g2 2>&1 | tee "$pair/result.txt"
done 
python3 /home/porro/telecho/telecho.py salvo done

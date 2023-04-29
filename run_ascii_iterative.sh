#!/bin/bash
for pair in /home/licata/McSplit/ascii_edgelists/* ; do
    echo "Processing $pair"
    g1="$pair/g1.txt"
    g2="$pair/g2.txt"
    if [ -f "$pair/iterations.txt" ]; then
        iterations=$(cat "$pair/iterations.txt")
        /home/licata/McSplit/iter min_max classic -A -c -s $iterations $g1 $g2 2>&1 | tee "$pair/result.txt"
    elif [ -f "$pair/timeout.txt" ]; then
        timeout=$(cat "$pair/timeout.txt")
        /home/licata/McSplit/iter min_max classic -A -c -t $timeout $g1 $g2 2>&1 | tee "$pair/result.txt"
    fi
done 
python3 /home/porro/telecho/telecho.py salvo done

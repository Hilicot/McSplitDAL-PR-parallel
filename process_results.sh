for pair in ./ascii_edgelists/* ; do
    res_count=$(ls $pair | grep result | wc -l)
    if [ $res_count -ne 0 ]; then
        for result_file in $pair/result* ; do
            match="Solution size "
            res=$(cat $result_file | grep "$match")
            res=${res//$match/}
            desc=$(tail -n 1 $result_file) 
            echo "Dataset: $pair | Solution size: $res | Description: $desc"
        done
    fi
done
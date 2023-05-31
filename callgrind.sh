valgrind --tool=callgrind --callgrind-out-file=profiler.out ./build/mcsplit-dal -AI 10000 -s pagerank -p 1 -B 1000000 min_max ./ascii_edgelists/test/g1.txt ./ascii_edgelists/test/g2.txt 
callgrind_annotate --tree=both --inclusive=yes --auto=yes --show-percs=yes profiler.out > profiler.txt
rm -rf profiler.out
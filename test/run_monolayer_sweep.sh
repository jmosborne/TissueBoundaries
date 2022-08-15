#!/bin/bash
#
# Script to illustrate running batch jobs and passing in arguments.
#
# 
# This script assumes that the following has been run successfully:
# scons co=1 b=GccOpt ts=projects/TissueBoundaries/test/TestGrowingMonolayer.hpp
#

initial_sim=0;
offset=0;
num_runs=1;
num_sweeps=3;

for (( i=2 ; i<${num_sweeps} ; i++))
do
	start_sim=`expr $i \* $num_runs + $initial_sim`;
	offset_num_runs=`expr $num_runs-$offset`;
	echo $start_sim
    	# NB "nice -20" gives the jobs low priority (good if they are going to dominate the server and no slower if nothing else is going on)
    	# ">" directs std::cout to the file.
    	# "2>&1" directs std::cerr to the same place.
    	# "&" on the end lets the script carry on and not wait until this has finished.
    	nice -20 ../build/optimised/TestGrowingMonolayerRunner -run_index $start_sim -num_runs $offset_num_runs > output/GM_run_${i}_output.txt 2>&1 &
done

echo "Jobs submitted"

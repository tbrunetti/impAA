#!/bin/bash

for i in {1..22};
	do
	head -n1 chr"$i".cut00.results.txt > chr"$i"_results.txt
	for fname in chr"$i".cut*.results.txt
	do
	tail -n+2 $fname >> chr"$i"_results.txt
	done
	head -n1 chr"$i".cut00.mach.info > chr"$i".mach.info
	for fname2 in chr"$i".cut*.mach.info
	do
	tail -n+2 $fname2 >> chr"$i".mach.info
	done
	paste chr"$i"_results.txt chr"$i".mach.info > chr"$i"_results_info.txt
	echo $?
done
cat chr*_results_info.txt | sort -r | tail -n+22 > ./allchr_results_info.txt
echo $?



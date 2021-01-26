#!/bin/bash
#export PERL5LIB=~/TOOLS/vcftools-vcftools-1d27c24/src/perl            ### PAth to the perl library of vcftools
#export PATH=~/TOOLS/htslib-1.5:$PATH			###PATH to the tabix library
module load R
module load torque
seq 0 22 > /gpfs/share/barnescaapa/CAAPA_MEGA/PERU_TGP_MEGA_Analysis/lung_function_asfprefevpp_cases_and_controls_model/chr_list.txt
i=0
while read line; do
((i++))
echo $line
varname="var$i"
printf -v $varname "$line"
done < /gpfs/share/barnescaapa/CAAPA_MEGA/PERU_TGP_MEGA_Analysis/lung_function_asfprefevpp_cases_and_controls_model/chr_list.txt     ### reading the list of Illumina SNP only vcf files with entire path 
for j in `seq 2 $i`; do
curr_var=var$j
echo $curr_var
eval curr_var=\$$curr_var
##echo item: $curr_var
if [ "$curr_var" != "" ]; then
id_name=`echo $curr_var | awk 'END {print $1}'` 		#### Assigning the varaible with each file name 
#echo ${id_name}
name="chr"$id_name
sbatch --mem 8000 --output=${name}.out --error=${name}.err /gpfs/share/barnescaapa/CAAPA_MEGA/PERU_TGP_MEGA_Analysis/lung_function_asfprefevpp_cases_and_controls_model/test_script.sh ${name} ${id_name}                                 ### Submitting the job scripts for each file 
sleep 1 # pause 
fi
done

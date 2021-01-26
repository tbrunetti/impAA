#!/bin/bash
#SBATCH -p defq # Partition
#SBATCH -n 1              # one CPU
#SBATCH -N 1              # on one node   Running time of 1 hours
#SBATCH --share
#SBATCH -t 0-12:00 
#SBATCH --account=barnescaapa-brunettt
##Read the filenames for each set
/gpfs/share/TICR/Software/DosageConvertor/bin/DosageConvertor --vcfDose ${1}.dose${2}.vcf --info ${1}.cut${2}.info --prefix ${1}.cut${2} --type mach --format DS
rm -f ${1}.cut${2}.info ${1}.dose${2}.vcf
gunzip ${1}.cut${2}.mach.dose.gz
cut -f 1,2 -d ":" --output-delimiter=$'\t' ${1}.cut${2}.mach.info > ${1}.cut${2}.test.txt
awk -F" " '{print $1}' ${1}.cut${2}.mach.info > ${1}.cut${2}.test1.txt
awk -F" " '{print $2}' ${1}.cut${2}.test.txt > ${1}.cut${2}.test2.txt
paste ${1}.cut${2}.test1.txt ${1}.cut${2}.test2.txt > ${1}.cut${2}.mach.txt
sed -i -e '1s/REF(0)/position/' ${1}.cut${2}.mach.txt
rm -f ${1}.cut${2}.test.txt ${1}.cut${2}.test2.txt ${1}.cut${2}.test1.txt
name_chr=$1"."
dose=${1}.cut${2}.mach.dose
echo $dose
info=${1}.cut${2}.mach.info
echo $info
posfile=${1}.cut${2}.mach.txt
echo $posfile
chr_num=$(echo $1|sed 's/^chr//')
module load R
module load torque
Rscript /gpfs/share/barnescaapa/CAAPA_MEGA/PERU_TGP_MEGA_Analysis/lung_function_asfprefevpp_cases_and_controls_model/GENESIS_analysis_cases.R --args $name_chr --args $chr_num --args $dose --args $info --args $posfile --args ${2} 

#echo "finally ended, yay"

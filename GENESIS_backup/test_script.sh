#!/bin/bash
#SBATCH -p defq # Partition
#SBATCH -n 1              # one CPU
#SBATCH -N 1              # on one node
#SBATCH -t 0-12:00         # Running time of 12 hours
#SBATCH --share
#SBATCH --account=barnescaapa-brunettt
##Read the filenames for each set
cd /gpfs/share/barnescaapa/CAAPA_MEGA/PERU_TGP_MEGA_Analysis/lung_function_asfprefevpp_cases_and_controls_model #######location where the results of each splitted chromosome is stored
module load gcc-6.1.0-gcc-5.1.0-bglhpmfcjer5a67ideug4bj2ofzzfqyz
  # unzip the file
#gunzip -c /gpfs/barnes_share/dcl01_data_aniket/data/CAAPA_jhuGRAAD_BDOS_032416/WASHINGTON/imputed/${1}.dose.vcf.gz > ${1}.dose.vcf
#gunzip -c /gpfs/barnes_share/dcl01_data_aniket/data/CAAPA_jhuGRAAD_BDOS_032416/WASHINGTON/imputed/${1}.info.gz  > ${1}.info
header1="header1_"${1}
header2="header2_"${1}
variants1="variants1_"${1}
variants2="variants2_"${1}
gunzip -c /gpfs/share/barnescaapa/CAAPA_MEGA/virtualenv-15.0.0/analysis_env/CAAPA_MEGA_analysis_output/imputation_CAAPA_resolved_10042017/imputation_results/PERU_TGP/${1}.dose.vcf.gz|head -n 10000 | grep "^#" > $header1			###grab the headers
gunzip -c /gpfs/share/barnescaapa/CAAPA_MEGA/virtualenv-15.0.0/analysis_env/CAAPA_MEGA_analysis_output/imputation_CAAPA_resolved_10042017/imputation_results/PERU_TGP/${1}.dose.vcf.gz|grep -v "^#" > $variants1				###grab the non-headers	
gunzip -c /gpfs/share/barnescaapa/CAAPA_MEGA/virtualenv-15.0.0/analysis_env/CAAPA_MEGA_analysis_output/imputation_CAAPA_resolved_10042017/imputation_results/PERU_TGP/${1}.info.gz |head -n 10000 | grep "^SNP" > $header2

gunzip -c /gpfs/share/barnescaapa/CAAPA_MEGA/virtualenv-15.0.0/analysis_env/CAAPA_MEGA_analysis_output/imputation_CAAPA_resolved_10042017/imputation_results/PERU_TGP/${1}.info.gz |grep -v "^SNP" > $variants2
  #split into chunks with 100000 lines
split -d -l 100000 $variants1 ${1}.dose
split -d -l 100000 $variants2 ${1}.cut
  #reattach the header to each and clean up
for j in ${1}.dose*;
  do 
    cat $header1 $j >$j.vcf && rm -f $j 
done
  #reattach the header to each and clean up
for j in ${1}.cut*;
  do cat $header2 $j >$j.info && rm -f $j;
done
rm -f $header1 $variants1 $header2 $variants2
numdirs=(${1}.*.info)
numdirs=${#numdirs[@]}
END=$(($numdirs-1))
echo $numdirs
echo $END
for j in $(seq -w 00 $END);
  do
    #/gpfs/barnes_share/Software/DosageConvertor/bin/DosageConvertor --vcfDose ${1}.dose${j}.vcf --info ${1}.cut${j}.info --prefix ${1}.cut${j} --type mach --format DS
    sbatch --mem 60000 --output=${1}.${j}.out --error=${1}.${j}.err dosage_converter.sh ${1} ${j} 
   # rm -f ${1}.cut${j}.info ${1}.dose${j}.vcf
   # gunzip ${1}.cut${j}.mach.dose.gz
   # cut -f 1,2 -d ":" --output-delimiter=$'\t' ${1}.cut${j}.mach.info > ${1}.cut${j}.test.txt
   # awk -F" " '{print $1}' ${1}.cut${j}.mach.info > ${1}.cut${j}.test1.txt
   # awk -F" " '{print $2}' ${1}.cut${j}.test.txt > ${1}.cut${j}.test2.txt
   # paste ${1}.cut${j}.test1.txt ${1}.cut${j}.test2.txt > ${1}.cut${j}.mach.txt
   # sed -i -e '1s/REF(0)/position/' ${1}.cut${j}.mach.txt
    #rm -f ${1}.cut${j}.test.txt ${1}.cut${j}.test2.txt ${1}.cut${j}.test1.txt
  done

#!/bin/bash

#**STEP 2: VARIABLES THAT NEED UPDATING **#
#********* START *********#
dosageConvertor="/full/path/to/executable/DosageConvertor"
doseVcf="/full/path/to/imputed/dose/vcf/myChromosome.dose.vcf.gz"
doseInfo="/full/path/to/imputed/info/file/myChromosome.info.gz"
chrID="myChromsomeNumber"
#********** END **********#


#-----------------------------DO NOT CHANGE CODE PAST THIS LINE-----------------------------#


#**STEP 3: DECOMPRESS FILES AND RUN DOSAGE CONVERTOR**#
#********* START *********#
gunzip -c "${doseVcf}" > ${doseVcf::-3}
gunzip -c "${doseInfo}" > ${doseInfo::-3}
${dosageConvertor} --vcfDose ${doseVcf::-3} --info ${doseInfo::-3} --prefix $chrID --type mach --format "1" DS
gunzip ${chrID}".mach.dose.gz"
#********** END **********#



#**STEP 4: GENERATE SNP POSITION FILE FOR GENESIS**#
#********* START *********#
positionHeader=$(grep -v "^##" ${doseVcf::-3} | head -n 1 | awk -F'\t' '{for (i==1; i<=NF; i++) if($i=="POS") print(i)}')
snpHeader=$(grep -v "^##" ${doseVcf::-3} | head -n 1 | awk -F'\t' '{for (i==1; i<=NF; i++) if($i=="ID") print(i)}')
echo -e "SNP\tposition" > ${chrID}".posFile.txt"
awk -v snp="$snpHeader" -v position="$positionHeader" '{print $snp"\t"$position}' <(grep -v "^#" ${doseVcf::-3}) >> ${chrID}".posFile.txt"
echo "Finished running dosageCovertor!"
#********** END **********#



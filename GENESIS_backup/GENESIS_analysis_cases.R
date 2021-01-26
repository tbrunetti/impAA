#! Rscript --vanilla --default-packages=utils
# load GENESIS file
load("/gpfs/share/barnescaapa/CAAPA_MEGA/PERU_TGP_MEGA_Analysis/PERU_Ayo_Rasika_GENESIS_PC_pheno_covar_data_updated_sex_covars_normalized.Rdata")
library("gdsfmt", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3/")
library("SNPRelate", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3/")
library("GWASTools", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3/")
library("SeqVarTools", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3/")
library("GENESIS", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
args <- commandArgs()
file_name=args[7]
chr_num=args[9]
dosefile=args[11]
markfile=args[13]
posfile=args[15]
subset_val=args[17]
gdsfile <- tempfile()
scanfile <- tempfile()
snpfile <- tempfile()
###association testing
test <- list.files(pattern=glob2rx("chr*.cut*.mach.dose"))
filenames1=test[grep(pattern=file_name,test,fixed=TRUE)]
setwd("/gpfs/share/barnescaapa/CAAPA_MEGA/PERU_TGP_MEGA_Analysis/lung_function_asfprefevpp_cases_and_controls_model")
imputedDosageFile(input.files=c(dosefile,markfile,posfile),filename=gdsfile,chromosome=chr_num,input.type="MaCH",input.dosage=T,file.type="gds")
gds <- GdsGenotypeReader(gdsfile)
genoData <- GenotypeData(gds)
#geno <- getGenotype(genoData)
nullmod.bin <- fitNullMM(scanData = scanAnnot, outcome = "asfprefevpp", covars = c("pc1", "pc2", "pc3", "pc4", "status", "age", "sex", "ses", "log10_asfsfbmi"), covMatList = covMatList, family=gaussian)
myassoc2 <- assocTestMM(genoData = genoData, nullMMobj = nullmod.bin, test="Wald")
data_file=read.table(markfile,header=T,stringsAsFactors=F)
myassoc=cbind(myassoc2,data_file)
write.table(myassoc,file=paste(file_name,"cut",subset_val,".results.txt",sep=''),sep="\t",row.names=F,col.names=T,quote=F)

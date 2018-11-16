library("gdsfmt")
library("SNPRelate")
library("GWASTools")
library("SeqVarTools")
library("GENESIS")


#STEP 2: Variable Assignments
#---------------------<start> VARIABLES THAT NEED TO BE UPDATED <start>------------------#
dataObj <- "/full/path/to/my/R/data/object/containing/scanAnnot/and/coVarMatrix/and/GRM"
phenotype <- "myPhenotype"
varType <- "logistic or gaussian"
chr_num <- "1" # options: string numbers 1-22,"X", "XY", "Y", "M", or "MT"
dosefile <- "/full/path/to/my/converted/chr/dose/vcf/file.txt"
markfile <- "/full/path/to/my/converted/chr/info/file.txt"
posfile <- "/full/path/to/my/positionFile.txt"
covariates <- c("myCov1", "myCov2", "...") # Note these must be available in scanAnnot with exact same header names, comma-separated
#---------------------<end> VARIABLES THAT NEED TO BE UPDATED <end>----------------------#



#STEP 3: File Conversions and Extractions
#-----------------------<start> DATA FORMATTING AND EXTRACTION <start>--------------------#
load(dataObj)
gdsfile <- tempfile()
scanfile <- tempfile()
snpfile <- tempfile()
imputedDosageFile(input.files=c(dosefile,markfile,posfile),filename=gdsfile,chromosome=chr_num,input.type="MaCH",input.dosage=T,file.type="gds")
gds <- GdsGenotypeReader(gdsfile)
genoData <- GenotypeData(gds)
#-------------------------<end> DATA FORMATTING AND EXTRACTION <end>----------------------#



#STEP 4: Association Analysis
#-----------------------<start> GENERATE MODELS, EXTRACT RESULTS, SAVE RESULTS <start>--------------------#
if (tolower(varType) == "logistic"){
  print("Running logistic mixed model...")
  nullmod.bin <- fitNullMM(scanData = scanAnnot, outcome = phenotype, covars = covariates, covMatList = covMatList, family=binomial(link="logit"))
  myassoc2 <- assocTestMM(genoData = genoData, nullMMobj = nullmod.bin, test="Score")
  data_file=read.table(markfile,header=T,stringsAsFactors=F)
  myassoc=cbind(myassoc2,data_file)
  print("Writing results and saving to file...")
  write.table(myassoc,file=paste("chr",chr_num, "_",  phenotype, "_", varType, "_association_results.txt",sep=''),sep="\t",row.names=F,col.names=T,quote=F)
  print("Finshed!")
}else if (tolower(varType) == "gaussian"){
  print("Running Gaussian mixed model...")
  nullmod.bin <- fitNullMM(scanData = scanAnnot, outcome = phenotype, covars = covariates, covMatList = covMatList, family=gaussian)
  myassoc2 <- assocTestMM(genoData = genoData, nullMMobj = nullmod.bin, test="Wald")
  data_file=read.table(markfile,header=T,stringsAsFactors=F)
  myassoc=cbind(myassoc2,data_file)
  print("Writing results and saving to file...")
  write.table(myassoc,file=paste("chr", chr_num, "_",  phenotype, "_", varType, "_association_results.txt",sep=''),sep="\t",row.names=F,col.names=T,quote=F)
  print("Finshed!")
}else{
  stop("ERROR! varType is not logistic or Gaussian, please correct!")
}
#-------------------------<end> GENERATE MODELS, EXTRACT RESULTS, SAVE RESULTS <end>----------------------#


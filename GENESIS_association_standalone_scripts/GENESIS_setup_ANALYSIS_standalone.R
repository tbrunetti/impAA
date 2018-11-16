#source("https://bioconductor.org/biocLite.R")
#biocLite("GWASTools")

library("BiocGenerics")
library("Biobase")
library("GWASTools")
library("SNPRelate")
library("GENESIS")


#STEP 2: Variable Assignments
#---------------------<start> VARIABLES THAT NEED TO BE UPDATED <start>------------------#
#filePrefix <- "/path/to/myPlinkFileName"
phenoFile <- "/path/to/my/phenoFile.txt"
pcmat_num <- as.integer(1) 
fullKin <- TRUE
KINGsoftware <- "/path/to/my/king/executable"
#---------------------<end> VARIABLES THAT NEED TO BE UPDATED <end>----------------------#




#STEP 3: Generating KING output
#-------------------<start> GENERATE KING FILES START <start>------------------#
if (fullKin != TRUE){
  system(paste(KINGsoftware, '-b', paste(filePrefix, '.bed', sep=''), '--prefix', filePrefix), wait = TRUE)
}else{
  system(paste(KINGsoftware, '-b', paste(filePrefix, '.bed', sep=''), '--prefix', filePrefix, '--kinship'), wait = TRUE)
}
#---------------------<end> GENERATE KING FILES START <end>--------------------#




#STEP 4: File format conversions and sample ID generation
#---------------<start> FILE CONVERSIONS AND SAMPLE NAMES <start>---------------#
snpgdsBED2GDS(bed.fn = paste(filePrefix, ".bed", sep=""), bim.fn = paste(filePrefix, ".bim", sep=""), fam.fn = paste(filePrefix, ".fam", sep=""), 
              out.gdsfn = paste(filePrefix, ".gds", sep=""))

file.kin0 <- paste(filePrefix, ".kin0", sep="")
file.kin <- paste(filePrefix, ".kin", sep="")

checkFilePopulation <- nrow(read.table(file.kin)) # make sure there is actual within family structure values

geno <- GdsGenotypeReader(filename = paste(filePrefix, ".gds", sep=""))
genoData <- GenotypeData(geno)
iids <- getScanID(genoData)

fam <- read.table(paste(filePrefix, ".fam", sep=""), sep = "", 
                  col.names = c("FID", "IID", "PAT", "MAT", "SEX", "AFF"))
fam$GENESIS_ID <- seq(1, length(fam$IID))
write.table(fam, file = paste(filePrefix, "_GENESIS_ID_key_file.txt", sep=""), quote = FALSE, col.names = TRUE, row.names = FALSE)
#-----------------<end> FILE CONVERSIONS AND SAMPLE NAMES <end>-----------------#





#STEP 5: Calculate PCs
#--------------<start> PC CALCULATION USING KINSHIP MATRIX <start>--------------#
if((fullKin == TRUE) & (checkFilePopulation > 1)){
  print("Both .kin0 and .kin will be used to generate pcair and pcrelate metrics")
  Kingmat <- king2mat(file.kin0=file.kin0,file.kin=file.kin,type="kinship",iids = iids)
  mypcair <- pcair(genoData = genoData, kinMat = Kingmat,divMat = Kingmat, v=200)
  mypcrel <- pcrelate(genoData = genoData, pcMat = mypcair$vectors[,1:pcmat_num],training.set = mypcair$unrels)
  #pcMat is not the number of PCs you have but instead the number of different admixture populations
}else{
  print("Only .kin0 will be used to generate pcair and pcrelate metrics")
  Kingmat <- king2mat(file.kin0=file.kin0,file.kin=NULL,type="kinship",iids = iids)
  mypcair <- pcair(genoData = genoData, kinMat = Kingmat,divMat = Kingmat, v=200)
  mypcrel <- pcrelate(genoData = genoData, pcMat = mypcair$vectors[,1:pcmat_num],training.set = mypcair$unrels)
  #pcMat is not the number of PCs you have but instead the number of different admixture populations 
}
#----------------<end> PC CALCULATION USING KINSHIP MATRIX <end>----------------#




#STEP 6: Add all phenotypes from phenoFile and other covariates
#-----------<start> CREATE PHENOTYPE COLUMNS FOR N PHENOTYPES/COVARS <start>-----------#
phenoData <- read.table(phenoFile, header = TRUE, na.strings = "NA")
#-------------<end> CREATE PHENOTYPE COLUMNS FOR N PHENOTYPES/COVARS <end>-------------#


#STEP 7: Create the final scanAnnot object to be used in downstream PC plots and GENESIS association analysis
#-----------<start> ADD ALL N PHENOTYPES/COVARS TO DATA OBJECT FROM PHENODATA <start>-----------#
scanAnnot <- ScanAnnotationDataFrame(data.frame(scanID = mypcrel$sample.id,pc1 = mypcair$vectors[,1],pc2 = mypcair$vectors[,2],
	pc3 = mypcair$vectors[,3],pc4 = mypcair$vectors[,4],pc5 = mypcair$vectors[,5],pc6 = mypcair$vectors[,6],
	pc7 = mypcair$vectors[,7],pc8 = mypcair$vectors[,8], pc9 = mypcair$vectors[,9],
	pc10 = mypcair$vectors[,10],pc11 = mypcair$vectors[,11],pc12 = mypcair$vectors[,12],
	pc13 = mypcair$vectors[,13],pc14 = mypcair$vectors[,14],pc15 = mypcair$vectors[,15],
	pc16 = mypcair$vectors[,16],pc17 = mypcair$vectors[,17],pc18 = mypcair$vectors[,18],
	pc19 = mypcair$vectors[,19],pc20 = mypcair$vectors[,20],pc21 = mypcair$vectors[,21],
	pc22 = mypcair$vectors[,22],pc23 = mypcair$vectors[,23],pc24 = mypcair$vectors[,24],
	pc25 = mypcair$vectors[,25],pc26 = mypcair$vectors[,26],pc27 = mypcair$vectors[,27],
	pc28 = mypcair$vectors[,28],pc29 = mypcair$vectors[,29],pc30 = mypcair$vectors[,30],
	pc31 = mypcair$vectors[,31],pc32 = mypcair$vectors[,32],pc33 = mypcair$vectors[,33],
	pc34 = mypcair$vectors[,34],pc35 = mypcair$vectors[,35],pc36 = mypcair$vectors[,36],
	pc37 = mypcair$vectors[,37],pc38 = mypcair$vectors[,38],pc39 = mypcair$vectors[,39],
	pc40 = mypcair$vectors[,40],pc41 = mypcair$vectors[,41],pc42 = mypcair$vectors[,42],
	pc43 = mypcair$vectors[,43],pc44 = mypcair$vectors[,44],pc45 = mypcair$vectors[,45],
	pc46 = mypcair$vectors[,46],pc47 = mypcair$vectors[,47],pc48 = mypcair$vectors[,48],
	pc49 = mypcair$vectors[,49],pc50 = mypcair$vectors[,50],
	age = as.vector(as.matrix(phenoData['age'])),
	sample = as.vector(as.matrix(phenoData['name'])),
	pheno1 = as.vector(as.matrix(phenoData['affected'])),
	pheno2 = as.vector(as.matrix(phenoData['concentration']))
	))
#-------------<end> ADD ALL N PHENOTYPES/COVARS TO DATA OBJECT FROM PHENODATA <end>-------------#



#STEP 8: Create the Genetic Relationships Matrix
#-----------<start> CREATE GENETIC RELATIONSHIPS MATRIX FROM PCRELATE VALUES <start>-----------#
covMatList <- list("Kin" = pcrelateMakeGRM(mypcrel))
#-------------<end> CREATE GENETIC RELATIONSHIPS MATRIX FROM PCRELATE VALUES <end>-------------#



#STEP 9: Save data object and finish
#----------------<start> CREATE AND SAVE BINARY RDATA OBJECT <start>----------------#
# creates a binary file -- can open in R with load()
save(scanAnnot, covMatList, mypcair, file = paste(filePrefix, '_GENESIS_PC_pheno_covar_data.Rdata', sep=""))
#------------------<end> CREATE AND SAVE BINARY RDATA OBJECT <end>------------------#

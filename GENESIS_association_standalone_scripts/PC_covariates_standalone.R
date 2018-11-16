library("BiocGenerics")
library("Biobase")
library("GENESIS")
library("GWASTools")


#STEP 2: Variable Assignments
#---------------------<start> VARIABLES THAT NEED TO BE UPDATED <start>------------------#
dataObj <- "/full/path/to/my/R/data/object/containing/scanAnnot/and/coVarMatrix/and/GRM"
phenotype <- "myPhenotype"
varType <- "logisitc or gaussian"
#---------------------<end> VARIABLES THAT NEED TO BE UPDATED <end>----------------------#




#STEP 3: Confirm Mixed Model with Random Effects Selection Type
#--------------------------<start> PRINT MODEL TYPE TO USE <start>-----------------------#
print(paste("You are running a model based on the phenotype column of scanAnnot called: ", phenotype, sep=""))
print(paste("This is being run on the following model type: ", varType, sep=""))
#----------------------------<end> PRINT MODEL TYPE TO USE <end>-------------------------#




#STEP 4: Check PC Significance
#-----------<start> DETERMINE IF ANY OF FIRST 20 PCs ARE SIGNIFICANTLY ASSCOCIATED WITH PHENOTYPE <start>--------#
load(dataObj)
pc.list <- c()
if(tolower(varType) == "logistic"){
  print(paste("Running logistic model with ", phenotype, " as the outcome/response/dependent variable and each of the first 20 PCs as the predictor/explanatory/independent variable", sep=""))
  for (i in 1:20) {
    pc.str <- paste0("pc", i)
    nullmod <- fitNullMM(scanData = scanAnnot, 
                         outcome = phenotype, 
                         covars = c(pc.str), 
                         covMatList = covMatList, 
                         family=binomial(link = "logit"))
    p.val <- nullmod$fixef[2, "pval"]
    if (p.val < 0.05) {
      pc.list <- c(pc.list, pc.str)
    }
  }
}else if(tolower(varType) == "gaussian"){
  print(paste("Running Gaussian model with ", phenotype, " as the outcome/response/dependent variable and each of the first 20 PCs as the predictor/explanatory/independent variable", sep=""))
  for (i in 1:20) {
    pc.str <- paste0("pc", i)
    nullmod <- fitNullMM(scanData = scanAnnot, 
                         outcome = phenotype, 
                         covars = c(pc.str), 
                         covMatList = covMatList, 
                         family=gaussian)
    p.val <- nullmod$fixef[2, "pval"]
    if (p.val < 0.05) {
      pc.list <- c(pc.list, pc.str)
    }
  }
}else{
  print(paste("ERROR! The following model", varType, "does not exist!", sep = " "))
  stop("Please confirm varType in STEP 2 is set to logistic or gaussian!!")
}
#-------------<end> DETERMINE IF ANY OF FIRST 20 PCs ARE SIGNIFICANTLY ASSCOCIATED WITH PHENOTYPE <end>----------#




#STEP 5: Report PC Significance
#-----------<start> LIST FIRST 20 PCs ARE SIGNIFICANTLY ASSCOCIATED WITH PHENOTYPE <start>--------#
print(paste("There are length ", length(pc.list), " PCs that are singificantly asscoiated (p<0.05) with your phenotype of interest"))
if (length(pc.list) != 0){
  print("The following PCs should be included as covariates in your association analysis model:")
  for (i in pc.list){
    print(i)
    }
  }else if (length(pc.list) == 0){
    print("No additional PCs need to be added as covariates in your model.")
  }
#-------------<end> LIST FIRST 20 PCs ARE SIGNIFICANTLY ASSCOCIATED WITH PHENOTYPE <end>----------#

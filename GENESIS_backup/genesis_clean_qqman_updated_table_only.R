#install.packages("data.table", repos="cran.r-project.org")
library("qqman", lib.loc='/homelink/brunettt/R/x86_64-pc-linux-gnu-library/3.3/')
#library("data.table", lib.loc = '/homelink/brunettt/R/x86_64-pc-linux-gnu-library/3.3/')

# name of input file
dat.info <- read.table('allchr_results_info.txt', header=T)
colnames(dat.info)
dim(dat.info)

data_output_prefix = 'Peru_log10_case_control_asthma_model'
#number of samples in analysis
total_samples = 727

#mark flag 1
dat.info$flag1<-ifelse((as.numeric(as.character(dat.info$Rsq)) <= 0.5 & dat.info$MAF <= 0.005),"flag1","clean")
summary(as.factor(dat.info$flag1))

#mark flag 2
dat.info$flag2<-ifelse((as.numeric(as.character(dat.info$Rsq)) <= 0.3 & dat.info$MAF > 0.005),"flag2","clean")
summary(as.factor(dat.info$flag2))

# create dataset where flag1 and flag2 are "clean" -- use this to generate plots
dat.info.out2<-dat.info[which(dat.info$flag1 == "clean" & dat.info$flag2 == "clean"),]
dim(dat.info.out2)

# format SNP and POS columns
dat.info.out <- dat.info.out2
dat.info.out$SNP<-as.character(dat.info.out$SNP)
dat.info.out$POS <- as.numeric(unlist(strsplit(dat.info.out$SNP, split=":"))[seq(2,dim(dat.info.out)[1]*2,2)])
dat.info.out<-dat.info.out[which(dat.info.out$Wald.pval > 1e-40),]
dat.info.out$MAC <- 2*as.numeric(total_samples)*dat.info.out$MAF
#write.table(dat.info.out, paste(data_output_prefix, "allchr_output_genesis_info_clean.txt", sep="_"),sep="\t",row.names=F,quote=F)

dat.info.out.mac<-dat.info.out[which(dat.info.out$MAC>6),]
write.table(dat.info.out.mac,paste(data_output_prefix, "allchr_output_genesis_info_clean_mac_filtered.txt", sep="_"),sep="\t",row.names=F,quote=F)
dim(dat.info.out.mac)

#
##common SNPs by MAF
##--------------------
#dat.info.out.common.maf<-dat.info.out[which(dat.info.out$MAF >= 0.05),]
#dim(dat.info.out.common.maf)
#summary(dat.info.out.common.maf$Wald.pval)
#
#
##common SNPs by MAC
##--------------------
#dat.info.out2<-dat.info.out[which(dat.info.out$MAC>6),]
#dat.info.out.common.mac<-dat.info.out2[which(dat.info.out2$MAF >= 0.05),]
#dim(dat.info.out.common.mac)
#summary(dat.info.out.common.mac$Wald.pval)
#
#
##rare SNPs by MAF
##-----------------
#dat.info.out.rare<-dat.info.out[which(dat.info.out$MAF < 0.05),]
#dim(dat.info.out.rare)
#summary(dat.info.out.rare$Wald.pval)
#
##rare SNPs by MAC
##-----------------
#dat.info.out2<-dat.info.out[which(dat.info.out$MAC>6),]
#dat.info.out.rare.mac<-dat.info.out2[which(dat.info.out2$MAF < 0.05),]
#dim(dat.info.out.rare.mac)
#summary(dat.info.out.common.mac$Wald.pval)
#
#
##all common and rare variants qq plot
#observed <- sort(dat.info.out$Wald.pval)
#lobs <- -(log10(observed))
#expected <- c(1:length(observed)) 
#lexp <- -(log10(expected / (length(expected)+1)))
#bmp(filename=paste(data_output_prefix, "qqplot_common_and_rare.bmp", sep="_"), width=800, height=800, bg="white", type="cairo")
#plot(c(0,10), c(0,10), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,10), ylim=c(0,10), las=1, xaxs="i", yaxs="i", bty="l")
#points(lexp, lobs, pch=23, cex=.4, bg="black") 
##inflation factor lambda
#chisq2 <- qchisq(1-dat.info.out$Wald.pval,1,lower.tail = T)
#lambda <- median(chisq2,na.rm=T)/qchisq(0.5,1)#lambda1
#mtext(bquote(lambda == .(lambda)), side=3, cex=2.0)
#dev.off()
#
##all common and rare variants manhattan plot
#library(qqman)
#bmp(filename=paste(data_output_prefix, "manhattan_common_and_rare.bmp", sep="_"), width=800, height=600, bg="white", type="cairo")
#par(font.axis = 2)
#manhattan(dat.info.out,chr = "chr", bp = "POS", p = "Wald.pval", snp = "SNP_hwe_all",col = c("gray60", "gray10"), chrlabs = NULL,highlight = NULL, logp = TRUE,suggestiveline = F, genomewideline = F,ylim=c(0,10),paste(main="Association analysis:Genesis", data_output_prefix, sep=" "))
#dev.off()
#
#
#
#
##common maf qq plot
#observed <- sort(dat.info.out.common.maf$Wald.pval)
#lobs <- -(log10(observed))
#expected <- c(1:length(observed)) 
#lexp <- -(log10(expected / (length(expected)+1)))
#bmp(filename=paste(data_output_prefix, "qqplot_common_variants.bmp", sep="_"), width=800, height=800, bg="white", type="cairo")
#plot(c(0,10), c(0,10), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,10), ylim=c(0,10), las=1, xaxs="i", yaxs="i", bty="l")
#points(lexp, lobs, pch=23, cex=.4, bg="black") 
##inflation factor lambda
#chisq2 <- qchisq(1-dat.info.out.common.maf$Wald.pval,1,lower.tail = T)
#lambda <- median(chisq2,na.rm=T)/qchisq(0.5,1)#lambda1
#mtext(bquote(lambda == .(lambda)), side=3, cex=2.0)
#dev.off()
#
##common maf manhattan plot
##library(qqman)
#bmp(filename=paste(data_output_prefix, "manhattan_common_variants.bmp", sep="_"), width=800, height=600, bg="white", type="cairo")
#par(font.axis = 2)
#manhattan(dat.info.out.common.maf,chr = "chr", bp = "POS", p = "Wald.pval", snp = "SNP_hwe_all",col = c("gray60", "gray10"), chrlabs = NULL,highlight = NULL, logp = TRUE,suggestiveline = F, genomewideline = F, ylim=c(0,10) ,main=paste("Association analysis:Genesis common variants", data_output_prefix, sep=" "))
#dev.off()
#
#
##common maf qq plot with MAC filter
#observed <- sort(dat.info.out.common.mac$Wald.pval)
#lobs <- -(log10(observed))
#expected <- c(1:length(observed)) 
#lexp <- -(log10(expected / (length(expected)+1)))
#bmp(filename=paste(data_output_prefix, "qqplot_common_variants_MAC_filtered.bmp", sep="_"), width=800, height=800, bg="white", type="cairo")
#plot(c(0,10), c(0,10), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,10), ylim=c(0,10), las=1, xaxs="i", yaxs="i", bty="l")
#points(lexp, lobs, pch=23, cex=.4, bg="black") 
##inflation factor lambda
#chisq2 <- qchisq(1-dat.info.out.common.mac$Wald.pval,1,lower.tail = T)
#lambda <- median(chisq2,na.rm=T)/qchisq(0.5,1)#lambda1
#mtext(bquote(lambda == .(lambda)), side=3, cex=2.0)
#dev.off()
#
#
##common maf manhattan plot with MAC filter
##library(qqman)
#bmp(filename=paste(data_output_prefix, "manhattan_common_variants_MAC_filtered.bmp", sep="_"), width=800, height=600, bg="white", type="cairo")
#par(font.axis = 2)
#manhattan(dat.info.out.common.mac,chr = "chr", bp = "POS", p = "Wald.pval", snp = "SNP_hwe_all",col = c("gray60", "gray10"), chrlabs = NULL,highlight = NULL, logp = TRUE,suggestiveline = F, genomewideline = F, ylim=c(0,10) ,main=paste("Association analysis:Genesis common MAC filtered variants", data_output_prefix, sep=" "))
#dev.off()
#
#
#
##rare maf qq plot
#observed <- sort(dat.info.out.rare$Wald.pval)
#lobs <- -(log10(observed))
#expected <- c(1:length(observed)) 
#lexp <- -(log10(expected / (length(expected)+1)))
#bmp(filename=paste(data_output_prefix, "qqplot_rare_variants.bmp", sep="_"), width=800, height=800, bg="white", type="cairo")
#plot(c(0,10), c(0,10), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,10), ylim=c(0,10), las=1, xaxs="i", yaxs="i", bty="l")
#points(lexp, lobs, pch=23, cex=.4, bg="black") 
##inflation factor lambda
#chisq2 <- qchisq(1-dat.info.out.rare$Wald.pval,1,lower.tail = T)
#lambda <- median(chisq2,na.rm=T)/qchisq(0.5,1)#lambda1
#mtext(bquote(lambda == .(lambda)), side=3, cex=2.0)
#dev.off()
#
#
##rare maf manhattan plot
##library(qqman)
#bmp(filename=paste(data_output_prefix, "manhattan_rare_variants.bmp", sep="_"), width=800, height=600, bg="white", type="cairo")
#par(font.axis = 2)
#manhattan(dat.info.out.rare,chr = "chr", bp = "POS", p = "Wald.pval", snp = "SNP_hwe_all",col = c("gray60", "gray10"), chrlabs = NULL,highlight = NULL, logp = TRUE,suggestiveline = F, genomewideline = F, ylim=c(0,10) ,main=paste("Association analysis:Genesis rare variants", data_output_prefix, sep=" "))
#dev.off()
#
#
##rare maf qq plot with MAC filter
#observed <- sort(dat.info.out.rare.mac$Wald.pval)
#lobs <- -(log10(observed))
#expected <- c(1:length(observed)) 
#lexp <- -(log10(expected / (length(expected)+1)))
#bmp(filename=paste(data_output_prefix, "qqplot_rare_variants_MAC_filtered.bmp", sep="_"), width=800, height=800, bg="white", type="cairo")
#plot(c(0,10), c(0,10), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,10), ylim=c(0,10), las=1, xaxs="i", yaxs="i", bty="l")
#points(lexp, lobs, pch=23, cex=.4, bg="black") 
##inflation factor lambda
#chisq2 <- qchisq(1-dat.info.out.rare.mac$Wald.pval,1,lower.tail = T)
#lambda <- median(chisq2,na.rm=T)/qchisq(0.5,1)#lambda1
#mtext(bquote(lambda == .(lambda)), side=3, cex=2.0)
#dev.off()
#
#
##rare maf manhattan plot with MAC filter
##library(qqman)
#bmp(filename=paste(data_output_prefix, "manhattan_rare_variants_MAC_filtered.bmp", sep="_"), width=800, height=600, bg="white", type="cairo")
#par(font.axis = 2)
#manhattan(dat.info.out.rare.mac,chr = "chr", bp = "POS", p = "Wald.pval", snp = "SNP_hwe_all",col = c("gray60", "gray10"), chrlabs = NULL,highlight = NULL, logp = TRUE,suggestiveline = F, genomewideline = F, ylim=c(0,10) ,main=paste("Association analysis:Genesis rare MAC filtered variants", data_output_prefix, sep=" "))
#dev.off()
#
##all common and rare variants qq plot with MAC filter
#observed <- sort(dat.info.out.mac$Wald.pval)
#lobs <- -(log10(observed))
#expected <- c(1:length(observed)) 
#lexp <- -(log10(expected / (length(expected)+1)))
#bmp(filename=paste(data_output_prefix, "qqplot_common_and_rare_MAC_filter.bmp", sep="_"), width=800, height=800, bg="white", type="cairo")
#plot(c(0,10), c(0,10), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,10), ylim=c(0,10), las=1, xaxs="i", yaxs="i", bty="l")
#points(lexp, lobs, pch=23, cex=.4, bg="black") 
##inflation factor lambda
#chisq2 <- qchisq(1-dat.info.out.mac$Wald.pval,1,lower.tail = T)
#lambda <- median(chisq2,na.rm=T)/qchisq(0.5,1)#lambda1
#mtext(bquote(lambda == .(lambda)), side=3, cex=2.0)
#dev.off()
#
##all common and rare variants manhattan plot with MAC fileter
##library(qqman)
#bmp(filename=paste(data_output_prefix, "manhattan_common_and_rare_MAC_filter.bmp", sep="_"), width=800, height=600, bg="white", type="cairo")
#par(font.axis = 2)
#manhattan(dat.info.out.mac,chr = "chr", bp = "POS", p = "Wald.pval", snp = "SNP_hwe_all",col = c("gray60", "gray10"), chrlabs = NULL,highlight = NULL, logp = TRUE,suggestiveline = F, genomewideline = F,ylim=c(0,10),main=paste("Association analysis:Genesis common and rare MAC filtered variants", data_output_prefix, sep=" "))
#dev.off()
#

library("qqman", lib.loc='/homelink/brunettt/R/x86_64-pc-linux-gnu-library/3.3/')
dat.info<-read.delim("allchr_results_info.txt")
colnames(dat.info)
dim(dat.info)

dat.info$flag1<-ifelse((as.numeric(as.character(dat.info$Rsq)) <= 0.5 & dat.info$MAF <= 0.005),"flag1","clean")
summary(as.factor(dat.info$flag1))

dat.info$flag2<-ifelse((as.numeric(as.character(dat.info$Rsq)) <= 0.3 & dat.info$MAF > 0.005),"flag2","clean")
summary(as.factor(dat.info$flag2))

dat.info.out2<-dat.info[which(dat.info$flag1 == "clean" & dat.info$flag2 == "clean"),]
dim(dat.info.out2)
write.table(dat.info.out2,"allchr_output_genesis_info_clean_PERU_asthma_model.txt",sep="\t",row.names=F,quote=F)
write.table(dat.info,"allchr_output_genesis_info_PERU_asthma_model.txt",sep="\t",row.names=F,quote=F)

#
#
dat.info.out <- dat.info.out2
dat.info.out$SNP<-as.character(dat.info.out$SNP)
dat.info.out$POS <- as.numeric(unlist(strsplit(dat.info.out$SNP, split=":"))[seq(2,dim(dat.info.out)[1]*2,2)])

dat.info.out<-dat.info.out[which(dat.info.out$Score.pval > 1e-40),]

#fmiss and hwe clean
#dat.info.out<-dat.info.out[which(dat.info.out$F_MISS < 0.05),]
#dat.info.out<-dat.info.out[which(dat.info.out$P > 0.000001),]
#all cleaned SNPs
dim(dat.info.out)
summary(dat.info.out$Score.pval)


#common SNPs
dat.info.out.common<-dat.info.out[which(dat.info.out$MAF.1 >= 0.05),]
dim(dat.info.out.common)
summary(dat.info.out.common$Score.pval)

#rare SNPs
dat.info.out.rare<-dat.info.out[which(dat.info.out$MAF.1 < 0.05),]
dim(dat.info.out.rare)
summary(dat.info.out.rare$Score.pval)

#all qq plot
observed <- sort(dat.info.out$Score.pval)
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
bmp(filename="qqplot_typed_overlap_allcohort_new_PERU_asthma_model.bmp", width=800, height=800, bg="white", type="cairo")
#pdf("qqplot_typed_overlap_allcohort.pdf", width=6, height=6)
plot(c(0,10), c(0,10), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,10), ylim=c(0,10), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black") 
#inflation factor lambda
chisq2 <- qchisq(1-dat.info.out$Score.pval,1,lower.tail = T)
lambda <- median(chisq2,na.rm=T)/qchisq(0.5,1)#lambda1
mtext(bquote(lambda == .(lambda)), side=3, cex=2.0)
dev.off()

#manhattan plot
library(qqman)
bmp(filename="manhattan_typed_overlap_allcohort_new_PERU_asthma_model.bmp", width=800, height=600, bg="white", type="cairo")
#pdf("manhattan_typed_overlap_allcohort.pdf",width=21,height=10)
par(font.axis = 2)
manhattan(dat.info.out,chr = "chr", bp = "POS", p = "Score.pval", snp = "SNP_hwe_all",col = c("gray60", "gray10"), chrlabs = NULL,highlight = NULL, logp = TRUE,suggestiveline = F, genomewideline = F,ylim=c(0,10),main="Association analysis:Genesis PERU asthma model")
dev.off()

#common qq plot
observed <- sort(dat.info.out.common$Score.pval)
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
bmp(filename="qqplot_typed_overlap_allcohort_common_new_PERU_asthma_model.bmp", width=800, height=800, bg="white", type="cairo")
#pdf("qqplot_typed_overlap_allcohort_common.pdf", width=6, height=6)
plot(c(0,10), c(0,10), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,10), ylim=c(0,10), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black") 
#inflation factor lambda
chisq2 <- qchisq(1-dat.info.out.common$Score.pval,1,lower.tail = T)
lambda <- median(chisq2,na.rm=T)/qchisq(0.5,1)#lambda1
mtext(bquote(lambda == .(lambda)), side=3, cex=2.0)
dev.off()

#manhattan plot common
library(qqman)
bmp(filename="manhattan_typed_overlap_allcohort_common_new_PERU_asthma_model.bmp", width=800, height=600, bg="white", type="cairo")
#pdf("manhattan_typed_overlap_allcohort_common_new.pdf",width=21,height=10)
par(font.axis = 2)
manhattan(dat.info.out.common,chr = "chr", bp = "POS", p = "Score.pval", snp = "SNP_hwe_all",col = c("gray60", "gray10"), chrlabs = NULL,highlight = NULL, logp = TRUE,suggestiveline = F, genomewideline = F, ylim=c(0,10) ,main="Association analysis:Genesis PERU asthma model common variants")
dev.off()

#rare
observed <- sort(dat.info.out.rare$Score.pval)
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
bmp(filename="qqplot_typed_overlap_allcohort_rare_new_PERU_asthma_model.bmp", width=800, height=800, bg="white", type="cairo")
#pdf("qqplot_typed_overlap_allcohort_rare.pdf", width=6, height=6)
plot(c(0,10), c(0,10), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,10), ylim=c(0,10), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black") 
#inflation factor lambda
chisq2 <- qchisq(1-dat.info.out.rare$Score.pval,1,lower.tail = T)
lambda <- median(chisq2,na.rm=T)/qchisq(0.5,1)#lambda1
dev.off()

#manhattan plot rare
library(qqman)
bmp(filename="manhattan_typed_overlap_allcohort_rare_new_PERU_asthma_model.bmp", width=800, height=600, bg="white", type="cairo")
#pdf("manhattan_typed_overlap_allcohort_rare.pdf",width=21,height=10)
par(font.axis = 2)
manhattan(dat.info.out.rare,chr = "chr", bp = "POS", p = "Score.pval", snp = "SNP_hwe_all",col = c("gray60", "gray10"), chrlabs = NULL,highlight = NULL, logp = TRUE,suggestiveline = F, genomewideline = F, ylim=c(0,10) ,main="Association analysis:Genesis PERU asthma model rare variants")
dev.off()


#qq plot all snps with CI

bmp(filename="qqplot_typed_overlap_allcohort_new_with_CI_PERU_asthma_model.bmp", width=800, height=800, bg="white", type="cairo")

## obs <- readfile; p-values only
## read in your p-values,
## here I generated some
obs<- dat.info.out$Score.pval
N <- 1000000 ## number of p-values
## create the null distribution
## (-log10 of the uniform)
null <- -log(1:N/N,10)
MAX <- max(c(obs,null))
## create the confidence intervals
c95 <- rep(0,N)
c05 <- rep(0,N)
## the jth order statistic from a
## uniform(0,1) sample
## has a beta(j,n-j+1) distribution
## (Casella & Berger, 2002,
## 2nd edition, pg 230, Duxbury)
for(i in 1:N){
  c95[i] <- qbeta(0.95,i,N-i+1)
  c05[i] <- qbeta(0.05,i,N-i+1)
}
## plot the two confidence lines
plot(null, -log(c95,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l",
     axes=FALSE, xlab="", ylab="")
par(new=T)
plot(null, -log(c05,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l",
     axes=FALSE, xlab="", ylab="")
## add the diagonal
abline(0,1,col="red")
par(new=T)

observed <- sort(dat.info.out$Score.pval)
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
#pdf("qqplot_typed_overlap_allcohort.pdf", width=6, height=6)
plot(c(0,10), c(0,10), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,10), ylim=c(0,10), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black") 
chisq2 <- qchisq(1-dat.info.out$Score.pval,1,lower.tail = T)
lamba <- median(chisq2,na.rm=T)/qchisq(0.5,1)#lambda1
mtext(bquote(lambda == .(lambda)), side=3, cex=2.0)
dev.off()

# qq common with CI
bmp(filename="qqplot_typed_overlap_allcohort_common_new_with_CI_PERU_asthma_model.bmp", width=800, height=800, bg="white", type="cairo")

obs<- dat.info.out.common$Score.pval
N <- 1000000 ## number of p-values
## create the null distribution
## (-log10 of the uniform)
null <- -log(1:N/N,10)
MAX <- max(c(obs,null))
## create the confidence intervals
c95 <- rep(0,N)
c05 <- rep(0,N)
## the jth order statistic from a
## uniform(0,1) sample
## has a beta(j,n-j+1) distribution
## (Casella & Berger, 2002,
## 2nd edition, pg 230, Duxbury)
for(i in 1:N){
  c95[i] <- qbeta(0.95,i,N-i+1)
  c05[i] <- qbeta(0.05,i,N-i+1)
}
## plot the two confidence lines
plot(null, -log(c95,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l",
     axes=FALSE, xlab="", ylab="")
par(new=T)
plot(null, -log(c05,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l",
     axes=FALSE, xlab="", ylab="")
## add the diagonal
abline(0,1,col="red")
par(new=T)


observed <- sort(dat.info.out.common$Score.pval)
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
#pdf("qqplot_typed_overlap_allcohort_common.pdf", width=6, height=6)
plot(c(0,10), c(0,10), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,10), ylim=c(0,10), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black") 
chisq2 <- qchisq(1-dat.info.out.common$Score.pval,1,lower.tail = T)
lambda <- median(chisq2,na.rm=T)/qchisq(0.5,1)#lambda1
mtext(bquote(lambda == .(lambda)), side=3, cex=2.0)
dev.off()
# qq plot rare with CI
bmp(filename="qqplot_typed_overlap_allcohort_rare_new_with_CI_PERU_asthma_model.bmp", width=800, height=800, bg="white", type="cairo")

obs<- dat.info.out.rare$Score.pval
N <- 1000000 ## number of p-values
## create the null distribution
## (-log10 of the uniform)
null <- -log(1:N/N,10)
MAX <- max(c(obs,null))
## create the confidence intervals
c95 <- rep(0,N)
c05 <- rep(0,N)
## the jth order statistic from a
## uniform(0,1) sample
## has a beta(j,n-j+1) distribution
## (Casella & Berger, 2002,
## 2nd edition, pg 230, Duxbury)
for(i in 1:N){
  c95[i] <- qbeta(0.95,i,N-i+1)
  c05[i] <- qbeta(0.05,i,N-i+1)
}
## plot the two confidence lines
plot(null, -log(c95,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l",
     axes=FALSE, xlab="", ylab="")
par(new=T)
plot(null, -log(c05,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l",
     axes=FALSE, xlab="", ylab="")
## add the diagonal
abline(0,1,col="red")
par(new=T)

observed <- sort(dat.info.out.rare$Score.pval)
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
#pdf("qqplot_typed_overlap_allcohort_rare.pdf", width=6, height=6)
plot(c(0,10), c(0,10), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,10), ylim=c(0,10), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black") 
chisq2 <- qchisq(1-dat.info.out.rare$Score.pval,1,lower.tail = T)
lambda <- median(chisq2,na.rm=T)/qchisq(0.5,1)#lambda1
mtext(bquote(lambda == .(lambda)), side=3, cex=2.0)
dev.off()

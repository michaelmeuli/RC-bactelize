#!/usr/bin/Rscript

sink("/home/mmeuli/batch/out/1-1-A_result-R.txt")
ImgData <- read.table("/home/mmeuli/batch/out/1-1-A_result.txt", header=TRUE)
ImgDataSub1 <- subset(ImgData, physicalSize > 1 & physicalSize < 9)
write.table(ImgDataSub1, file="/home/mmeuli/batch/out/1-1-A_result-subset1.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#https://cran.r-project.org/web/packages/VGAM/VGAM.pdf
#Estimates the five parameters of a mixture of two univariate normal distributions by maximum likelihood estimation.
#By default, the five linear/additive predictors are (logit(φ), µ1, log(σ1), µ2, log(σ2))T.
#Threshold to decide if bacteria is inside Macrophage:

macrophageLog <- log(ImgDataSub1$mean_macrophage) 
m2.estimation <- median(macrophageLog)  
m1.estimation <- m2.estimation/2
cat("Excluding bacteria outside of macrophages...\n")
cat("Estimated value vor m2:\t", m2.estimation, "\n")
cat("Estimated value for m1:\t", m1.estimation, "\n")
#install.packages("VGAM")
require(VGAM, lib.loc = "/home/mmeuli/R/x86_64-pc-linux-gnu-library/3.2")
init.params <- list(m1=m1.estimation, m2=m2.estimation)
macrophageData <- data.frame(macrophageLog)
sink()
fit <- vglm(macrophageLog ~ 1, mix2normal(eq.sd=FALSE, imu1=init.params$m1, imu2=init.params$m2), macrophageData)
sink("/home/mmeuli/batch/out/1-1-A_result-R.txt", append = TRUE)
pars <- as.vector(coef(fit))
m1 <- pars[2]
sd1 <- exp(pars[3])
m2 <- pars[4]
sd2 <- exp(pars[5])
threshold <- (m2-(2*sd2))
cat("Calculated value vor m2:\t", m2, "\n")
cat("Calculated value vor m1:\t", m1, "\n")
cat("Threshold to exclude bacteria outside of macrophages:\t", threshold, "\n\n")

sink()
png(filename="/home/mmeuli/batch/out/1-1-A_result-hist-macrophages.png")
hist(macrophageLog, breaks=100)
abline(v=threshold, col=2)
abline(v=m2, col=2)
dev.off()
sink("/home/mmeuli/batch/out/1-1-A_result-R.txt", append = TRUE)

ImgDataSub2 <- subset(ImgDataSub1, log(mean_macrophage) > threshold)
write.table(ImgDataSub2, file="/home/mmeuli/batch/out/1-1-A_result-subset2.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
ImgDataSubBac <- subset(ImgDataSub2, (maxDiameter/minDiameter) > 1.7 & roundness > 0.3)
write.table(ImgDataSubBac, file="/home/mmeuli/batch/out/1-1-A_result-subsetBac.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
#ImgDataSubDeb <- ImgDataSub2[!ImgDataSubBac]
ImgDataSubDeb <- subset(ImgDataSub2, roundness <= 0.3)
write.table(ImgDataSubDeb, file="/home/mmeuli/batch/out/1-1-A_result-subsetDeb.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

C1 <- subset(ImgDataSub2, coverslipNr >= 1 & coverslipNr <= 3)
C2 <- subset(ImgDataSub2, coverslipNr >= 4 & coverslipNr <= 6)
C3 <- subset(ImgDataSub2, coverslipNr >= 7 & coverslipNr <= 9)
cat("Median of mean_lysosome of wt condition Sub2", median(C1$mean_lysosome), "\n")
cat("Median of mean_lysosome of ko condition Sub2", median(C3$mean_lysosome), "\n")

C1log <- log(C1$mean_lysosome)   
C3log <- log(C3$mean_lysosome)
wilcox.test(C1log, C3log, paired=F)


C1Bac <- subset(ImgDataSubBac, coverslipNr >= 1 & coverslipNr <= 3)
C2Bac <- subset(ImgDataSubBac, coverslipNr >= 4 & coverslipNr <= 6)
C3Bac <- subset(ImgDataSubBac, coverslipNr >= 7 & coverslipNr <= 9)
cat("Median of mean_lysosome of wt condition SubBac", median(C1Bac$mean_lysosome), "\n")
cat("Median of mean_lysosome of ko condition SubBac", median(C3Bac$mean_lysosome), "\n")

C1Baclog <- log(C1Bac$mean_lysosome)   
C3Baclog <- log(C3Bac$mean_lysosome)
wilcox.test(C1Baclog, C3Baclog, paired=F)


C1Deb <- subset(ImgDataSubDeb, coverslipNr >= 1 & coverslipNr <= 3)
C2Deb <- subset(ImgDataSubDeb, coverslipNr >= 4 & coverslipNr <= 6)
C3Deb <- subset(ImgDataSubDeb, coverslipNr >= 7 & coverslipNr <= 9)
cat("Median of mean_lysosome of wt condition SubDeb", median(C1Deb$mean_lysosome), "\n")
cat("Median of mean_lysosome of ko condition SubDeb", median(C3Deb$mean_lysosome), "\n")

C1Deblog <- log(C1Deb$mean_lysosome)   
C3Deblog <- log(C3Deb$mean_lysosome)
wilcox.test(C1Deblog, C3Deblog, paired=F)


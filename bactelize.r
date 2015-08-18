#!/usr/bin/Rscript

sink("/home/mmeuli/batch/out/1-1-A_result-R.txt")

ImgData <- read.table("/home/mmeuli/batch/out/1-1-A_result.txt", header=TRUE)
ImgDataSub1 <- subset(ImgData, physicalSize > 1 & physicalSize < 9 & (maxDiameter/minDiameter) > 1.7)
write.table(ImgDataSub1, file="/home/mmeuli/batch/out/1-1-A_result-subset1.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

median.macro <- median(ImgDataSub1$mean_macrophage)
ImgDataSub2 <- subset(ImgDataSub1, mean_macrophage > (median.macro/2))
ImgDataSubEx <- subset(ImgDataSub1, mean_macrophage <= (median.macro/2))
write.table(ImgDataSub2, file="/home/mmeuli/batch/out/1-1-A_result-subset2.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(ImgDataSubEx, file="/home/mmeuli/batch/out/1-1-A_result-subsetEx.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

C1 <- subset(ImgDataSub2, coverslipNr >= 1 & coverslipNr <= 3)
C2 <- subset(ImgDataSub2, coverslipNr >= 4 & coverslipNr <= 6)
C3 <- subset(ImgDataSub2, coverslipNr >= 7 & coverslipNr <= 9)
cat("Median mean_lysosome of C1\n")
median(C1$mean_lysosome)
cat("Median mean_lysosome of C3\n")
median(C3$mean_lysosome)

C1log <- log(C1$mean_lysosome)   
C3log <- log(C3$mean_lysosome)
wilcox.test(C1log, C3log, paired=F)


#LANG=C sh -c 'R_DEFAULT_PACKAGES="$R_DEFAULT_PACKAGES Rcmdr" R "$@"'

#sudo apt-get install r-cran-rjava (required by XLConnect, which is needed for Rcmdr)
#install.packages("Rcmdr")
#install.packages("relimp")
#install.packages("VGAM")

ImgData0 <- read.table("/home/mmeuli/batch/out/1-1-A_result.txt", header=TRUE, sep="", na.strings="-nan", dec=".", strip.white=TRUE)

length.ImgData0.lyso <- length(ImgData0$mean_lysosome)
require(relimp)
options(max.print = 99999999)
#showData(ImgData0[!complete.cases(ImgData0),], placement='-0-0', maxwidth=500, maxheight=200, suppress.X11.warnings=FALSE)
length.na.maxDiameter <- length(ImgData0$maxDiameter[is.na(ImgData0$maxDiameter)])
ImgData1 <- na.omit(ImgData0)
length.ImgData1.lyso <- length(ImgData1$mean_lysosome)

ImgData1$condition[ImgData1$coverslipNr >= 1 & ImgData1$coverslipNr <= 3] <- "live"
ImgData1$condition[ImgData1$coverslipNr >= 4 & ImgData1$coverslipNr <= 6] <- "dead"
ImgData1$condition[ImgData1$coverslipNr >= 7 & ImgData1$coverslipNr <= 9] <- "zmp1"
ImgData1$condition <- factor(ImgData1$condition, levels = c("live", "dead", "zmp1"))
condition.freq1 <- table(ImgData1$condition)

ImgData1$log.mean_macrophage <- log(ImgData1$mean_macrophage)
ImgData1$log.mean_lysosome <- log(ImgData1$mean_lysosome)

ImgData2 <- subset(ImgData1, physicalSize < 3)
condition.freq2 <- table(ImgData2$condition)


#hist(ImgData2$physicalSize, freq=TRUE, breaks=100, xlab="physical Size", ylab="Frequency", main="Frequency Histogram of Age")
#plot(ImgData2$physicalSize, ImgData2$log.mean_lysosome)
#plot(ImgData2$physicalSize, ImgData2$log.mean_lysosome)
#plot(ImgData2$physicalSize, ImgData2$roundness)
#plot(ImgData2$roundness, ImgData2$log.mean_lysosome)
#cor(ImgData2$physicalSize, ImgData2$roundness)
#showData(ImgData1, placement='-0-0', maxwidth=500, maxheight=200, suppress.X11.warnings=FALSE)

ImgData3 <- subset(ImgData2, physicalSize > 1 & log.mean_lysosome > 0)
condition.freq3 <- table(ImgData3$condition)

require(Rcmdr)
png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData1-plotMeans-lyso.png")
plotMeans(ImgData1$log.mean_lysosome, ImgData1$condition, error.bars="se", xlab="Condition", ylab="log of mean_lysosome values", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")))
dev.off()
png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData2-plotMeans-lyso.png")
plotMeans(ImgData2$log.mean_lysosome, ImgData2$condition, error.bars="se", xlab="Condition", ylab="log of mean_lysosome values", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")))
dev.off()
png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData3-plotMeans-lyso.png")
plotMeans(ImgData3$log.mean_lysosome, ImgData3$condition, error.bars="se", xlab="Condition", ylab="log of mean_lysosome values", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")))
dev.off()

png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData1-boxplot-lyso.png")
boxplot(log.mean_lysosome ~ condition, data=ImgData1, xlab = "Condition", ylab="log of mean_lysosome values", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")))
dev.off()
png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData2-boxplot-lyso.png")
boxplot(log.mean_lysosome ~ condition, data=ImgData2, xlab = "Condition", ylab="log of mean_lysosome values", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")))
dev.off()
png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData3-boxplot-lyso.png")
boxplot(log.mean_lysosome ~ condition, data=ImgData3, xlab = "Condition", ylab="log of mean_lysosome values", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")))
dev.off()

png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData1-plotMeans-roundness.png")
plotMeans(ImgData1$roundness, ImgData1$condition, error.bars="se", xlab="Condition", ylab="roundness", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")))
dev.off()
png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData2-plotMeans-roundness.png")
plotMeans(ImgData2$roundness, ImgData2$condition, error.bars="se", xlab="Condition", ylab="roundness", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")))
dev.off()
png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData3-plotMeans-roundness.png")
plotMeans(ImgData3$roundness, ImgData3$condition, error.bars="se", xlab="Condition", ylab="roundness", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")))
dev.off()

png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData1-plotMeans-size.png")
plotMeans(ImgData1$physicalSize, ImgData1$condition, error.bars="se", xlab="Condition", ylab="physicalSize", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")))
dev.off()
png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData2-plotMeans-size.png")
plotMeans(ImgData2$physicalSize, ImgData2$condition, error.bars="se", xlab="Condition", ylab="physicalSize", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")))
dev.off()
png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData3-plotMeans-size.png")
plotMeans(ImgData3$physicalSize, ImgData3$condition, error.bars="se", xlab="Condition", ylab="physicalSize", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")))
dev.off()

png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData1-plot-size-lyso.png")
plot(ImgData1$physicalSize, ImgData1$log.mean_lysosome, col=c("red","blue","green")[ImgData1$condition], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")))
dev.off()
png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData2-plot-size-lyso.png")
plot(ImgData2$physicalSize, ImgData2$log.mean_lysosome, col=c("red","blue","green")[ImgData2$condition], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")))
dev.off()
png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData3-plot-size-lyso.png")
plot(ImgData3$physicalSize, ImgData3$log.mean_lysosome, col=c("red","blue","green")[ImgData3$condition], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")))
dev.off()

png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData1-plot-roundness-log.mean_lysosome.png")
plot(ImgData1$roundness, ImgData1$log.mean_lysosome, col=c("red","blue","green")[ImgData1$condition])
dev.off()
png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData2-plot-roundness-log.mean_lysosome.png")
plot(ImgData2$roundness, ImgData2$log.mean_lysosome, col=c("red","blue","green")[ImgData2$condition])
dev.off()
png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData3-plot-roundness-log.mean_lysosome.png")
plot(ImgData3$roundness, ImgData3$log.mean_lysosome, col=c("red","blue","green")[ImgData3$condition])
dev.off()

png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData1-plot-physicalSize-roundness.png")
plot(ImgData1$physicalSize, ImgData1$roundness, col=c("red","blue","green")[ImgData1$condition])
dev.off()
png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData2-plot-physicalSize-roundness.png")
plot(ImgData2$physicalSize, ImgData2$roundness, col=c("red","blue","green")[ImgData2$condition])
dev.off()
png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData3-plot-physicalSize-roundness.png")
plot(ImgData3$physicalSize, ImgData3$roundness, col=c("red","blue","green")[ImgData3$condition])
dev.off()


length.ImgData0.lyso
length.na.maxDiameter
length.ImgData1.lyso
length.ImgData0.lyso - length.ImgData1.lyso
levels(ImgData1$condition)
names(ImgData1)
summary(ImgData1$physicalSize)

condition.freq1
sum(condition.freq1)
condition.freq2
sum(condition.freq2)
condition.freq3
sum(condition.freq3)

by(ImgData1$log.mean_lysosome, ImgData1$condition, summary)
by(ImgData1$log.mean_lysosome, ImgData1$condition, sd)
by(ImgData2$log.mean_lysosome, ImgData2$condition, summary)
by(ImgData2$log.mean_lysosome, ImgData2$condition, sd)
by(ImgData3$log.mean_lysosome, ImgData3$condition, summary)
by(ImgData3$log.mean_lysosome, ImgData3$condition, sd)

cor.matrix <- cor(ImgData1[, c("physicalSize", "maxDiameter", "minDiameter", "roundness", "log.mean_lysosome")])
round(cor.matrix, 2)
cor.matrix <- cor(ImgData2[, c("physicalSize", "maxDiameter", "minDiameter", "roundness", "log.mean_lysosome")])
round(cor.matrix, 2)
cor.matrix <- cor(ImgData3[, c("physicalSize", "maxDiameter", "minDiameter", "roundness", "log.mean_lysosome")])
round(cor.matrix, 2)


#require(VGAM)









if(FALSE) {
condition.ind1 <- which(ImgData1$condition=="live")
condition.ind2 <- which(ImgData1$condition=="dead")
condition.ind3 <- which(ImgData1$condition=="zmp1")
png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData1-plot-size-lyso-1-live.png")
plot(ImgData1$physicalSize[condition.ind1], ImgData1$log.mean_lysosome[condition.ind1], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")))
mtext("Live bacteria")
dev.off()
png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData1-plot-size-lyso-2-dead.png")
plot(ImgData1$physicalSize[condition.ind2], ImgData1$log.mean_lysosome[condition.ind2], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")))
mtext("Dead bacteria")
dev.off()
png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData1-plot-size-lyso-3-ko.png")
plot(ImgData1$physicalSize[condition.ind3], ImgData1$log.mean_lysosome[condition.ind3], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")))
mtext("Zmp1 ko bacteria")
dev.off()

condition.ind1 <- which(ImgData2$condition=="live")
condition.ind2 <- which(ImgData2$condition=="dead")
condition.ind3 <- which(ImgData2$condition=="zmp1")
png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData2-plot-size-lyso-1-live.png")
plot(ImgData2$physicalSize[condition.ind1], ImgData2$log.mean_lysosome[condition.ind1], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")))
mtext("Live bacteria")
dev.off()
png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData2-plot-size-lyso-2-dead.png")
plot(ImgData2$physicalSize[condition.ind2], ImgData2$log.mean_lysosome[condition.ind2], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")))
mtext("Dead bacteria")
dev.off()
png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData2-plot-size-lyso-3-ko.png")
plot(ImgData2$physicalSize[condition.ind3], ImgData2$log.mean_lysosome[condition.ind3], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")))
mtext("Zmp1 ko bacteria")
dev.off()

condition.ind1 <- which(ImgData3$condition=="live")
condition.ind2 <- which(ImgData3$condition=="dead")
condition.ind3 <- which(ImgData3$condition=="zmp1")
png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData3-plot-size-lyso-1-live.png")
plot(ImgData3$physicalSize[condition.ind1], ImgData3$log.mean_lysosome[condition.ind1], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")))
mtext("Live bacteria")
dev.off()
png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData3-plot-size-lyso-2-dead.png")
plot(ImgData3$physicalSize[condition.ind2], ImgData3$log.mean_lysosome[condition.ind2], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")))
mtext("Dead bacteria")
dev.off()
png(filename="/home/mmeuli/batch/out/1-1-A_result-ImgData3-plot-size-lyso-3-ko.png")
plot(ImgData3$physicalSize[condition.ind3], ImgData3$log.mean_lysosome[condition.ind3], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")))
mtext("Zmp1 ko bacteria")
dev.off()
}

#LANG=C sh -c 'R_DEFAULT_PACKAGES="$R_DEFAULT_PACKAGES Rcmdr" R "$@"'

#sudo apt-get install r-cran-rjava (required by XLConnect, which is needed for Rcmdr)
#install.packages("Rcmdr")
#install.packages("relimp")
#install.packages("VGAM")

setwd("/media/mmeuli/WD-HD-ext4/20160408_BCG_Pasteur-Aeras_in_THP-1/data-deconvolved/out")

min.size      <- 0
max.size      <- 1000000
min.red       <- -1000000000
min.roundness <- 0
min.ration    <- 0
max.ration    <- 100

ImgData0 <- read.table("1-1-A_result.txt", header=TRUE, sep="", na.strings="-nan", dec=".", strip.white=TRUE)

length.ImgData0.lyso <- length(ImgData0$lysosome)
require(relimp)
options(max.print = 99999999)
#showData(ImgData0[!complete.cases(ImgData0),], placement='-0-0', maxwidth=500, maxheight=200, suppress.X11.warnings=FALSE)
length.na.maxDiameter <- length(ImgData0$maxDiameter[is.na(ImgData0$maxDiameter)])
ImgData1 <- na.omit(ImgData0)
length.ImgData1.lyso <- length(ImgData1$lysosome)

ImgData1$condition[ImgData1$cNr >= 1 & ImgData1$cNr <= 3] <- "wt"
ImgData1$condition[ImgData1$cNr >= 4 & ImgData1$cNr <= 6] <- "dead"
ImgData1$condition[ImgData1$cNr >= 7 & ImgData1$cNr <= 9] <- "zmp1"
ImgData1$condition <- factor(ImgData1$condition, levels = c("wt", "dead", "zmp1"))
condition.freq1 <- table(ImgData1$condition)

ImgData1$log.macrophage <- log(ImgData1$macrophage)
ImgData1$log.lysosome <- log(ImgData1$lysosome)

ImgData2 <- subset(ImgData1, pSize > min.size & pSize < max.size & (maxDiameter/minDiameter) > min.ration & (maxDiameter/minDiameter) < max.ration)
condition.freq2 <- table(ImgData2$condition)

#hist(ImgData2$pSize, freq=TRUE, breaks=100, xlab="physical Size", ylab="Frequency", main="Frequency Histogram of Age")
#showData(ImgData1, placement='-0-0', maxwidth=500, maxheight=200, suppress.X11.warnings=FALSE)

ImgData3 <- subset(ImgData2, log.lysosome > min.red & roundness > min.roundness) 
condition.freq3 <- table(ImgData3$condition)

ImgData3w <- ImgData3
ImgData3w$lysosome <- ImgData3w$log.lysosome
ImgData3w$macrophage <- ImgData3w$log.macrophage
ImgData3w$pixels <- ImgData3w$condition <- ImgData3w$log.lysosome <- ImgData3w$log.macrophage <- NULL

write.table(format(ImgData3w, digits=3), file="1-1-A_result-ImgData3.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

C1D3 <- subset(ImgData3, condition == "wt")
C2D3 <- subset(ImgData3, condition == "dead")
C3D3 <- subset(ImgData3, condition == "zmp1")
wilcox.C1D3.C3D3 <- wilcox.test(C1D3$log.lysosome, C3D3$log.lysosome, paired=F)



require(Rcmdr)

png(filename="1-1-A_result-ImgData1-plotMeans-lyso.png")
plotMeans(ImgData1$log.lysosome, ImgData1$condition, error.bars="se", xlab="Condition", ylab="log of lysosome values", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData1")
dev.off()
png(filename="1-1-A_result-ImgData2-plotMeans-lyso.png")
plotMeans(ImgData2$log.lysosome, ImgData2$condition, error.bars="se", xlab="Condition", ylab="log of lysosome values", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData2")
dev.off()
png(filename="1-1-A_result-ImgData3-plotMeans-lyso.png")
plotMeans(ImgData3$log.lysosome, ImgData3$condition, error.bars="se", xlab="Condition", ylab="log of lysosome values", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData3")
dev.off()

png(filename="1-1-A_result-ImgData1-boxplot-lyso.png")
boxplot(log.lysosome ~ condition, data=ImgData1, xlab = "Condition", ylab="log of lysosome values", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData1")
dev.off()
png(filename="1-1-A_result-ImgData2-boxplot-lyso.png")
boxplot(log.lysosome ~ condition, data=ImgData2, xlab = "Condition", ylab="log of lysosome values", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData2")
dev.off()
png(filename="1-1-A_result-ImgData3-boxplot-lyso.png")
boxplot(log.lysosome ~ condition, data=ImgData3, xlab = "Condition", ylab="log of lysosome values", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData3")
dev.off()

png(filename="1-1-A_result-ImgData1-plotMeans-roundness.png")
plotMeans(ImgData1$roundness, ImgData1$condition, error.bars="se", xlab="Condition", ylab="roundness", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData1")
dev.off()
png(filename="1-1-A_result-ImgData2-plotMeans-roundness.png")
plotMeans(ImgData2$roundness, ImgData2$condition, error.bars="se", xlab="Condition", ylab="roundness", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData2")
dev.off()
png(filename="1-1-A_result-ImgData3-plotMeans-roundness.png")
plotMeans(ImgData3$roundness, ImgData3$condition, error.bars="se", xlab="Condition", ylab="roundness", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData3")
dev.off()

png(filename="1-1-A_result-ImgData1-plotMeans-size.png")
plotMeans(ImgData1$pSize, ImgData1$condition, error.bars="se", xlab="Condition", ylab="pSize", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData1")
dev.off()
png(filename="1-1-A_result-ImgData2-plotMeans-size.png")
plotMeans(ImgData2$pSize, ImgData2$condition, error.bars="se", xlab="Condition", ylab="pSize", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData2")
dev.off()
png(filename="1-1-A_result-ImgData3-plotMeans-size.png")
plotMeans(ImgData3$pSize, ImgData3$condition, error.bars="se", xlab="Condition", ylab="pSize", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData3")
dev.off()

png(filename="1-1-A_result-ImgData1-plot-size-lyso.png")
plot(ImgData1$pSize, ImgData1$log.lysosome, col=c("red","blue","green")[ImgData1$condition], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData1")
dev.off()
png(filename="1-1-A_result-ImgData2-plot-size-lyso.png")
plot(ImgData2$pSize, ImgData2$log.lysosome, col=c("red","blue","green")[ImgData2$condition], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData2")
dev.off()
png(filename="1-1-A_result-ImgData3-plot-size-lyso.png")
plot(ImgData3$pSize, ImgData3$log.lysosome, col=c("red","blue","green")[ImgData3$condition], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData3")
dev.off()

png(filename="1-1-A_result-ImgData1-plot-roundness-log.lysosome.png")
plot(ImgData1$roundness, ImgData1$log.lysosome, col=c("red","blue","green")[ImgData1$condition], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData1")
dev.off()
png(filename="1-1-A_result-ImgData2-plot-roundness-log.lysosome.png")
plot(ImgData2$roundness, ImgData2$log.lysosome, col=c("red","blue","green")[ImgData2$condition], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData2")
dev.off()
png(filename="1-1-A_result-ImgData3-plot-roundness-log.lysosome.png")
plot(ImgData3$roundness, ImgData3$log.lysosome, col=c("red","blue","green")[ImgData3$condition], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData3")
dev.off()

png(filename="1-1-A_result-ImgData1-plot-pSize-roundness.png")
plot(ImgData1$pSize, ImgData1$roundness, col=c("red","blue","green")[ImgData1$condition], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData1")
dev.off()
png(filename="1-1-A_result-ImgData2-plot-pSize-roundness.png")
plot(ImgData2$pSize, ImgData2$roundness, col=c("red","blue","green")[ImgData2$condition], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData2")
dev.off()
png(filename="1-1-A_result-ImgData3-plot-pSize-roundness.png")
plot(ImgData3$pSize, ImgData3$roundness, col=c("red","blue","green")[ImgData3$condition], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData3")
dev.off()



png(filename="1-1-A-A-result-graphs.png", width=210, height=297, units="mm", res=300, pointsize=12)
par(mfrow=c(4,3))

plotMeans(ImgData1$log.lysosome, ImgData1$condition, error.bars="se", xlab="Condition", ylab="log of lysosome values", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData1")
plotMeans(ImgData2$log.lysosome, ImgData2$condition, error.bars="se", xlab="Condition", ylab="log of lysosome values", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData2")
plotMeans(ImgData3$log.lysosome, ImgData3$condition, error.bars="se", xlab="Condition", ylab="log of lysosome values", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData3")

boxplot(log.lysosome ~ condition, data=ImgData1, xlab = "Condition", ylab="log of lysosome values", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData1")
boxplot(log.lysosome ~ condition, data=ImgData2, xlab = "Condition", ylab="log of lysosome values", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData2")
boxplot(log.lysosome ~ condition, data=ImgData3, xlab = "Condition", ylab="log of lysosome values", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData3")

plotMeans(ImgData1$roundness, ImgData1$condition, error.bars="se", xlab="Condition", ylab="roundness", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData1")
plotMeans(ImgData2$roundness, ImgData2$condition, error.bars="se", xlab="Condition", ylab="roundness", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData2")
plotMeans(ImgData3$roundness, ImgData3$condition, error.bars="se", xlab="Condition", ylab="roundness", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData3")

plotMeans(ImgData1$pSize, ImgData1$condition, error.bars="se", xlab="Condition", ylab="pSize", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData1")
plotMeans(ImgData2$pSize, ImgData2$condition, error.bars="se", xlab="Condition", ylab="pSize", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData2")
plotMeans(ImgData3$pSize, ImgData3$condition, error.bars="se", xlab="Condition", ylab="pSize", main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData3")

dev.off()


png(filename="1-1-A-A-result-plots.png", width=210, height=297, units="mm", res=300, pointsize=12)
par(mfrow=c(3,3))

plot(ImgData1$pSize, ImgData1$log.lysosome, col=c("red","blue","green")[ImgData1$condition], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData1")
plot(ImgData2$pSize, ImgData2$log.lysosome, col=c("red","blue","green")[ImgData2$condition], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData2")
plot(ImgData3$pSize, ImgData3$log.lysosome, col=c("red","blue","green")[ImgData3$condition], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData3")

plot(ImgData1$roundness, ImgData1$log.lysosome, col=c("red","blue","green")[ImgData1$condition], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData1")
plot(ImgData2$roundness, ImgData2$log.lysosome, col=c("red","blue","green")[ImgData2$condition], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData2")
plot(ImgData3$roundness, ImgData3$log.lysosome, col=c("red","blue","green")[ImgData3$condition], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData3")

plot(ImgData1$pSize, ImgData1$roundness, col=c("red","blue","green")[ImgData1$condition], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData1")
plot(ImgData2$pSize, ImgData2$roundness, col=c("red","blue","green")[ImgData2$condition], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData2")
plot(ImgData3$pSize, ImgData3$roundness, col=c("red","blue","green")[ImgData3$condition], main=expression(paste("Colocalization of BCG Pasteur-Aeras ",Delta,"zmp1")), sub="ImgData3")

dev.off()


options(width=as.integer(160))
sink("1-1-A-A-result-text.txt")

cat("Wilcox.test with ImgData3 and C1 and C3:\n")
wilcox.C1D3.C3D3
cat("\n\n")
cat("Number of objects in ImgData0:           \t", length.ImgData0.lyso, "\n")
cat("Objects with maxDiameter not a number:   \t",length.na.maxDiameter, "\n")
cat("Number of objects in ImgData1:           \t", length.ImgData1.lyso, "\n")
omitted.from.ImgData0 <- length.ImgData0.lyso - length.ImgData1.lyso
cat("Number of objects omitted from ImgData0: \t", omitted.from.ImgData0, "\n\n")
cat("Minimal size:                \t", min.size,      "\n")
cat("Maximal size:                \t", max.size,      "\n")
cat("Minimal red:                 \t", min.red,       "\n")
cat("Minimal roundness:           \t", min.roundness, "\n")
cat("Minimal ratio of diameters   \t", min.ratio,     "\n")
cat("Maximal ratio of diameters   \t", max.ratio,     "\n")
cat("Levels of ImgData1$condition:\t", levels(ImgData1$condition), "\n\n")
cat("ImgData2 <- subset(ImgData1, pSize > min.size & pSize < max.size & (maxDiameter/minDiameter) > min.ration & (maxDiameter/minDiameter) < max.ration)\n\n")
cat("ImgData3 <- subset(ImgData2, log.lysosome > min.red & roundness > min.roundness)\n\n") 
cat("Names of ImgData1:\n", names(ImgData1), "\n\n")
cat("Names of ImgData2:\n", names(ImgData2), "\n\n")
cat("Names of ImgData3:\n", names(ImgData3), "\n\n")
cat("\nSummary of ImgData1$pSize:\n")
summary(ImgData1$pSize)
cat("\nSummary of ImgData2$pSize:\n")
summary(ImgData2$pSize)
cat("\nSummary of ImgData3$pSize:\n")
summary(ImgData3$pSize)
cat("\nFrequencies of conditions in ImgData1: (sum:", sum(condition.freq1), ")")
condition.freq1
cat("\nFrequencies of conditions in ImgData2: (sum:", sum(condition.freq2), ")")
condition.freq2
cat("\nFrequencies of conditions in ImgData3: (sum:", sum(condition.freq3), ")")
condition.freq3
cat("\n\n")
cat("Summaries of ImgData1$log.lysosome for the conditions:\n\n")
by(ImgData1$log.lysosome, ImgData1$condition, summary)
cat("\n")
cat("Stardard deviations of ImgData1$log.lysosome for the conditions:\n\n")
by(ImgData1$log.lysosome, ImgData1$condition, sd)
cat("\n")
cat("Summaries of ImgData2$log.lysosome for the conditions:\n\n")
by(ImgData2$log.lysosome, ImgData2$condition, summary)
cat("\n")
cat("Stardard deviations of ImgData2$log.lysosome for the conditions:\n\n")
by(ImgData2$log.lysosome, ImgData2$condition, sd)
cat("\n")
cat("Summaries of ImgData3$log.lysosome for the conditions:\n\n")
by(ImgData3$log.lysosome, ImgData3$condition, summary)
cat("\n")
cat("Stardard deviations of ImgData3$log.lysosome for the conditions:\n\n")
by(ImgData3$log.lysosome, ImgData3$condition, sd)
cat("\n\n")

cat("Correlation matrix of ImgData1:\n")
cor.matrix <- cor(ImgData1[, c("pSize", "maxDiameter", "minDiameter", "roundness", "log.lysosome")])
round(cor.matrix, 2)
cat("\n")
cat("Correlation matrix of ImgData2:\n")
cor.matrix <- cor(ImgData2[, c("pSize", "maxDiameter", "minDiameter", "roundness", "log.lysosome")])
round(cor.matrix, 2)
cat("\n")
cat("Correlation matrix of ImgData3:\n")
cor.matrix <- cor(ImgData3[, c("pSize", "maxDiameter", "minDiameter", "roundness", "log.lysosome")])
round(cor.matrix, 2)
cat("\n")

sink()


#could be done: use a threshold to exclude bacteria wich are probably outside the macrophages
#threshold <- mean(ImgData3$log.macrophage)-(2*sd(ImgData3$log.lysosome))







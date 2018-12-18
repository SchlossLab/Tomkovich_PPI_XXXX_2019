#Load packages
library(tidyverse)
library(dtw)
library(vegan)

#Define input files
cfu <- 'data/08_20_18_Omep.#2/08_20_18_Omep._#2_exp.xlsx'

# Define output files
plot_file <- 'results/figures/ppi2_cfu.pdf'

#Load in data
ppi2_data <- readxl::read_excel("data/08_20_18_Omep.#2/08_20_18_Omep._#2_exp.xlsx", sheet = "cfu_final")

#Create data frame of CFU data by selecting group column & columns ending with CFU/g
ppi2_cfu <- select(ppi2_data, Group, ends_with("CFU/g")) 
#make sure Group column is a factor vector and all CFU/g are treated as numeric vectors:
ppi2_cfu$Group <- as.factor(ppi2_cfu$Group)
str(ppi2_cfu)

#Format CFU over time
ppi2_cfu[,2:ncol(ppi2_cfu)] <- log10(ppi2_cfu[,2:ncol(ppi2_cfu)] + 1)
cfu_median <- aggregate(ppi2_cfu[,2:ncol(ppi2_cfu)], by=list(ppi2_cfu$Group), FUN=quantile, probs=0.5, na.rm = T)
rownames(cfu_median) <- cfu_median$Group.1
cfu_median$Group.1 <- NULL
cfu_median <- as.data.frame(t(cfu_median))
cfu_q25 <- aggregate(ppi2_cfu[,2:ncol(ppi2_cfu)], by=list(ppi2_cfu$Group), FUN=quantile, probs=0.25, na.rm = T)
rownames(cfu_q25) <- cfu_q25$Group.1
cfu_q25$Group.1 <- NULL
cfu_q25 <- as.data.frame(t(cfu_q25))
cfu_q75 <- aggregate(ppi2_cfu[,2:ncol(ppi2_cfu)], by=list(ppi2_cfu$Group), FUN=quantile, probs=0.75, na.rm = T)
rownames(cfu_q75) <- cfu_q75$Group.1
cfu_q75$Group.1 <- NULL
cfu_q75 <- as.data.frame(t(cfu_q75))

#Plot of CFU over time
pdf(file=plot_file)
par(mar=c(4,4,1,1), las=1, mgp=c(2.5, 0.75, 0)) #http://rfunction.com/archives/1302 
plot(0, type='n', xlab='Days Post-Infection', ylab='Total cfu/g Feces', xaxt='n', yaxt='n', xlim=c(1,7), ylim=c(0,10))
lines(cfu_median$`C+`, lwd=2.5, col="blue", type='b', pch=19) 
segments(x0=c(1:7), y0=cfu_q25$`C+`, x1=c(1:7), y1=cfu_q75$`C+`, col="blue", lwd=2.5)
lines(cfu_median$`CO+`, lwd=2.5, col="grey3", type='b', pch=19) 
segments(x0=c(1:7), y0=cfu_q25$`CO+`, x1=c(1:7), y1=cfu_q75$`CO+`, col="grey3", lwd=2.5)
lines(cfu_median$`O+`, lwd=2.5, col="purple", type='b', pch=19) 
segments(x0=c(1:7), y0=cfu_q25$`O+`, x1=c(1:7), y1=cfu_q75$`O+`, col="purple", lwd=2.5)
abline(h=2, lwd=2, lty=2)
axis(side=1, at=c(0:11), labels=c(-1:10))
axis(side=2, at=seq(0,10,1), labels=c(0, parse(text=paste(rep(10,10), '^', seq(1,10,1), sep=''))), las=1)
legend(x=5, y=10, legend=c('Clindamycin', 'Clind. + Omep.', 'Omeprazole'),
       pch=16, col=c("blue", "grey3", "purple"), cex=0.9, pt.cex=1.5)
text(x=7, y=2.2, 'LOD', cex=0.9)
dev.off()


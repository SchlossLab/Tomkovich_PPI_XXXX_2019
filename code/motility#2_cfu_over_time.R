#Load packages
library(tidyverse)
library(dtw)
library(vegan)

#Define input files
cfu <- 'data/Motility #2/11_14_18_motility_exp_2.xlsx'

# Define output files
plot_file <- 'results/figures/motility2_cfu.pdf'

#Load in data
motility2_data <- readxl::read_excel("data/Motility #2/11_14_18_motility_exp_2.xlsx", sheet = "cfu_final")

#Create data frame of CFU data by selecting group column & columns ending with CFU/g
motility2_cfu <- select(motility2_data, Group, ends_with("CFU/g")) 
#make sure Group column is a factor vector and all CFU/g are treated as numeric vectors:
motility2_cfu$Group <- as.factor(motility2_cfu$Group)
motility2_cfu$`D0 C. difficile CFU/g` <- as.numeric(motility2_cfu$`D0 C. difficile CFU/g`)
motility2_cfu$`D3 C. difficile CFU/g` <- as.numeric(motility2_cfu$`D3 C. difficile CFU/g`)
motility2_cfu$`D5 C. difficile CFU/g` <- as.numeric(motility2_cfu$`D5 C. difficile CFU/g`)
str(motility2_cfu)


#Format CFU over time
motility2_cfu[,2:ncol(motility2_cfu)] <- log10(motility2_cfu[,2:ncol(motility2_cfu)] + 1)
cfu_median <- aggregate(motility2_cfu[,2:ncol(motility2_cfu)], by=list(motility2_cfu$Group), FUN=quantile, probs=0.5, na.rm = T)
rownames(cfu_median) <- cfu_median$Group.1
cfu_median$Group.1 <- NULL
cfu_median <- as.data.frame(t(cfu_median))
cfu_q25 <- aggregate(motility2_cfu[,2:ncol(motility2_cfu)], by=list(motility2_cfu$Group), FUN=quantile, probs=0.25, na.rm = T)
rownames(cfu_q25) <- cfu_q25$Group.1
cfu_q25$Group.1 <- NULL
cfu_q25 <- as.data.frame(t(cfu_q25))
cfu_q75 <- aggregate(motility2_cfu[,2:ncol(motility2_cfu)], by=list(motility2_cfu$Group), FUN=quantile, probs=0.75, na.rm = T)
rownames(cfu_q75) <- cfu_q75$Group.1
cfu_q75$Group.1 <- NULL
cfu_q75 <- as.data.frame(t(cfu_q75))

#Plot of CFU over time
pdf(file=plot_file)
par(mar=c(4,4,1,1), las=1, mgp=c(2.5, 0.75, 0)) #http://rfunction.com/archives/1302 
plot(0, type='n', xlab='Days Post-Infection', ylab='Total cfu/g Feces', xaxt='n', yaxt='n', xlim=c(1,6), ylim=c(0,10))
lines(cfu_median$`L`, lwd=2.5, col="turquoise4", type='b', pch=19) 
segments(x0=c(1:7), y0=cfu_q25$`L`, x1=c(1:7), y1=cfu_q75$`L`, col="turquoise4", lwd=2.5)
lines(cfu_median$`LC`, lwd=2.5, col="turquoise", type='b', pch=19) 
segments(x0=c(1:7), y0=cfu_q25$`LC`, x1=c(1:7), y1=cfu_q75$`LC`, col="turquoise", lwd=2.5)
lines(cfu_median$`P`, lwd=2.5, col="purple", type='b', pch=19) 
segments(x0=c(1:7), y0=cfu_q25$`P`, x1=c(1:7), y1=cfu_q75$`P`, col="purple", lwd=2.5)
lines(cfu_median$`PC`, lwd=2.5, col="mediumpurple1", type='b', pch=19) 
segments(x0=c(1:7), y0=cfu_q25$`PC`, x1=c(1:7), y1=cfu_q75$`PC`, col="mediumpurple1", lwd=2.5)
lines(cfu_median$`C`, lwd=2.5, col="grey3", type='b', pch=19) 
segments(x0=c(1:7), y0=cfu_q25$`C`, x1=c(1:7), y1=cfu_q75$`C`, col="grey3", lwd=2.5)
abline(h=2, lwd=2, lty=2)
axis(side=1, at=c(0:11), labels=c(-1:10))
axis(side=2, at=seq(0,10,1), labels=c(0, parse(text=paste(rep(10,10), '^', seq(1,10,1), sep=''))), las=1)
legend(x=4.5, y=10, legend=c('Loperamide', 'Loperamide + Clind.', 'Miralax', 'Miralax + Clind.', 'Clindamycin'),
       pch=16, col=c("turquoise4", "turquoise", "purple", "mediumpurple1", "grey3"), cex=0.9, pt.cex=1.5)
text(x=6, y=2.2, 'LOD', cex=0.9)
dev.off()

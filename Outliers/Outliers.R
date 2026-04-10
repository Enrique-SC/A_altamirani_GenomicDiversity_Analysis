library(pcadapt)
library(vegan)
library(qvalue)
library(vcfR)
library(dplyr)
library(tidyverse)

#VCF
Alt.vcf <-read.vcfR("alt.variantsprunnedready.org.noindel.final.chrom.vcf")

#PCADAPT
pca_genotype <- read.pcadapt("alt.variantsprunnedready.org.noindel.final.chrom.vcf", type = "vcf")
K <- 3
x <- pcadapt(pca_genotype, K = K, min.maf = 0.05)

padj <- p.adjust(x$pvalues,method="bonferroni")
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.05

outliers_pcadapt <- which(padj < alpha)

#Generar lista de SNPs a eliminar con bcftools
chr.pos <- as.data.frame(Alt.vcf@fix[,1:2])
chr.pos$number <- 1:nrow(chr.pos)

snp.pcadapt <- as.data.frame(outliers_pcadapt)

snps.neutral <- chr.pos %>% filter(!number %in% snp.pcadapt$outliers_pcadapt)
snps.outlier <- chr.pos %>% filter(number %in% snp.pcadapt$outliers_pcadapt)

write.csv(snps.outlier, "snp.pcadapt.csv", row.names = T)

#RDA
alta.gt <- read.table("alta.final.chrom.raw.raw", header=TRUE, sep=" ", row.names=1)
alta.gt <- alta.gt[,6:112399]

env.alta <- read.csv("metadata_alta_bd_final.csv")
rownames(env.alta)<- env.alta$sample.id
env.alta <- env.alta[3:9]
head(env.alta)
# Confirm that genotypes and environmental data are in the same order
rownames(alta.gt)<- rownames(env.alta)
identical(rownames(alta.gt), rownames(env.alta))

gt.env.rda <- rda(alta.gt ~ ., data=env.rda, scale=T)
gt.env.rda

load.rda <- scores(gt.env.rda, choices=c(1:3), display="species")  # Species scores for the first three constrained axes

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

cand1 <- outliers(load.rda[,1],3) # 64
cand2 <- outliers(load.rda[,2],3) # 563
cand3 <- outliers(load.rda[,3],3) # 1063

ncand <- length(cand1) + length(cand2) + length(cand3)

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

foo <- matrix(nrow=(ncand), ncol=7)  # 8 columns for 8 predictors
colnames(foo) <- c("Elevation", "WaterTemp","WaterPH","DO","MeanMaxTemp ","MeanMinTemp","Precipitation")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- alta.gt[,nam]
  foo[i,] <- apply(env.rda,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)  

cand$snp_1<-gsub(cand$snp, pattern = '\\.',replacement = "_")
cand.final <- separate(cand,snp_1,c('CHROM', 'Cop','POS','VAR'))

snps.outlier.pca <- chr.pos %>% filter(number %in% snp.pcadapt$outliers_pcadapt)
snps.outlier.rda <- chr.pos %>% filter(POS %in% cand.final$POS)

pos.snp.pca <- snps.outlier.pca$POS
pos.snp.rda <- snps.outlier.rda$POS

snps.outlier.final <- as.data.frame(rbind(snps.outlier.pca,snps.outlier.rda))

snps.outlier.final <- snps.outlier.final[!duplicated(snps.outlier.final$POS),]
snps.neutral.final <- chr.pos %>% filter(!number %in% snps.outlier.final$number)

write.csv(snps.neutral.final, "snp.neutrals.final1.csv", row.names = F)
write.csv(snps.outlier.final, "snp.outlier.final1.csv", row.names = F)
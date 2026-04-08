######################################################
###Figure 2 Script - Paper Soto-Cortés et al., 2026###
######################################################
#Loading packages

library(vcfR)
library(ggplot2)
library(tidyverse)
library(ggsci)
library(agricolae)
library(inbreedR)
library(dartR)
library(ggpubr)
library(nlme)
library(lmodel2)
library(FSA)


#Calculating multilocus heterozygosity (MLH)
# Eead vcf
vcf.alta <- read.vcfR("alta.neutral.final.chroms.vcf", verbose = FALSE )
# Extract genotypes
gt.alta <- extract.gt(vcf.alta)
# Transpose and data.frame
gt.alta <- as.data.frame(t(gt.alta), stringsAsFactors = FALSE)
# NA handling
gt.alta[gt.alta == "."] <- NA
# Convert to raw format 
alta_snp_genotypes <- convert_raw(gt.alta)


gen.alt <- as.data.frame(MLH(alta_snp_genotypes))
colnames(gen.alt) <- c("MLH")

### Estimating inbreeding coefficient 

# Fis was estimated with plink with command: ./plink --vcf alta.neutral.final.chroms.vcf --allow-extra-chr --het

div.stacks.alt <- read.table("plink.het.txt",header = T) #Previously calculated in plink
gen.alt$Fis <- div.stacks.alt$F

#Adding locality to the dataframe
localities <- read.table("pop.alt2_new.txt")
gen.alt$Locality <- localities$V2

#Scaled Mass Index (SMI) estimation

data_body <- read.csv("svl_weight_alta.csv")
SMAreg_opum<- lmodel2(log10(data_body$weight) ~ log10(data_body$svl), data=data_body, nperm=1000)

### Getting the parameters to apply the formula of Peig & Green (2009).

bSMA_m<-  2.466789 
Lmean_m<-mean(data_body$svl)

### To calculate the Scaled Mass Index for every frog and save it in an extra column in the datasheet.
data_body$SMI<-data_body$weight*(Lmean_m/data_body$svl)^bSMA_m
write.csv(data_body,"svl_weight_SMI_alta.csv",row.names = F)

#Adding SMI colum to gen.alt dataframe 

gen.alt$SMI <- data_body$SMI

#Remove outlier sample values 

#Outliers MLH values
gen.alt.out<- gen.alt %>%  filter(!row.names(gen.alt)=='AlTe02_AlTe02')
gen.alt.out<- gen.alt.out %>%  filter(!row.names(gen.alt.out)=='AlLS07_AlLS07')
gen.alt.out<- gen.alt.out %>%  filter(!row.names(gen.alt.out)=='AlO20_AlO20')
gen.alt.out<- gen.alt.out %>%  filter(!row.names(gen.alt.out)=='AlSe19_AlSe19')
#Outlier SMI value

gen.alt.out<- gen.alt.out %>%  filter(!row.names(gen.alt.out)=='AlTe08_AlTe08')

#Boxplots plotting 

#MHL boxplot
mlh.summarized = gen.alt.out %>% group_by(Locality) %>% summarize(MLH=max(MLH))
hsd.mlh=HSD.test(aov(MLH~Locality,data=gen.alt.out), "Locality", group=T)
hsd.mlh
hsd.df.mlh <- as.data.frame(hsd.mlh$groups)
hsd.df.mlh$sample <- row.names(hsd.df.mlh)
hsd.df.mlh <- hsd.df.mlh[order(hsd.df.mlh$sample),]

mlh.alt <- ggplot(gen.alt.out, aes(x=Locality, y=MLH, fill=Locality)) + 
  geom_boxplot(lwd=1,alpha=0.8) + 
  scale_fill_locuszoom() + theme_classic2() +
  geom_text(data=mlh.summarized,aes(x=Locality,y=0.007+MLH,label=hsd.df.mlh$groups),vjust=0,fontface="bold", size = 5)+
  theme(legend.position="none",
        legend.text = element_text(size=13), 
        axis.text.x=element_text(face = "bold",colour = "black", size = 20),
        axis.text.y=element_text(colour = "black", size = 18),
        axis.title.y=element_text(face = "bold",colour = "black", size = 20),
        axis.title.x=element_blank())+
  ylab("Multilocus Heterozygosity (MLH)")

#Fis boxplot
fis.summarized = gen.alt.out %>% group_by(Locality) %>% summarize(Fis=max(Fis))
hsd.fis=HSD.test(aov(Fis~Locality,data=gen.alt.out), "Locality", group=T)
hsd.fis
hsd.df.fis <- as.data.frame(hsd.fis$groups)
hsd.df.fis$sample <- row.names(hsd.df.fis)
hsd.df.fis <- hsd.df.fis[order(hsd.df.fis$sample),]

fis.alt <- ggplot(gen.alt.out, aes(x=Locality, y=Fis, fill=Locality)) + 
  geom_boxplot(lwd=1,alpha=0.8) +
  scale_fill_locuszoom() + theme_classic2()  +
  geom_text(data=fis.summarized,aes(x=Locality,y=0.007+Fis,label=hsd.df.fis$groups),vjust=0,fontface="bold", size = 5)+
  theme(legend.position="none",
        axis.text.x=element_text(face = "bold",colour = "black", size = 20),
        axis.text.y=element_text(colour = "black", size = 18),
        axis.title.y=element_text(face = "bold",colour = "black", size = 20),
        axis.title.x=element_blank())+
  ylab("Inbreeding coefficient (Fis)")
fis.alt

fig2 <- ggarrange(mlh.alt,fis.alt, ncol=2, nrow = 1,
                  labels = c("A","B"), font.label = list(size = 15, face = "bold"),
                  common.legend = F)
fig2
ggsave("Fig2.pdf", fig2, width = 13, height =7, device = "pdf",
       limitsize = F, units = "in", path = "./Plots/Finales")

#Kruskal Wallis tests 

kw.mlh <- kruskal.test(MLH ~ Locality, data = gen.alt.out)
kw.mlh #Kruskal-Wallis chi-squared = 46.296, df = 4, p-value = 2.137e-09

kw.fis <- kruskal.test(Fis ~ Locality, data = gen.alt.out)
kw.fis #Kruskal-Wallis chi-squared = 11.055, df = 4, p-value = 0.02596


#Dunn tests

dunn.mlh <- dunnTest(MLH ~ Locality, data = gen.alt.out)
dunn.mlh

dunn_df_mlh <- as.data.frame(dunn.mlh$res)
write.csv(dunn_df_mlh, "dunn_df_mlh.csv",row.names = F)

dunn.fis <- dunnTest(Fis ~ Locality, data = gen.alt.out)
dunn.fis

dunn_df_fis <- as.data.frame(dunn.fis$res)
dunn_df_fis

write.csv(dunn_df_fis, "dunn_df_fis.csv",row.names = F)

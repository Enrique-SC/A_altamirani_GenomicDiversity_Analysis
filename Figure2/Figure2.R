#Figure 2 Script - Paper Soto-Cortés et al..


library(ggplot2)
library(tidyverse)
library(ggsci)
library(agricolae)
library(reshape2)
library("LEA")
library(inbreedR)
library(dartR)
library(ggpubr)
library(gwscaR)
library(vegan)
library(ggmanh)
library(poppr)
library(lfmm)
library(smplot2)

gen.alt <- read.csv("Div.genet.alt.org.neutros_rmoutliers.csv")

#Boxplots
#MHL
mlh.summarized = gen.alt %>% group_by(Locality) %>% summarize(MLH=max(MLH))
hsd.mlh=HSD.test(aov(MLH~Locality,data=gen.alt), "Locality", group=T)
hsd.mlh
hsd.df.mlh <- as.data.frame(hsd.mlh$groups)
hsd.df.mlh$sample <- row.names(hsd.df.mlh)
hsd.df.mlh <- hsd.df.mlh[order(hsd.df.mlh$sample),]

mlh.alt <- ggplot(gen.alt, aes(x=Locality, y=MLH, fill=Locality)) + 
  geom_boxplot(lwd=1,alpha=0.8) + 
  scale_fill_locuszoom() + theme_classic2() +
  geom_text(data=mlh.summarized,aes(x=Locality,y=0.007+MLH,label=hsd.df.mlh$groups),vjust=0,fontface="bold", size = 5)+
  theme(legend.position="none",
        legend.text = element_text(size=13), 
        axis.text.x=element_text(face = "bold",colour = "black", size = 20),
        axis.text.y=element_text(colour = "black", size = 18),
        axis.title.y=element_text(face = "bold",colour = "black", size = 20),
        axis.title.x=element_blank())+
  ylab("Multilocus Heterocigozity (MLH)")
mlh.alt
ggsave("mlh.viol.png", mlh.alt, width = 14, height = 8,
       limitsize = F, units = "in", path = "./Plots/Finales/")

#Fis
fis.summarized = gen.alt %>% group_by(Locality) %>% summarize(Fis=max(Fis))
hsd.fis=HSD.test(aov(Fis~Locality,data=gen.alt), "Locality", group=T)
hsd.fis
hsd.df.fis <- as.data.frame(hsd.fis$groups)
hsd.df.fis$sample <- row.names(hsd.df.fis)
hsd.df.fis <- hsd.df.fis[order(hsd.df.fis$sample),]

fis.alt <- ggplot(gen.alt, aes(x=Locality, y=Fis, fill=Locality)) + 
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
ggsave("he.viol.png", he.alt, width = 10, height = 8,
       limitsize = F, units = "in", path = "./Plots/")



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

mean(gen.alt$MLH) #0.40
sd(gen.alt$MLH) #0.013

mean(gen.alt$Fis) #0.02
sd(gen.alt$Fis) #0.066

#Dunn test fig 2


dunn.mlh <- dunnTest(MLH ~ Locality, data = gen.alt.out)
dunn.mlh

dunn_df_mlh <- as.data.frame(dunn.mlh$res)
write.csv(dunn_df_mlh, "dunn_df_mlh.csv",row.names = F)

dunn.fis <- dunnTest(Fis ~ Locality, data = gen.alt.out)
dunn.fis

dunn_df_fis <- as.data.frame(dunn.fis$res)
dunn_df_fis
write.csv(dunn_df_fis, "dunn_df_fis.csv",row.names = F)

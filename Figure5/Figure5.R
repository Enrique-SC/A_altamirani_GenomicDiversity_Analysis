######################################################
###Figure 5 Script - Paper Soto-Cortés et al., 2026###
######################################################

library(vegan)
library(ggplo2)


#Generate .raw format of vcfile to perform RDA using vcftools and plink:
#vcftools --vcf alt.variantsprunnedready.org.noindel.final.chrom.vcf --plink-tped --out alta.org.fb
#./plink --tped alta.org.fb.tped --tfam alta.org.fb.tfam --recodeA --out alta.org.fb
#./plink --vcf alt.variantsprunnedready.org.noindel.final.chrom.vcf --allow-extra-chr --recode A --out alta.final.chrom.raw



alta.gt <- read.table("alta.final.chrom.raw.raw", header=TRUE, sep=" ", row.names=1)
dim(alta.gt)
alta.gt <- alta.gt[,6:112399]

env.alta <- read.csv("metadata_alta_bd_final.csv")
rownames(env.alta)<- env.alta$sample.id
env.alta <- env.alta[3:9]


# Confirm that genotypes and environmental data are in the same order
rownames(alta.gt)<- rownames(env.alta)
identical(rownames(alta.gt), rownames(env.alta))
## [1] TRUE

pairs.panels(env.alta, scale=T)

env.rda <- subset(env.alta, select=c(Elevation, WaterpH,MinTem,Precipitation,DO))

pairs.panels(env.rda, scale=T)
str(env.rda)

#RDA model
gt.env.rda <- rda(alta.gt ~ ., data=env.rda, scale=T)

#
perma.rda.all <-anova.cca(gt.env.rda, by = "terms", permu=999)
perma.rda.all
save(perma.rda.all,file = "perma.rda.all.rda")


#Plotting RDA

site_scores <- as.data.frame(scores(gt.env.rda, display = "sites"))
site_scores$Locality <- env.df$Locality
species_scores <- as.data.frame(scores(gt.env.rda, display = "species"))
vec_scores <- as.data.frame(scores(gt.env.rda, display = "bp"))
vec_scores$F_values <- c("7.6","3.4","3","5","1.9","1.2","1.1")
vec_scores$F_values <- as.numeric(vec_scores$F_values)

adj.vec <- 6
adj.txt <- 6

rda.gen <- ggplot(site_scores, aes(x=RDA1, y=RDA2)) +
  geom_point(size=5, pch=21, colour="black", aes(fill=Locality),stroke=1.5) +
  geom_segment(data=vec_scores, inherit.aes=F, linewidth = 1,
               mapping=aes(x=0, y=0, xend=adj.vec*RDA1, yend=adj.vec*RDA2,colour=F_values), 
               arrow=arrow(length=unit(0.4, 'cm'))) + scale_color_gradient(low = "lightgrey", high="black")+
  geom_text(data=vec_scores, inherit.aes=F, size=5, fontface="bold",
            mapping=aes(x=adj.txt*RDA1, y=adj.txt*RDA2, 
                        label=c('Elevation','WaterTemp','WaterPH','DO',"MeanMaxTemp", "MeanMinTemp", "Precipitation")))+
  theme_bw() + scale_fill_locuszoom() +
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 14,face = "bold"),
        axis.text.y=element_text(colour = "black", size = 14),
        axis.text.x=element_text(colour = "black", size = 14),
        axis.title.y=element_text(face = "bold",colour = "black", size = 15),
        axis.title.x=element_text(face = "bold",colour = "black", size = 15))
rda.gen$labels$x <- ord_labels(gt.env.rda)[1]
rda.gen$labels$y <- ord_labels(gt.env.rda)[2]

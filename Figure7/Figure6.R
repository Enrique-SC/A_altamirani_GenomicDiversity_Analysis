######################################################
###Figure 6 Script - Paper Soto-Cortés et al., 2026###
######################################################

#Loading packages

library(tidyverse)
library(ggvenn)
library(scales)
library(ggplo2)
library(ggpubr)

#Loding dataframe containing the SNPs significanly associated with abiotic and biotic factors. To see how this dataframe was generated please go to the Models scrip from this repository

lfmm.env.qv<- read.csv("lfmm.env.qv.filtered.csv")
row.names(lfmm.env.qv)<- lfmm.env.qv$X 

#Extracting SNPs associated to each particular factor

lfmm.env.qv.sub.elev <- lfmm.env.qv %>% filter(Elevation < 0.05)
lfmm.env.qv.sub.elev.only<- as.data.frame(lfmm.env.qv.sub.elev[,1])
colnames(lfmm.env.qv.sub.elev.only)<- "Elevation"
lfmm.env.qv.sub.elev.only$SNP<- row.names(lfmm.env.qv.sub.elev)

lfmm.env.qv.sub.temp <- lfmm.env.qv %>% filter(WaterTemp < 0.05)
lfmm.env.qv.sub.temp.only<- as.data.frame(lfmm.env.qv.sub.temp[,2])
colnames(lfmm.env.qv.sub.temp.only)<- "Water temperature"
lfmm.env.qv.sub.temp.only$SNP<- row.names(lfmm.env.qv.sub.temp)

lfmm.env.qv.sub.do <- lfmm.env.qv %>% filter(DO < 0.05)
lfmm.env.qv.sub.do.only<- as.data.frame(lfmm.env.qv.sub.do[,3])
colnames(lfmm.env.qv.sub.do.only)<- "Dissoved oxygen"
lfmm.env.qv.sub.do.only$SNP<- row.names(lfmm.env.qv.sub.do)

lfmm.biotic.qv.sub.bd <- lfmm.env.qv %>% filter(BdLog < 0.05)
lfmm.biotic.qv.sub.bd.only <- as.data.frame(lfmm.biotic.qv.sub.bd[,4])
colnames(lfmm.biotic.qv.sub.bd.only)<- "Bd intection"
lfmm.biotic.qv.sub.bd.only$SNP<- row.names(lfmm.biotic.qv.sub.bd)

lfmm.biotic.qv.sub.body <- lfmm.env.qv %>% filter(BodyCondition < 0.05)
lfmm.biotic.qv.sub.body.only<- as.data.frame(lfmm.biotic.qv.sub.body[,2])
colnames(lfmm.biotic.qv.sub.body.only)<- "Body condition"
lfmm.biotic.qv.sub.body.only$SNP<- row.names(lfmm.biotic.qv.sub.body)

#Generating character vectors 

bd.lfmm <- lfmm.biotic.qv.sub.bd.only$SNP
body.lfmm <- lfmm.biotic.qv.sub.body.only$SNP
elev.lfmm <- lfmm.env.qv.sub.elev.only$SNP
temp.lfmm <- lfmm.env.qv.sub.temp.only$SNP
do.lfmm <- lfmm.env.qv.sub.do.only$SNP


#Venn diagram plotting 

lfmmm = list(Elevation=elev.lfmm,Temp=temp.lfmm,DO=do.lfmm,Bd=bd.lfmm,BC=body.lfmm)

plot.venn2<-ggvenn(lfmmm, fill_color = c("#8569D5", "#5E738F","#D1A33D", "#8A7C64","#599861"),
                   fill_alpha = 0.8, set_name_size = 7, text_size = 6)
plot.venn2


#Barplot plotting

#Loading dataframe with annotations for each SNPs

snps.lfmm.ann <- read.csv("lfmm.env.qv.filtered.annotated.csv")



snps.lfmm.ann.genic <- snps.lfmm.ann[snps.lfmm.ann$Region_type == 'Genic region',]
snps.lfmm.ann.Bd <- snps.lfmm.ann[snps.lfmm.ann$BdLog < 0.05, ]
snps.lfmm.ann.body <- snps.lfmm.ann[snps.lfmm.ann$BodyCondition < 0.05, ]
snps.lfmm.ann.temp <- snps.lfmm.ann[snps.lfmm.ann$WaterTemp < 0.05, ]
snps.lfmm.ann.do <- snps.lfmm.ann[snps.lfmm.ann$DO < 0.05, ]
snps.lfmm.ann.elev <- snps.lfmm.ann[snps.lfmm.ann$Elevation < 0.05, ]

#Generating character vectors 
setBd <- snps.lfmm.ann.Bd$SNP
setBody <- snps.lfmm.ann.body$SNP
setTemp <- snps.lfmm.ann.temp$SNP
setDO <- snps.lfmm.ann.do$SNP
setElev <- snps.lfmm.ann.elev$SNP


## Identifyting the SNPs uniquely associated to each factor

# Find the union of all other sets
other_sets_union_Temp <- union(union(union(setBd, setBody), setDO), setElev)
other_sets_union_DO <- union(union(union(setBd, setBody), setTemp), setElev)
other_sets_union_Elev <- union(union(union(setBd, setBody), setDO), setTemp)
other_sets_union_Bd <- union(union(union(setTemp, setBody), setDO), setElev)
other_sets_union_body <- union(union(union(setBd, setTemp), setDO), setElev)

# Find elements in set1 that are NOT in the union of other sets
unique_to_setTemp <- setdiff(setTemp, other_sets_union_Temp)
unique_to_setDO <- setdiff(setDO, other_sets_union_DO)
unique_to_setElev <- setdiff(setElev, other_sets_union_Elev)
unique_to_setBody <- setdiff(setBody, other_sets_union_body)
unique_to_setBd <- setdiff(setBd, other_sets_union_Bd)


#Keeping unique SNPs for each factor in dataframe

snps.lfmm.ann.Bd.unique <- snps.lfmm.ann.Bd[snps.lfmm.ann.Bd$SNP %in% unique_to_setBd, ]
snps.lfmm.ann.body.unique <-  snps.lfmm.ann.body[snps.lfmm.ann.body$SNP %in% unique_to_setBody, ]
snps.lfmm.ann.temp.unique <-  snps.lfmm.ann.temp[snps.lfmm.ann.temp$SNP %in% unique_to_setTemp, ]
snps.lfmm.ann.do.unique <-  snps.lfmm.ann.do[snps.lfmm.ann.do$SNP %in% unique_to_setDO, ]
snps.lfmm.ann.elev.unique <-  snps.lfmm.ann.elev[snps.lfmm.ann.elev$SNP %in% unique_to_setElev, ]


#Obtaining frequencies and re-organizing dataframe to plot barplot

Bd.char <- snps.lfmm.ann.Bd.unique$Genic_classification
freq_table.bd <- table(Bd.char)
char_df_bd <- as.data.frame(freq_table.bd)
# Rename columns for clarity (optional)
names(char_df_bd) <- c("Region type", "Count")
char_df_bd$Factor <- c("Bd")

Body.char <- snps.lfmm.ann.body.unique$Genic_classification
freq_table.body <- table(Body.char)
char_df_body <- as.data.frame(freq_table.body)
# Rename columns for clarity (optional)
names(char_df_body) <- c("Region type", "Count")
char_df_body$Factor <- c("BC")

temp.char <- snps.lfmm.ann.temp.unique$Genic_classification
freq_table.temp <- table(temp.char)
char_df_temp <- as.data.frame(freq_table.temp)
# Rename columns for clarity (optional)
names(char_df_temp) <- c("Region type", "Count")
char_df_temp$Factor <- c("Temp")


do.char <- snps.lfmm.ann.do.unique$Genic_classification
freq_table.body <- table(do.char)
char_df_do <- as.data.frame(freq_table.body)
# Rename columns for clarity (optional)
names(char_df_do) <- c("Region type", "Count")
char_df_do$Factor <- c("DO")


elev.char <- snps.lfmm.ann.elev.unique$Genic_classification
freq_table.elev <- table(elev.char)
char_df_elev <- as.data.frame(freq_table.elev)
# Rename columns for clarity (optional)
names(char_df_elev) <- c("Region type", "Count")
char_df_elev$Factor <- c("Elevation")

df_factors <- rbind(char_df_bd, char_df_body,char_df_temp,char_df_do,char_df_elev)

## Plotting stacked barplot with percentage

desired_order <- c("Temp", "DO", "Elevation", "Bd","BC")
df_factors$Factor <- factor(df_factors$Factor, levels = desired_order)

barplot1<- ggplot(df_factors, aes(fill=`Region type`, y=Count, x=Factor)) + 
  geom_bar(position="fill", stat="identity")+
  scale_y_continuous(labels = percent)+
  theme(legend.text = element_text(colour = "black", size = 15),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text.x=element_text(face = "bold", colour = "black", size = 15),
        axis.text.y=element_text(colour = "black", size = 12),
        axis.title.y=element_text(face = "bold",colour = "black", size = 15),
        axis.title.x=element_blank())+ scale_fill_manual(values = c("#5a5255","#1b85b8","#BF9675","#FF2C2C"))+
  labs(y = "Percentage  of SNPs (%)",
       x = "Region Type")
barplot1

fig2.1 <- ggarrange(plot.venn2, barplot1, ncol=2, nrow = 1,
                    labels = c("A","B"), font.label = list(size = 15, face = "bold"),
                    common.legend = F)
fig2.1

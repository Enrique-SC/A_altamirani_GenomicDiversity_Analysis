######################################################
###Figure 4 Script - Paper Soto-Cortés et al., 2026###
######################################################


library(fsthet)
library(vcfR)
library(ape)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(ggsci)

#Genetic structure 

#The next command was ran in the Lab cluster, the 
# pca.alt4 <- glPca(Alt.genlight4, scale=T) 
alt.pca.scores4 <- read.csv("alta.pca.scores.meta.csv",header=T)
alt.pca.scores4 <- as.data.frame(pca.alt4$scores)
alt.pca.scores4$Locality <- pop(Alt.genlight4)
alt.pca.scores4$sample <- (Alt.genlight4$ind.names)


# Variance explained calculation  
k <- 100 * pca.alt4$eig/sum(pca.alt4$eig)
var_pc1 <- as.character(round(k[1], digits = 2))
var_pc2 <- as.character(round(k[2], digits = 2)) 
var_pc1 <- paste0("(", var_pc1, "%)")
var_pc2 <- paste0("(", var_pc2, "%)")

pc1 <- "PC1"
pc2 <- "PC2"
pca1_var <- paste(pc1, var_pc1) #Estas son las que van en la gr?fica
pca2_var <- paste(pc2, var_pc2)
desired_order <- c("A", "B", "C", "D","E")

# Convert the 'category' column to a factor with the specified levels
alt.pca.scores4$Locality <- factor(alt.pca.scores4$Locality, levels = desired_order)




plotpca.alt4 <- ggplot(alt.pca.scores4, aes(x=PC1, y=PC2, color=Locality)) +
  geom_point(size=5, pch=21, colour="black", aes(fill=Locality),stroke=1)+ scale_fill_locuszoom() + theme_classic()+
  theme(legend.text = element_text(size = 14,face = "bold"),
        legend.title = element_text(size = 15),
        axis.text.y=element_text(colour = "black", size = 15),
        axis.text.x=element_text(colour = "black", size = 15),
        axis.title.y=element_text(face = "bold",colour = "black", size = 15),
        axis.title.x=element_text(face = "bold",colour = "black", size = 15))
plotpca.alt4
plotpca.alt4 <- plotpca.alt4 + xlab(pca1_var) + ylab(pca2_var)
plotpca.alt4
ggsave("Alt.pca.final.png", plotpca.alt4, width = 10, height = 10,
       limitsize = F, units = "in", path = "./Plots/Finales/")

mapa.sites
plotpca.alt4
admx.5

#Fst heatmap
#alt_fst=genet.dist(Alt.genind4,method="WC84")
#alt_fst.df<- as.matrix(alt_fst)
#alt_fst.df<- as.data.frame(alt_fst.df)
#write.csv(alt_fst.df, "FST.alta.csv", row.names = T)

alt_fst%>%round(digits=3)
fst.matrix= as.matrix(alt_fst) 
ind=which( upper.tri(fst.matrix),arr.ind=TRUE) 
fst.df=data.frame(Site1=dimnames(fst.matrix)[[2]][ind[,2]], 
                  Site2= dimnames(fst.matrix)[[1]][ind[,1]], 
                  Fst=fst.matrix[ind] %>%round(digits=3))


# Keep the order of the levels in the data.frame for plotting 
fst.df$Site1 = factor(fst.df$Site1, levels = unique(fst.df$Site1))
fst.df$Site2 = factor(fst.df$Site2, levels = unique(fst.df$Site2))

#Convertminusvaluestozero 
fst.df$Fst[fst.df$Fst< 0]=0

fst.label=expression(italic("F")[ST])

mid=max(fst.df$Fst)/2

fst.heatmap<- ggplot(data = fst.df, aes(x = Site1, y = Site2, fill = Fst))+
  geom_tile(colour = "black")+
  geom_text(aes(label = Fst), color="black", size = 6)+
  scale_fill_gradient2(low = "#F5F5DC", mid = "#FFD580", high = "#B87333", midpoint = mid, name = fst.label, limits = c(0, max(fst.df$Fst)), breaks = c(0, 0.05, 0.10, 0.15))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), position = "right")+
  theme(axis.text = element_text(colour = "black", size = 16, face = "bold"),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14,face = "bold"))

fst.heatmap

#Admixture plots
#General
samplelist <- read.table("pop.alt2.admix.txt")
all_data <- tibble(sample=character(),
                   k=numeric(),
                   Q=character(),
                   value=numeric())

for (k in 1:5){
  data <- read_delim(paste0("alta.final.chroms.neutros.int.",k,".Q"),
                     col_names = paste0("Q",seq(1:k)),
                     delim=" ")
  data$sample <- samplelist$V3
  data$k <- k
  
  #This step converts from wide to long.
  data %>% gather(Q, value, -sample,-k) -> data
  all_data <- rbind(all_data,data)
}
all_data


admx.3 <- all_data %>%
  filter(k == 3) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(color="black",stat="identity",position="stack", width = 1,size=0.7) +
  xlab("Individuals") + ylab("Ancestry proportion") +
  theme_classic2() +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(colour = "black", face = "bold", size = 12),
        axis.title.y=element_text(face = "bold",colour = "black", size = 15),
        axis.title.x=element_blank(),
        legend.position="none") +
  scale_fill_manual(values = c("#357EBDFF","violet","#EEA236FF"))
admx.3
ggsave("admixture.3.png", admx.3,  width = 12, height = 5,
       limitsize = F, units = "in", path = "./Plots/")

admx.4 <- all_data %>%
  filter(k == 4) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(color="black",stat="identity",position="stack", width = 0.9,size=0.9) +
  xlab("Individuals") + ylab("Ancestry proportion") +
  theme_classic2() +
  theme(axis.text.x=element_text(colour = "black", size = 11,angle = 90, hjust = 1),
        axis.text.y=element_text(colour = "black", face = "bold", size = 12),
        axis.title.y=element_text(face = "bold",colour = "black", size = 15),
        axis.title.x=element_text(face = "bold",colour = "black", size = 12),
        legend.position="none") +
  scale_fill_manual(values = c("#EEA236FF","#357EBDFF","violet","#46B8DAFF"))
admx.4
ggsave("admixture.4.png", admx.4,  width = 12, height = 5,
       limitsize = F, units = "in", path = "./Plots/")

admx.5 <- all_data %>%
  filter(k == 5) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(color="black",stat="identity",position="stack", width = 1,size=0.7) +
  xlab("Individuals") + ylab("Ancestry proportion") +
  theme_classic2() +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(colour = "black", face = "bold", size = 12),
        axis.title.y=element_text(face = "bold",colour = "black", size = 15),
        axis.title.x=element_blank(),
        legend.position="none") +
  scale_fill_manual(values = c("#D43F3AFF","#5CB85CFF","#46B8DAFF","#357EBDFF","#EEA236FF"))
admx.5


fig4 <- ggarrange(ggarrange(plotpca.alt4, fst.heatmap, ncol = 2, labels = c("A", "B","C")), # Second row with box and dot plots
                  nrow = 3,
                  admx.3,admx.5) 
fig4
ggsave("Fig4.pdf", fig4, width = 13, height = 15,
       limitsize = F, units = "in", path = "./Plots/Finales")

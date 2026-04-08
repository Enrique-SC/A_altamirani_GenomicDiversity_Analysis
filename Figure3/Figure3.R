######################################################
###Figure 3 Script - Paper Soto-Cortés et al., 2026###
######################################################

#Loading packages

library(smplot2)

#Loading file with genetic diversity metrics obtained in Figure2 script

gen.alt.out <- read.csv("Div.genet_rmoutliers.csv",header = T)

mlh.smi.plot = ggplot(gen.alt.out, aes(y = MLH, x = SMI)) + 
  geom_point(size = 5) + sm_statCorr(show_text = T,color="#95C182",lwd=4, text_size = 8,corr_method = "spearman")+
  theme(legend.position="none",
        legend.text = element_text(size=13), 
        axis.text.x=element_text(colour = "black", size = 18),
        axis.text.y=element_text(colour = "black", size = 18),
        axis.title.y=element_text(face = "bold",colour = "black", size = 23),
        axis.title.x=element_text(face = "bold",colour = "black", size = 23))+
  xlab("Body condition (Scale Mass Index)")+ylab("Multilocus Heterozygosity (MLH)")
mlh.smi.plot


fis.smi.plot = ggplot(gen.alt.out, aes(y = Fis, x = SMI)) + 
  geom_point(size = 5) + sm_statCorr(show_text = T,color="#D8584E",lwd=4,text_size = 8,corr_method = "spearman")+
  theme(legend.position="none",
        legend.text = element_text(size=13), 
        axis.text.x=element_text(colour = "black", size = 18),
        axis.text.y=element_text(colour = "black", size = 18),
        axis.title.y=element_text(face = "bold",colour = "black", size = 23),
        axis.title.x=element_text(face = "bold",colour = "black", size = 23))+
  xlab("Body condition (Scale Mass Index)")+ylab("Inbreeding coefficient (Fis)")
fis.smi.plot



cor.smi.plot <- ggarrange(mlh.smi.plot, fis.smi.plot, ncol=2, nrow = 1,
                          labels = c("A","B"), font.label = list(size = 13, face = "bold"))

ggsave("Fig3.pdf", cor.smi.plot, width = 16, height = 9,
       limitsize = F, units = "in", path = "./Plots/Finales")

#Spearman correlation tests

cor.test(gen.alt.out$SMI, gen.alt.out$MLH, 
         method = "spearman")

cor.test(gen.alt.out$SMI, gen.alt.out$Fis, 
         method = "spearman")


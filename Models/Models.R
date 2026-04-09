######################################################
###Models script - Paper Soto-Cortés et al., 2026###
######################################################

#Loading packages
library(nlme)
library(lfmm)
library(qvalue)
library(vcfR) 
library(adegenet)

### Linear mixed models

gen.alt.models <- read.csv("Div.genet.alt.org.neutros_rmoutliers.csv",header = T)
gen.alt.models.infected<- gen.alt.models[gen.alt.models$BdLog != 0, ]


lme.mlh.loc <- lme(MLH ~ BdLog+SMI+Morphotype,  random=~1|Locality, data = gen.alt.models.infected)
summary(lme.mlh.loc) 

lme.loc.fis <- lme(Fis ~ BdLog+SMI+Morphotype,  random=~1|Locality, data = gen.alt.models.infected)
summary(lme.loc.fis) 
anova(lme.loc.fis)


### Latent factor mixed model (lfmm)

#Loading SNPS data and converting vcf file to genlight object
Alt.vcf <-read.vcfR("alt.variantsprunnedready.org.noindel.final.chrom.vcf")
Alt.genlight.chr <- vcfR2genlight(Alt.vcf)

#Loading metadata

env.df <- read.csv("metadata_alta_bd_final.csv", header = T)

env.df.lfmm<- env.df[3:4]
env.df.lfmm$DO <- env.df$DO
env.df.lfmm$BdLog <- env.df$BdLog
env.df.lfmm$BodyCondition <- env.df$BodyCondition

lfmm.env <- lfmm_ridge(Y=Alt.genlight.chr, X=env.df.lfmm, K=3) 
lfmm.env.pv <- lfmm_test(Y=Alt.genlight.chr, X=env.df.lfmm, lfmm=lfmm.env, calibrate="gif")
lfmm.env.qv <- qvalue(lfmm.env.pv$calibrated.pvalue)$qvalues
lfmm.env.qv <- as.data.frame(lfmm.env.qv)
write.csv(lfmm.env.qv, "lfmm.env.qv.csv")

#Final SNPs idenified by the model for all variables 
lfmm.env.qv.filtered<- lfmm.env.qv %>% filter(if_any(everything(), ~ . < 0.05))
write.csv(lfmm.env.qv.filtered, "lfmm.env.qv.filtered.csv")


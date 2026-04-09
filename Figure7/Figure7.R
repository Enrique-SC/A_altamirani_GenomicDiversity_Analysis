######################################################
###Figure 7 Script - Paper Soto-Cortés et al., 2026###
######################################################

#Statistical analyses between SNP genotypes
library(FSA)
library(dplyr)
library(rstatix)
library(tidyr)


gt.outliers.lfmm.coding.Bd_t_nozeros <- read.csv("gt.snps.coding.Bd.infected.csv")
gt.outliers.lfmm.coding.body_t <- read.csv("gt.snps.coding.body.csv")


#Bd - Krustal Wallis and pairwise Dun test among genotyes

grouping_vars <- colnames(gt.outliers.lfmm.coding.Bd_t_nozeros)[-c(1, ncol(gt.outliers.lfmm.coding.Bd_t_nozeros))]

results_list <- list()

for (g_var in grouping_vars) {
  # Kruskal-Wallis
  kruskal_res <- gt.outliers.lfmm.coding.Bd_t_nozeros %>%
    kruskal_test(as.formula(paste("BdLog ~", g_var)))
  
  if (kruskal_res$p < 0.05) {
    # Dunn’s post-hoc
    dunn_res <- gt.outliers.lfmm.coding.Bd_t_nozeros %>%
      dunn_test(as.formula(paste("BdLog ~", g_var)), 
                p.adjust.method = "bonferroni") %>%
      select(group1, group2, p, p.adj, p.adj.signif, statistic) %>%
      mutate(SNP = g_var,
             KW_chisq = kruskal_res$statistic,
             KW_df = kruskal_res$df,
             KW_p = kruskal_res$p)
    
    # Mean & SD per group
    summary_stats <- gt.outliers.lfmm.coding.Bd_t_nozeros %>%
      group_by(!!sym(g_var)) %>%
      summarise(
        mean = mean(BdLog, na.rm = TRUE),
        sd   = sd(BdLog, na.rm = TRUE),
        n    = n(),
        .groups = "drop"
      ) %>%
      rename_with(~"Group", !!sym(g_var)) %>%  # rename the SNP column to "Group"
      mutate(SNP = g_var) %>%
      select(SNP, Group, mean, sd, n)
    
    results_list[[g_var]] <- list(
      kruskal = kruskal_res,
      dunn = dunn_res,
      summary = summary_stats
    )
    
  } else {
    results_list[[g_var]] <- list(
      kruskal = kruskal_res,
      dunn = NULL,
      summary = gt.outliers.lfmm.coding.Bd_t_nozeros %>%
        group_by(!!sym(g_var)) %>%
        summarise(
          mean = mean(BdLog, na.rm = TRUE),
          sd   = sd(BdLog, na.rm = TRUE),
          n    = n(),
          .groups = "drop"
        ) %>%
        rename(Group = !!sym(g_var)) %>%
        mutate(SNP = g_var) %>%
        select(SNP, Group, mean, sd, n)
    )
  }
}

# Combine Dunn + Kruskal into one table
dunn_all.coding.Bd <- bind_rows(lapply(results_list, function(x) x$dunn)) %>%
  arrange(SNP)

dunn_all.coding.Bd <- dunn_all.coding.Bd[, c("SNP", names(dunn_all.coding.Bd)[names(dunn_all.coding.Bd) != "SNP"])]

# Combine summaries
summary_all_coding_Bd <- bind_rows(lapply(results_list, function(x) x$summary)) %>%
  arrange(SNP)

# Save outputs
write.csv(dunn_all.coding.Bd, "results_dunn_kruskal_coding_Bd.csv", row.names = FALSE)
write.csv(summary_all_coding_Bd, "results_summary_meansd_coding_Bd.csv", row.names = FALSE)


##Body condition - Krustal Wallis and pairwise Dun test among genotyes

grouping_vars <- colnames(gt.outliers.lfmm.coding.body_t)[-c(1, ncol(gt.outliers.lfmm.coding.body_t))]

results_list <- list()
kruskal_tes
for (g_var in grouping_vars) {
  # Kruskal-Wallis
  kruskal_res <- gt.outliers.lfmm.coding.body_t %>%
    kruskal_test(as.formula(paste("BodyCondition ~", g_var)))
  
  if (kruskal_res$p < 0.05) {
    # Dunn’s post-hoc
    dunn_res <- gt.outliers.lfmm.coding.body_t %>%
      dunn_test(as.formula(paste("BodyCondition ~", g_var)), 
                p.adjust.method = "bonferroni") %>%
      select(group1, group2, p, p.adj, p.adj.signif, statistic) %>%
      mutate(SNP = g_var,
             KW_chisq = kruskal_res$statistic,
             KW_df = kruskal_res$df,
             KW_p = kruskal_res$p)
    
    # Mean & SD per group
    summary_stats <- gt.outliers.lfmm.coding.body_t %>%
      group_by(!!sym(g_var)) %>%
      summarise(
        mean = mean(BodyCondition, na.rm = TRUE),
        sd   = sd(BodyCondition, na.rm = TRUE),
        n    = n(),
        .groups = "drop"
      ) %>%
      rename_with(~"Group", !!sym(g_var)) %>%  # rename the SNP column to "Group"
      mutate(SNP = g_var) %>%
      select(SNP, Group, mean, sd, n)
    
    results_list[[g_var]] <- list(
      kruskal = kruskal_res,
      dunn = dunn_res,
      summary = summary_stats
    )
    
  } else {
    results_list[[g_var]] <- list(
      kruskal = kruskal_res,
      dunn = NULL,
      summary = gt.outliers.lfmm.coding.body_t %>%
        group_by(!!sym(g_var)) %>%
        summarise(
          mean = mean(BodyCondition, na.rm = TRUE),
          sd   = sd(BodyCondition, na.rm = TRUE),
          n    = n(),
          .groups = "drop"
        ) %>%
        rename(Group = !!sym(g_var)) %>%
        mutate(SNP = g_var) %>%
        select(SNP, Group, mean, sd, n)
    )
  }
}

# Combine Dunn + Kruskal into one table
dunn_all.coding.body <- bind_rows(lapply(results_list, function(x) x$dunn)) %>%
  arrange(SNP)

dunn_all.coding.body <- dunn_all.coding.body[, c("SNP", names(dunn_all.coding.body)[names(dunn_all.coding.body) != "SNP"])]

unique(dunn_all.coding.body$SNP)

# Combine summaries
summary_all_coding_body <- bind_rows(lapply(results_list, function(x) x$summary)) %>%
  arrange(SNP)

# Save outputs
write.csv(dunn_all.coding.body, "results_dunn_kruskal_coding_body.csv", row.names = FALSE)
write.csv(summary_all_coding_body, "results_summary_meansd_coding_body.csv", row.names = FALSE)


# Keeping SNPs with Dunn test p < 0.05 #

#Bd SNPs
dunn_all.coding.Bd_0.05 <- dunn_all.coding.Bd[dunn_all.coding.Bd$p.adj.signif < 0.05, ]
snps_to_keep_bd<- unique(dunn_all.coding.Bd_0.05$SNP)
snps_to_keep_bd<- c(snps_to_keep_bd,"BdLog")
gt.snps.coding.bd_10 <- gt.outliers.lfmm.coding.Bd_t_nozeros[, snps_to_keep_bd]

#Body condition snps
dunn_all.coding.body_0.05 <- dunn_all.coding.body[dunn_all.coding.body$p.adj.signif < 0.05, ]
snps_to_keep_body<- unique(dunn_all.coding.body_0.05$SNP)
snps_to_keep_body<- c(snps_to_keep_body,"BodyCondition")
gt.snps.coding.body_17 <- gt.outliers.lfmm.coding.body_t[, snps_to_keep_body]
write.csv(gt.snps.coding.body_17, "gt.snps.coding.body.17.csv",row.names = T)


#Re-shaping dataframe for plotting boxplot

snp.bd.boxplot_long <- gt.snps.coding.bd_10 %>%
  pivot_longer(
    cols = starts_with("CM"), # Selects columns x1, x2, x3
    names_to = "variable",
    values_to = "value"
  )

snp.body.boxplot_long <- gt.snps.coding.body_17 %>%
  pivot_longer(
    cols = starts_with("CM"), # Selects columns x1, x2, x3
    names_to = "variable",
    values_to = "value"
  )

#Plotting boxplots

boxplot.snps.bd <- ggplot(snp.bd.boxplot_long, aes(x =value, y = BdLog)) +
  geom_boxplot(lwd=1,fill="grey",width=0.5) + theme_bw()+
  facet_wrap(~ variable, scales = "free_x") + # Creates separate plots for each x-variable
  labs(y = "Bd infection intensity (BdLoad)", x = "SNP genotype") +
  theme(legend.position="none",
        axis.text.x=element_text(face = "bold",colour = "black", size = 16),
        axis.text.y=element_text(colour = "black", size = 13),
        axis.title.y=element_text(face = "bold",colour = "black", size = 17),
        axis.title.x=element_text(face = "bold",colour = "black", size = 17),
        strip.text = element_text(face = "bold", size = 16, color = "black"))
boxplot.snps.bd
ggsave("Fig7.png", boxplot.snps.bd, width = 13, height = 10,
       limitsize = F, units = "in", path = "./Plots/Finales/")

boxplot.snps.body<- ggplot(snp.body.boxplot_long, aes(x = value, y = BodyCondition)) +
  geom_boxplot(lwd=1,fill="grey") + theme_bw()+
  facet_wrap(~ variable, scales = "free_x") + # Creates separate plots for each x-variable
  labs(y = "Body condition (Scale Mass Index)", x = "SNP genotype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position="none",
        axis.text.x=element_text(face = "bold",colour = "black", size = 16),
        axis.text.y=element_text(colour = "black", size = 13),
        axis.title.y=element_text(face = "bold",colour = "black", size = 17),
        axis.title.x=element_text(face = "bold",colour = "black", size = 17),
        strip.text = element_text(face = "bold", size = 16, color = "black"))
boxplot.snps.body
ggsave("FigS6.png", boxplot.snps.body, width = 15, height = 12,
       limitsize = F, units = "in", path = "./Plots/Finales/")

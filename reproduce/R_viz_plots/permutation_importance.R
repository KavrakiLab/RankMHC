library(data.table)
library(ggplot2)
library(viridis)
library(stringr)

permutation_importance <- fread("./permutation_importance_v2_0.75_50.csv")
permutation_importance[, V1 := NULL]
setDT(permutation_importance)[, pred := str_replace(pred, "pep_", "")]
setDT(permutation_importance)[, pred := str_replace(pred, "rsa_", "")]
setDT(permutation_importance)[, pred := str_replace(pred, "global_", "")]
permutation_importance[, p1 := -p1*100]

# Spearman importance
spearman_ordering <- head(permutation_importance[order(spearman, decreasing = TRUE)], 15)[, .(pred, spearman)]
spearman_ordering[, pred := factor(pred, levels = unlist(spearman_ordering[, pred]), ordered = T)]
ggplot(data = spearman_ordering, aes(x = pred, y = spearman)) +
  geom_bar(position="dodge", stat="identity", width = 0.65, fill = '#287C8EFF') +
  xlab("Feature") +
  ylab("Difference in ρ") +
  theme_bw() +
  theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12), 
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
        strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16),
        axis.title.y = element_text(size = rel(1.6)), axis.title.x = element_text(size = rel(1.6)), 
        axis.text.y = element_text(angle = 0, hjust = 1, size = 16), plot.margin = margin(l = 0 + 23))

# LRMSD importance
lrmsd_ordering <- head(permutation_importance[order(p1, decreasing = TRUE)], 15)[, .(pred, p1)]
lrmsd_ordering[, pred := factor(pred, levels = unlist(lrmsd_ordering[, pred]), ordered = T)]
ggplot(data = lrmsd_ordering, aes(x = pred, y = p1)) +
  geom_bar(position="dodge", stat="identity", width = 0.65, fill = '#287C8EFF') +
  xlab("Feature") +
  ylab("Difference in P@1 (%) ") +
  theme_bw() +
  theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12), 
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
        strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16),
        axis.title.y = element_text(size = rel(1.6)), axis.title.x = element_text(size = rel(1.6)), 
        axis.text.y = element_text(angle = 0, hjust = 1, size = 16), plot.margin = margin(l = 0 + 23))

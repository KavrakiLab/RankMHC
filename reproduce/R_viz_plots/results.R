library(data.table)
library(magrittr)
library(ggplot2)
library(scales)
library(viridis)
library(DescTools)
source("CorrAggBias.R")
library(dplyr)

ndcg <- function(pred, true) {
  true_sorted_by_pred <- true[order(pred, decreasing=TRUE)]
  true_worst <- sort(true, decreasing=TRUE)
  true <- sort(true, decreasing=FALSE)
  log_array <- sapply(seq(1, length(true)), function(x) {1/log2(1 + x)})
  dcg <- sum(true_sorted_by_pred*log_array) - sum(true_worst*log_array)
  idcg <- sum(true*log_array) - sum(true_worst*log_array)
  dcg/idcg
}

reciprocal_rank <- function(true, pred) {
  pred_sorted_by_true <- pred[order(true, decreasing=TRUE)]
  ranks <- seq(1, length(pred_sorted_by_true))
  ranks_sorted_by_pred <- ranks[order(pred_sorted_by_true)]
  1/match(1, ranks_sorted_by_pred)
}

rmsds <- fread("rmsd_scores_original.csv") 
names(rmsds) <- c("pdb_code", "filename", "CA_rmsd", "BB_rmsd", "HA_rmsd", "dscore")
rmsds <- rmsds[, .(pdb_code, filename, HA_rmsd)]
non_redundant <- fread("./Non_rendundant_structs.csv")
rmsds <- merge(rmsds, non_redundant, by = c("pdb_code", "filename"))
precision_denominator <- rmsds[, pdb_code] %>% unique %>% length - 1

rosetta_score <- fread("Rosetta_score_full.csv")[, .(pdb_code, filename, ref2015)]
rmsds <- merge(rmsds, rosetta_score, by = c("pdb_code", "filename"))

gm_score <- fread("gm_modified.csv")
gm_score[, gradock := ifelse(gradock == 0.0, 1000000, gradock)]
rmsds <- merge(rmsds, gm_score, by = c("pdb_code", "filename"))

vv_score <- fread("vv_modified.csv")
rmsds <- merge(rmsds, vv_score, by = c("pdb_code", "filename"))

Anja_score <- fread("3pHLA.csv")
rmsds <- merge(rmsds, Anja_score, by = c("pdb_code", "filename"))

Abella_score <- fread("Abella_results.csv")
rmsds <- merge(rmsds, Abella_score, by = c("pdb_code", "filename"))

peptide_only_SVC_LKPO <- fread("./res_all/SVR_LKPO_Rosetta_peptide_only_HA_rmsd_redundant_res.csv")
peptide_only_SVC_LKPO <- merge(peptide_only_SVC_LKPO, rmsds, by = c("pdb_code", "filename"))
peptide_only_SVC_LKPO <- melt.data.table(peptide_only_SVC_LKPO, measure.vars = c("Score", "HA_rmsd", "ref2015", "gradock", "molPDF", "vina", "vinardo", "3pHLA", "RF_pred"))
peptide_only_SVC_LKPO[variable == "Score", variable := "LinearSVC_Regressor"]
peptide_only_SVC_LKPO[variable == "HA_rmsd", variable := "Best"]
peptide_only_SVC_LKPO[, CV_type := "LKPO"]
peptide_only_regressor_LKPO <- fread("./res_all/LKPO_regr_Rosetta_peptide_only_HA_rmsd_redundant_res.csv")
peptide_only_regressor_LKPO <- merge(peptide_only_regressor_LKPO, rmsds, by = c("pdb_code", "filename"))
peptide_only_regressor_LKPO <- melt.data.table(peptide_only_regressor_LKPO, measure.vars = c("Score", "HA_rmsd", "ref2015", "gradock", "molPDF", "vina", "vinardo", "3pHLA", "RF_pred"))
peptide_only_regressor_LKPO[variable == "Score", variable := "XGB_Regressor"]
peptide_only_regressor_LKPO[variable == "HA_rmsd", variable := "Best"]
peptide_only_regressor_LKPO[, CV_type := "LKPO"]
peptide_only_pairwise_LKPO <- fread("./res_all/LKPO_rank:pairwise_Rosetta_peptide_only_HA_rmsd_redundant_res.csv")
peptide_only_pairwise_LKPO <- merge(peptide_only_pairwise_LKPO, rmsds, by = c("pdb_code", "filename"))
peptide_only_pairwise_LKPO <- melt.data.table(peptide_only_pairwise_LKPO, measure.vars = c("Score", "HA_rmsd", "ref2015", "gradock", "molPDF", "vina", "vinardo", "3pHLA", "RF_pred"))
peptide_only_pairwise_LKPO[variable == "Score", variable := "XGB_Ranker:pairwise"]
peptide_only_pairwise_LKPO[variable == "HA_rmsd", variable := "Best"]
peptide_only_pairwise_LKPO[, CV_type := "LKPO"]
peptide_only_ndcg_LKPO <- fread("./res_all/LKPO_rank:ndcg_Rosetta_peptide_only_HA_rmsd_redundant_res.csv")
peptide_only_ndcg_LKPO <- merge(peptide_only_ndcg_LKPO, rmsds, by = c("pdb_code", "filename"))
peptide_only_ndcg_LKPO <- melt.data.table(peptide_only_ndcg_LKPO, measure.vars = c("Score", "HA_rmsd", "ref2015", "gradock", "molPDF", "vina", "vinardo", "3pHLA", "RF_pred"))
peptide_only_ndcg_LKPO[variable == "Score", variable := "XGB_Ranker:ndcg"]
peptide_only_ndcg_LKPO[variable == "HA_rmsd", variable := "Best"]
peptide_only_ndcg_LKPO[, CV_type := "LKPO"]
peptide_only_pairwise_129_LKPO <- fread("./res_all/LKPO_rank:pairwise_Rosetta_peptide_only_129_HA_rmsd_redundant_129_res.csv")
peptide_only_pairwise_129_LKPO <- merge(peptide_only_pairwise_129_LKPO, rmsds, by = c("pdb_code", "filename"))
peptide_only_pairwise_129_LKPO <- melt.data.table(peptide_only_pairwise_129_LKPO, measure.vars = c("Score", "HA_rmsd", "ref2015", "gradock", "molPDF", "vina", "vinardo", "3pHLA", "RF_pred"))
peptide_only_pairwise_129_LKPO[variable == "Score", variable := "XGB_Ranker:pairwise_129"]
peptide_only_pairwise_129_LKPO[variable == "HA_rmsd", variable := "Best"]
peptide_only_pairwise_129_LKPO[, CV_type := "LKPO"]
peptide_only_pairwise_345678_LKPO <- fread("./res_all/LKPO_rank:pairwise_Rosetta_peptide_only_345678_HA_rmsd_redundant_345678_res.csv")
peptide_only_pairwise_345678_LKPO <- merge(peptide_only_pairwise_345678_LKPO, rmsds, by = c("pdb_code", "filename"))
peptide_only_pairwise_345678_LKPO <- melt.data.table(peptide_only_pairwise_345678_LKPO, measure.vars = c("Score", "HA_rmsd", "ref2015", "gradock", "molPDF", "vina", "vinardo", "3pHLA", "RF_pred"))
peptide_only_pairwise_345678_LKPO[variable == "Score", variable := "XGB_Ranker:pairwise_345678"]
peptide_only_pairwise_345678_LKPO[variable == "HA_rmsd", variable := "Best"]
peptide_only_pairwise_345678_LKPO[, CV_type := "LKPO"]
peptide_only_SVC_LKAO <- fread("./res_all/SVR_LKAO_Rosetta_peptide_only_HA_rmsd_redundant_res.csv")
peptide_only_SVC_LKAO <- merge(peptide_only_SVC_LKAO, rmsds, by = c("pdb_code", "filename"))
peptide_only_SVC_LKAO <- melt.data.table(peptide_only_SVC_LKAO, measure.vars = c("Score", "HA_rmsd", "ref2015", "gradock", "molPDF", "vina", "vinardo", "3pHLA", "RF_pred"))
peptide_only_SVC_LKAO[variable == "Score", variable := "LinearSVC_Regressor"]
peptide_only_SVC_LKAO[variable == "HA_rmsd", variable := "Best"]
peptide_only_SVC_LKAO[, CV_type := "LKAO"]
peptide_only_regressor_LKAO <- fread("./res_all/LKAO_regr_Rosetta_peptide_only_HA_rmsd_redundant_res.csv")
peptide_only_regressor_LKAO <- merge(peptide_only_regressor_LKAO, rmsds, by = c("pdb_code", "filename"))
peptide_only_regressor_LKAO <- melt.data.table(peptide_only_regressor_LKAO, measure.vars = c("Score", "HA_rmsd", "ref2015", "gradock", "molPDF", "vina", "vinardo", "3pHLA", "RF_pred"))
peptide_only_regressor_LKAO[variable == "Score", variable := "XGB_Regressor"]
peptide_only_regressor_LKAO[variable == "HA_rmsd", variable := "Best"]
peptide_only_regressor_LKAO[, CV_type := "LKAO"]
peptide_only_pairwise_LKAO <- fread("./res_all/LKAO_rank:pairwise_Rosetta_peptide_only_HA_rmsd_redundant_res.csv")
peptide_only_pairwise_LKAO <- merge(peptide_only_pairwise_LKAO, rmsds, by = c("pdb_code", "filename"))
peptide_only_pairwise_LKAO <- melt.data.table(peptide_only_pairwise_LKAO, measure.vars = c("Score", "HA_rmsd", "ref2015", "gradock", "molPDF", "vina", "vinardo", "3pHLA", "RF_pred"))
peptide_only_pairwise_LKAO[variable == "Score", variable := "XGB_Ranker:pairwise"]
peptide_only_pairwise_LKAO[variable == "HA_rmsd", variable := "Best"]
peptide_only_pairwise_LKAO[, CV_type := "LKAO"]
peptide_only_ndcg_LKAO <- fread("./res_all/LKAO_rank:ndcg_Rosetta_peptide_only_HA_rmsd_redundant_res.csv")
peptide_only_ndcg_LKAO <- merge(peptide_only_ndcg_LKAO, rmsds, by = c("pdb_code", "filename"))
peptide_only_ndcg_LKAO <- melt.data.table(peptide_only_ndcg_LKAO, measure.vars = c("Score", "HA_rmsd", "ref2015", "gradock", "molPDF", "vina", "vinardo", "3pHLA", "RF_pred"))
peptide_only_ndcg_LKAO[variable == "Score", variable := "XGB_Ranker:ndcg"]
peptide_only_ndcg_LKAO[variable == "HA_rmsd", variable := "Best"]
peptide_only_ndcg_LKAO[, CV_type := "LKAO"]
peptide_only_pairwise_129_LKAO <- fread("./res_all/LKAO_rank:pairwise_Rosetta_peptide_only_129_HA_rmsd_redundant_129_res.csv")
peptide_only_pairwise_129_LKAO <- merge(peptide_only_pairwise_129_LKAO, rmsds, by = c("pdb_code", "filename"))
peptide_only_pairwise_129_LKAO <- melt.data.table(peptide_only_pairwise_129_LKAO, measure.vars = c("Score", "HA_rmsd", "ref2015", "gradock", "molPDF", "vina", "vinardo", "3pHLA", "RF_pred"))
peptide_only_pairwise_129_LKAO[variable == "Score", variable := "XGB_Ranker:pairwise_129"]
peptide_only_pairwise_129_LKAO[variable == "HA_rmsd", variable := "Best"]
peptide_only_pairwise_129_LKAO[, CV_type := "LKAO"]
peptide_only_pairwise_345678_LKAO <- fread("./res_all/LKAO_rank:pairwise_Rosetta_peptide_only_345678_HA_rmsd_redundant_345678_res.csv")
peptide_only_pairwise_345678_LKAO <- merge(peptide_only_pairwise_345678_LKAO, rmsds, by = c("pdb_code", "filename"))
peptide_only_pairwise_345678_LKAO <- melt.data.table(peptide_only_pairwise_345678_LKAO, measure.vars = c("Score", "HA_rmsd", "ref2015", "gradock", "molPDF", "vina", "vinardo", "3pHLA", "RF_pred"))
peptide_only_pairwise_345678_LKAO[variable == "Score", variable := "XGB_Ranker:pairwise_345678"]
peptide_only_pairwise_345678_LKAO[variable == "HA_rmsd", variable := "Best"]
peptide_only_pairwise_345678_LKAO[, CV_type := "LKAO"]
peptide_only_SVC_LOLO <- fread("./res_all/SVR_LOLO_Rosetta_peptide_only_HA_rmsd_redundant_res.csv")
peptide_only_SVC_LOLO <- merge(peptide_only_SVC_LOLO, rmsds, by = c("pdb_code", "filename"))
peptide_only_SVC_LOLO <- melt.data.table(peptide_only_SVC_LOLO, measure.vars = c("Score", "HA_rmsd", "ref2015", "gradock", "molPDF", "vina", "vinardo", "3pHLA", "RF_pred"))
peptide_only_SVC_LOLO[variable == "Score", variable := "LinearSVC_Regressor"]
peptide_only_SVC_LOLO[variable == "HA_rmsd", variable := "Best"]
peptide_only_SVC_LOLO[, CV_type := "LOLO"]
peptide_only_regressor_LOLO <- fread("./res_all/LOLO_regr_Rosetta_peptide_only_HA_rmsd_redundant_res.csv")
peptide_only_regressor_LOLO <- merge(peptide_only_regressor_LOLO, rmsds, by = c("pdb_code", "filename"))
peptide_only_regressor_LOLO <- melt.data.table(peptide_only_regressor_LOLO, measure.vars = c("Score", "HA_rmsd", "ref2015", "gradock", "molPDF", "vina", "vinardo", "3pHLA", "RF_pred"))
peptide_only_regressor_LOLO[variable == "Score", variable := "XGB_Regressor"]
peptide_only_regressor_LOLO[variable == "HA_rmsd", variable := "Best"]
peptide_only_regressor_LOLO[, CV_type := "LOLO"]
peptide_only_pairwise_LOLO <- fread("./res_all/LOLO_rank:pairwise_Rosetta_peptide_only_HA_rmsd_redundant_res.csv")
peptide_only_pairwise_LOLO <- merge(peptide_only_pairwise_LOLO, rmsds, by = c("pdb_code", "filename"))
peptide_only_pairwise_LOLO <- melt.data.table(peptide_only_pairwise_LOLO, measure.vars = c("Score", "HA_rmsd", "ref2015", "gradock", "molPDF", "vina", "vinardo", "3pHLA", "RF_pred"))
peptide_only_pairwise_LOLO[variable == "Score", variable := "XGB_Ranker:pairwise"]
peptide_only_pairwise_LOLO[variable == "HA_rmsd", variable := "Best"]
peptide_only_pairwise_LOLO[, CV_type := "LOLO"]
peptide_only_ndcg_LOLO <- fread("./res_all/LOLO_rank:ndcg_Rosetta_peptide_only_HA_rmsd_redundant_res.csv")
peptide_only_ndcg_LOLO <- merge(peptide_only_ndcg_LOLO, rmsds, by = c("pdb_code", "filename"))
peptide_only_ndcg_LOLO <- melt.data.table(peptide_only_ndcg_LOLO, measure.vars = c("Score", "HA_rmsd", "ref2015", "gradock", "molPDF", "vina", "vinardo", "3pHLA", "RF_pred"))
peptide_only_ndcg_LOLO[variable == "Score", variable := "XGB_Ranker:ndcg"]
peptide_only_ndcg_LOLO[variable == "HA_rmsd", variable := "Best"]
peptide_only_ndcg_LOLO[, CV_type := "LOLO"]
peptide_only_pairwise_129_LOLO <- fread("./res_all/LOLO_rank:pairwise_Rosetta_peptide_only_129_HA_rmsd_redundant_129_res.csv")
peptide_only_pairwise_129_LOLO <- merge(peptide_only_pairwise_129_LOLO, rmsds, by = c("pdb_code", "filename"))
peptide_only_pairwise_129_LOLO <- melt.data.table(peptide_only_pairwise_129_LOLO, measure.vars = c("Score", "HA_rmsd", "ref2015", "gradock", "molPDF", "vina", "vinardo", "3pHLA", "RF_pred"))
peptide_only_pairwise_129_LOLO[variable == "Score", variable := "XGB_Ranker:pairwise_129"]
peptide_only_pairwise_129_LOLO[variable == "HA_rmsd", variable := "Best"]
peptide_only_pairwise_129_LOLO[, CV_type := "LOLO"]
peptide_only_pairwise_345678_LOLO <- fread("./res_all/LOLO_rank:pairwise_Rosetta_peptide_only_345678_HA_rmsd_redundant_345678_res.csv")
peptide_only_pairwise_345678_LOLO <- merge(peptide_only_pairwise_345678_LOLO, rmsds, by = c("pdb_code", "filename"))
peptide_only_pairwise_345678_LOLO <- melt.data.table(peptide_only_pairwise_345678_LOLO, measure.vars = c("Score", "HA_rmsd", "ref2015", "gradock", "molPDF", "vina", "vinardo", "3pHLA", "RF_pred"))
peptide_only_pairwise_345678_LOLO[variable == "Score", variable := "XGB_Ranker:pairwise_345678"]
peptide_only_pairwise_345678_LOLO[variable == "HA_rmsd", variable := "Best"]
peptide_only_pairwise_345678_LOLO[, CV_type := "LOLO"]
df <- rbindlist(list(peptide_only_SVC_LKPO, peptide_only_regressor_LKPO, peptide_only_pairwise_LKPO, peptide_only_ndcg_LKPO, peptide_only_pairwise_129_LKPO, peptide_only_pairwise_345678_LKPO,
                     peptide_only_SVC_LKAO, peptide_only_regressor_LKAO, peptide_only_pairwise_LKAO, peptide_only_ndcg_LKAO, peptide_only_pairwise_129_LKAO, peptide_only_pairwise_345678_LKAO,
                     peptide_only_SVC_LOLO, peptide_only_regressor_LOLO, peptide_only_pairwise_LOLO, peptide_only_ndcg_LOLO, peptide_only_pairwise_129_LOLO, peptide_only_pairwise_345678_LOLO))
df <- merge(df, rmsds, by = c("pdb_code", "filename"))
df[, qid := NULL]
df[, V1 := NULL]
df <- unique(df)
df[, variable := factor(variable, levels = c("ref2015", "gradock", "molPDF", "vina", "vinardo", "3pHLA", "RF_pred", "LinearSVC_Regressor", "XGB_Regressor", "XGB_Ranker:pairwise", "XGB_Ranker:ndcg", "XGB_Ranker:pairwise_129", "XGB_Ranker:pairwise_345678", "Best"), ordered = T)]
df[, CV_type := factor(CV_type, levels = c("LKPO", "LKAO", "LOLO"), ordered = T)]
df[, value := ifelse(variable == "XGB_Ranker:pairwise", -value, value)]
df[, value := ifelse(variable == "XGB_Ranker:ndcg", -value, value)]
df[, value := ifelse(variable == "XGB_Ranker:pairwise_129", -value, value)]
df[, value := ifelse(variable == "XGB_Ranker:pairwise_345678", -value, value)]
df[, value := ifelse(variable == "3pHLA", -value, value)]
df[, value := ifelse(variable == "RF_pred", -value, value)]
df_temp <- df[CV_type == "LKPO" & pdb_code == "7LG3.pdb" & variable %in% c("ref2015", "XGB_Ranker:pairwise"), .(pdb_code, filename, variable, value, HA_rmsd)]

N <- df[variable == "XGB_Ranker:ndcg" & CV_type == "LOLO", .N, by = pdb_code][, N]
df_corr <- df[, cor(value, HA_rmsd, method = "spearman"), by = .(pdb_code, peptide, allele, variable, CV_type)]
df_corr <- df_corr[variable != "Best"]
df_corr[, metric := "Spearman"]
df_corr_ref2015_LKPO <- df_corr[variable == "ref2015" & CV_type == "LKPO", V1]
df_corr_ref2015_LKPO <- MeanR(df_corr_ref2015_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_ref2015_LKAO <- df_corr[variable == "ref2015" & CV_type == "LKAO", V1]
df_corr_ref2015_LKAO <- MeanR(df_corr_ref2015_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_ref2015_LOLO <- df_corr[variable == "ref2015" & CV_type == "LOLO", V1]
df_corr_ref2015_LOLO <- MeanR(df_corr_ref2015_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_gradock_LKPO <- df_corr[variable == "gradock" & CV_type == "LKPO", V1]
df_corr_gradock_LKPO <- MeanR(df_corr_gradock_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_gradock_LKAO <- df_corr[variable == "gradock" & CV_type == "LKAO", V1]
df_corr_gradock_LKAO <- MeanR(df_corr_gradock_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_gradock_LOLO <- df_corr[variable == "gradock" & CV_type == "LOLO", V1]
df_corr_gradock_LOLO <- MeanR(df_corr_gradock_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_3pHLA_LKPO <- df_corr[variable == "3pHLA" & CV_type == "LKPO", V1]
df_corr_3pHLA_LKPO <- MeanR(df_corr_3pHLA_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_3pHLA_LKAO <- df_corr[variable == "3pHLA" & CV_type == "LKAO", V1]
df_corr_3pHLA_LKAO <- MeanR(df_corr_3pHLA_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_3pHLA_LOLO <- df_corr[variable == "3pHLA" & CV_type == "LOLO", V1]
df_corr_3pHLA_LOLO <- MeanR(df_corr_3pHLA_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_RF_LKPO <- df_corr[variable == "RF_pred" & CV_type == "LKPO", V1]
df_corr_RF_LKPO[is.na(df_corr_RF_LKPO)] <- 0
df_corr_RF_LKPO <- MeanR(df_corr_RF_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_RF_LKAO <- df_corr[variable == "RF_pred" & CV_type == "LKAO", V1]
df_corr_RF_LKAO[is.na(df_corr_RF_LKAO)] <- 0
df_corr_RF_LKAO <- MeanR(df_corr_RF_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_RF_LOLO <- df_corr[variable == "RF_pred" & CV_type == "LOLO", V1]
df_corr_RF_LOLO[is.na(df_corr_RF_LOLO)] <- 0
df_corr_RF_LOLO <- MeanR(df_corr_RF_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_molPDF_LKPO <- df_corr[variable == "molPDF" & CV_type == "LKPO", V1]
df_corr_molPDF_LKPO <- MeanR(df_corr_molPDF_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_molPDF_LKAO <- df_corr[variable == "molPDF" & CV_type == "LKAO", V1]
df_corr_molPDF_LKAO <- MeanR(df_corr_molPDF_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_molPDF_LOLO <- df_corr[variable == "molPDF" & CV_type == "LOLO", V1]
df_corr_molPDF_LOLO <- MeanR(df_corr_molPDF_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_vina_LKPO <- df_corr[variable == "vina" & CV_type == "LKPO", V1]
df_corr_vina_LKPO <- MeanR(df_corr_vina_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_vina_LKAO <- df_corr[variable == "vina" & CV_type == "LKAO", V1]
df_corr_vina_LKAO <- MeanR(df_corr_vina_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_vina_LOLO <- df_corr[variable == "vina" & CV_type == "LOLO", V1]
df_corr_vina_LOLO <- MeanR(df_corr_vina_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_vinardo_LKPO <- df_corr[variable == "vinardo" & CV_type == "LKPO", V1]
df_corr_vinardo_LKPO <- MeanR(df_corr_vinardo_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_vinardo_LKAO <- df_corr[variable == "vinardo" & CV_type == "LKAO", V1]
df_corr_vinardo_LKAO <- MeanR(df_corr_vinardo_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_vinardo_LOLO <- df_corr[variable == "vinardo" & CV_type == "LOLO", V1]
df_corr_vinardo_LOLO <- MeanR(df_corr_vinardo_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_SVCregr_LKPO <- df_corr[variable == "LinearSVC_Regressor" & CV_type == "LKPO", V1]
df_corr_SVCregr_LKPO <- MeanR(df_corr_SVCregr_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_SVCregr_LKAO <- df_corr[variable == "LinearSVC_Regressor" & CV_type == "LKAO", V1]
df_corr_SVCregr_LKAO <- MeanR(df_corr_SVCregr_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_SVCregr_LOLO <- df_corr[variable == "LinearSVC_Regressor" & CV_type == "LOLO", V1]
df_corr_SVCregr_LOLO <- MeanR(df_corr_SVCregr_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_XGBregr_LKPO <- df_corr[variable == "XGB_Regressor" & CV_type == "LKPO", V1]
df_corr_XGBregr_LKPO <- MeanR(df_corr_XGBregr_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_XGBregr_LKAO <- df_corr[variable == "XGB_Regressor" & CV_type == "LKAO", V1]
df_corr_XGBregr_LKAO <- MeanR(df_corr_XGBregr_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_XGBregr_LOLO <- df_corr[variable == "XGB_Regressor" & CV_type == "LOLO", V1]
df_corr_XGBregr_LOLO <- MeanR(df_corr_XGBregr_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_XGBpairwise_LKPO <- df_corr[variable == "XGB_Ranker:pairwise" & CV_type == "LKPO", V1]
df_corr_XGBpairwise_LKPO <- MeanR(df_corr_XGBpairwise_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_XGBpairwise_LKAO <- df_corr[variable == "XGB_Ranker:pairwise" & CV_type == "LKAO", V1]
df_corr_XGBpairwise_LKAO <- MeanR(df_corr_XGBpairwise_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_XGBpairwise_LOLO <- df_corr[variable == "XGB_Ranker:pairwise" & CV_type == "LOLO", V1]
df_corr_XGBpairwise_LOLO <- MeanR(df_corr_XGBpairwise_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_XGBndcg_LKPO <- df_corr[variable == "XGB_Ranker:ndcg" & CV_type == "LKPO", V1]
df_corr_XGBndcg_LKPO <- MeanR(df_corr_XGBndcg_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_XGBndcg_LKAO <- df_corr[variable == "XGB_Ranker:ndcg" & CV_type == "LKAO", V1]
df_corr_XGBndcg_LKAO <- MeanR(df_corr_XGBndcg_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_XGBndcg_LOLO <- df_corr[variable == "XGB_Ranker:ndcg" & CV_type == "LOLO", V1]
df_corr_XGBndcg_LOLO <- MeanR(df_corr_XGBndcg_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_XGBpairwise_129_LKPO <- df_corr[variable == "XGB_Ranker:pairwise_129" & CV_type == "LKPO", V1]
df_corr_XGBpairwise_129_LKPO <- MeanR(df_corr_XGBpairwise_129_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_XGBpairwise_129_LKAO <- df_corr[variable == "XGB_Ranker:pairwise_129" & CV_type == "LKAO", V1]
df_corr_XGBpairwise_129_LKAO <- MeanR(df_corr_XGBpairwise_129_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_XGBpairwise_129_LOLO <- df_corr[variable == "XGB_Ranker:pairwise_129" & CV_type == "LOLO", V1]
df_corr_XGBpairwise_129_LOLO <- MeanR(df_corr_XGBpairwise_129_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_XGBpairwise_345678_LKPO <- df_corr[variable == "XGB_Ranker:pairwise_345678" & CV_type == "LKPO", V1]
df_corr_XGBpairwise_345678_LKPO <- MeanR(df_corr_XGBpairwise_345678_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_XGBpairwise_345678_LKAO <- df_corr[variable == "XGB_Ranker:pairwise_345678" & CV_type == "LKAO", V1]
df_corr_XGBpairwise_345678_LKAO <- MeanR(df_corr_XGBpairwise_345678_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_XGBpairwise_345678_LOLO <- df_corr[variable == "XGB_Ranker:pairwise_345678" & CV_type == "LOLO", V1]
df_corr_XGBpairwise_345678_LOLO <- MeanR(df_corr_XGBpairwise_345678_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_mean <- data.table(variable=rep(c("ref2015", "gradock", "molPDF", "vina", "vinardo", "3pHLA", "RF_pred", "LinearSVC_Regressor", "XGB_Regressor", "XGB_Ranker:pairwise", "XGB_Ranker:ndcg", "XGB_Ranker:pairwise_129", "XGB_Ranker:pairwise_345678"), 3),
                           CV_type=c(rep("LKPO", 13), rep("LKAO", 13), rep("LOLO", 13)),
                           metric=rep("Mean rho", 39),
                           value=c(df_corr_ref2015_LKPO, df_corr_gradock_LKPO, df_corr_molPDF_LKPO, df_corr_vina_LKPO, df_corr_vinardo_LKPO, df_corr_3pHLA_LKPO, df_corr_RF_LKPO, df_corr_SVCregr_LKPO, df_corr_XGBregr_LKPO, df_corr_XGBpairwise_LKPO, df_corr_XGBndcg_LKPO, df_corr_XGBpairwise_129_LKPO, df_corr_XGBpairwise_345678_LKPO,
                                   df_corr_ref2015_LKAO, df_corr_gradock_LKAO, df_corr_molPDF_LKAO, df_corr_vina_LKAO, df_corr_vinardo_LKAO, df_corr_3pHLA_LKAO, df_corr_RF_LKAO, df_corr_SVCregr_LKAO, df_corr_XGBregr_LKAO, df_corr_XGBpairwise_LKAO, df_corr_XGBndcg_LKAO, df_corr_XGBpairwise_129_LKAO, df_corr_XGBpairwise_345678_LKAO,
                                   df_corr_ref2015_LOLO, df_corr_gradock_LOLO, df_corr_molPDF_LOLO, df_corr_vina_LOLO, df_corr_vinardo_LOLO, df_corr_3pHLA_LOLO, df_corr_RF_LOLO, df_corr_SVCregr_LOLO, df_corr_XGBregr_LOLO, df_corr_XGBpairwise_LOLO, df_corr_XGBndcg_LOLO, df_corr_XGBpairwise_129_LOLO, df_corr_XGBpairwise_345678_LOLO))
df_ndcg <- df[, ndcg(value, -HA_rmsd), by = .(pdb_code, peptide, allele, variable, CV_type)]
df_ndcg <- df_ndcg[variable != "Best"]
df_ndcg[, metric := "Mean NDCG"]
df_rr <- df[, reciprocal_rank(value, -HA_rmsd), by = .(pdb_code, peptide, allele, variable, CV_type)]
df_rr <- df_rr[variable != "Best"]
df_rr[, metric := "MRR"]
df_p1 <- df[order(pdb_code, value), head(.SD, 1), by = .(pdb_code, variable, CV_type), .SDcols = names(df)]
df_p1 [, min_rmsd := min(HA_rmsd), by = .(pdb_code, variable, CV_type)]
df_p1[, `:=`(pdb_code = NULL, variable = NULL, CV_type = NULL, HA_rmsd = NULL)]
df_precision <- merge(rmsds, df_p1, by = c("pdb_code", "filename"))
df_precision[, `:=`(value = NULL, min_rmsd = NULL, ref2015.x = NULL, ref2015.y = NULL)]
df_p1 <- df_p1[variable != "Best"]
df_p1[, metric := "Mean LRMSD@1"]
df_p1[, `:=`(filename = NULL, V1 = min_rmsd, value = NULL, min_rmsd = NULL, ref2015 = NULL, gradock = NULL, molPDF = NULL, vina = NULL, vinardo = NULL, `3pHLA` = NULL, `RF_pred` = NULL)]
df_precision <- dcast.data.table(df_precision, ... ~ variable, value.var = c("HA_rmsd"))
df_precision <- df_precision[!is.na(Best)]
df_precision <- df_precision[, lapply(.SD, function(x) {sum(!is.na(x))/566}), by = CV_type, .SDcols = c("ref2015", "gradock", "molPDF", "vina", "vinardo", "3pHLA", "RF_pred", "LinearSVC_Regressor", "XGB_Regressor", "XGB_Ranker:pairwise", "XGB_Ranker:ndcg", "XGB_Ranker:pairwise_129", "XGB_Ranker:pairwise_345678", "Best")]
df_precision[, Best := NULL]
df_precision <- melt.data.table(df_precision, id.vars = c("CV_type"))
df_precision[, metric := "Mean P@1"]
df_p2 <- df[order(pdb_code, value), head(.SD, 2), by = .(pdb_code, variable, CV_type), .SDcols = names(df)]
df_p2 [, min_rmsd := min(HA_rmsd), by = .(pdb_code, variable, CV_type)]
df_p2[, `:=`(pdb_code = NULL, variable = NULL, CV_type = NULL, HA_rmsd = NULL)]
df_p2 <- unique(df_p2, by = c("pdb_code", "peptide", "allele", "variable", "CV_type"))
df_p2 <- df_p2[variable != "Best"]
df_p2[, metric := "Best_RMSD@2"]
df_p2[, `:=`(filename = NULL, V1 = min_rmsd, value = NULL, min_rmsd = NULL, ref2015 = NULL)]
df_p3 <- df[order(pdb_code, value), head(.SD, 3), by = .(pdb_code, variable, CV_type), .SDcols = names(df)]
df_p3 [, min_rmsd := min(HA_rmsd), by = .(pdb_code, variable, CV_type)]
df_p3[, `:=`(pdb_code = NULL, variable = NULL, CV_type = NULL, HA_rmsd = NULL)]
df_p3 <- unique(df_p3, by = c("pdb_code", "peptide", "allele", "variable", "CV_type"))
df_p3 <- df_p3[variable != "Best"]
df_p3[, metric := "Best_RMSD@3"]
df_p3[, `:=`(filename = NULL, V1 = min_rmsd, value = NULL, min_rmsd = NULL, ref2015 = NULL)]
df_metric <- rbindlist(list(df_ndcg, df_p1, df_rr), use.names=TRUE)
df_metric[, metric := factor(metric, levels = c("Mean NDCG", "Mean LRMSD@1", "MRR"), ordered = T)]
df_metric[, value := V1]
df_metric[, V1 := NULL]
mean_data <- df_metric[, mean(`value`), by = .(variable, CV_type, metric)]
mean_data[, value := round(V1, digits = 3)]
mean_data[, V1 := NULL]
mean_data <- rbindlist(list(mean_data, df_corr_mean, df_precision), use.names = TRUE)

#mean_df_HA <- df[variable == "ref2015" & CV_type == "LKPO", mean(HA_rmsd), by = allele]
#sd_df_HA <- df[variable == "ref2015" & CV_type == "LKPO", sd(HA_rmsd), by = allele]
#mean_df_HA <- merge(mean_df_HA, sd_df_HA, by = "allele")
#names(mean_df_HA) <- c("allele", "Mean RMSD", "SD")
#mean_df_HA[, label := "Ground Truth"]
#mean_df_Regr <- df[variable == "XGB_Regressor" & CV_type == "LKPO", mean(value), by = allele]
#sd_df_Regr <- df[variable == "XGB_Regressor" & CV_type == "LKPO", sd(value), by = allele]
#mean_df_Regr <- merge(mean_df_Regr, sd_df_Regr, by = "allele")
#names(mean_df_Regr) <- c("allele", "Mean RMSD", "SD")
#mean_df_Regr[, label := "RankMHC (pointwise)"]
#mean_df <- rbindlist(list(mean_df_HA, mean_df_Regr))
#sd_df <- rbindlist(list(sd_df_HA, sd_df_Regr))
#ggplot(mean_df, aes(x = allele, y = `Mean RMSD`, group = label, color = label)) +
#  geom_line() +
#  geom_point() +
#  geom_errorbar(aes(ymin=`Mean RMSD`- SD, ymax=`Mean RMSD`+ SD), width=.2,
#                position=position_dodge(0.05)) +
#  scale_color_viridis(discrete = TRUE) +
#  theme(
#    legend.position="none",
#    plot.title = element_text(size=11)
#  ) +
#  xlab("Alleles") +
#  ylab("Mean L-RMSD") +
#  theme_bw() +
#  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14), plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
#        strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16),
#        axis.title.y = element_text(size = rel(1.6)), axis.title.x = element_text(size = rel(1.6)), 
#        axis.text.y = element_text(angle = 0, hjust = 1, size = 16), legend.text = element_text(size=12))

peptide_only_ranker_LKPO <- fread("./res_all/LKPO_rank:pairwise_Rosetta_peptide_only_HA_rmsd_res.csv")
peptide_only_ranker_LKPO <- merge(peptide_only_ranker_LKPO, rmsds, by = c("pdb_code", "filename"))
peptide_only_ranker_LKPO[, `:=`(gradock = NULL, vina = NULL, vinardo = NULL, molPDF = NULL, `3pHLA` = NULL, `RF_pred` = NULL)]
peptide_only_ranker_LKPO <- melt.data.table(peptide_only_ranker_LKPO, measure.vars = c("Score", "HA_rmsd", "ref2015"))
peptide_only_ranker_LKPO[, variable := ifelse(variable == "Score", "Ranker_nonredundant", ifelse(variable == "HA_rmsd", "Best", "ref2015"))]
peptide_only_ranker_LKPO[, CV_type := "LKPO"]
peptide_only_ranker_LKAO <- fread("./res_all/LKAO_rank:pairwise_Rosetta_peptide_only_HA_rmsd_res.csv")
peptide_only_ranker_LKAO <- merge(peptide_only_ranker_LKAO, rmsds, by = c("pdb_code", "filename"))
peptide_only_ranker_LKAO[, `:=`(gradock = NULL, vina = NULL, vinardo = NULL, molPDF = NULL, `3pHLA` = NULL, `RF_pred` = NULL)]
peptide_only_ranker_LKAO <- melt.data.table(peptide_only_ranker_LKAO, measure.vars = c("Score", "HA_rmsd", "ref2015"))
peptide_only_ranker_LKAO[, variable := ifelse(variable == "Score", "Ranker_nonredundant", ifelse(variable == "HA_rmsd", "Best", "ref2015"))]
peptide_only_ranker_LKAO[, CV_type := "LKAO"]
peptide_only_ranker_LOLO <- fread("./res_all/LOLO_rank:pairwise_Rosetta_peptide_only_HA_rmsd_res.csv")
peptide_only_ranker_LOLO <- merge(peptide_only_ranker_LOLO, rmsds, by = c("pdb_code", "filename"))
peptide_only_ranker_LOLO[, `:=`(gradock = NULL, vina = NULL, vinardo = NULL, molPDF = NULL, `3pHLA` = NULL, `RF_pred` = NULL)]
peptide_only_ranker_LOLO <- melt.data.table(peptide_only_ranker_LOLO, measure.vars = c("Score", "HA_rmsd", "ref2015"))
peptide_only_ranker_LOLO[, variable := ifelse(variable == "Score", "Ranker_nonredundant", ifelse(variable == "HA_rmsd", "Best", "ref2015"))]
peptide_only_ranker_LOLO[, CV_type := "LOLO"]

peptide_only_ranker_redundant_reduced_LKPO <- fread("./res_all/LKPO_rank:pairwise_Rosetta_peptide_only_reduced_HA_rmsd_redundant_res.csv")
peptide_only_ranker_redundant_reduced_LKPO <- merge(peptide_only_ranker_redundant_reduced_LKPO, rmsds, by = c("pdb_code", "filename"))
peptide_only_ranker_redundant_reduced_LKPO[, `:=`(gradock = NULL, vina = NULL, vinardo = NULL, molPDF = NULL, `3pHLA` = NULL, `RF_pred` = NULL)]
peptide_only_ranker_redundant_reduced_LKPO <- melt.data.table(peptide_only_ranker_redundant_reduced_LKPO, measure.vars = c("Score", "HA_rmsd", "ref2015"))
peptide_only_ranker_redundant_reduced_LKPO[, variable := ifelse(variable == "Score", "Ranker_reduced", ifelse(variable == "HA_rmsd", "Best", "ref2015"))]
peptide_only_ranker_redundant_reduced_LKPO[, CV_type := "LKPO"]
peptide_only_ranker_redundant_reduced_LKAO <- fread("./res_all/LKAO_rank:pairwise_Rosetta_peptide_only_reduced_HA_rmsd_redundant_res.csv")
peptide_only_ranker_redundant_reduced_LKAO <- merge(peptide_only_ranker_redundant_reduced_LKAO, rmsds, by = c("pdb_code", "filename"))
peptide_only_ranker_redundant_reduced_LKAO[, `:=`(gradock = NULL, vina = NULL, vinardo = NULL, molPDF = NULL, `3pHLA` = NULL, `RF_pred` = NULL)]
peptide_only_ranker_redundant_reduced_LKAO <- melt.data.table(peptide_only_ranker_redundant_reduced_LKAO, measure.vars = c("Score", "HA_rmsd", "ref2015"))
peptide_only_ranker_redundant_reduced_LKAO[, variable := ifelse(variable == "Score", "Ranker_reduced", ifelse(variable == "HA_rmsd", "Best", "ref2015"))]
peptide_only_ranker_redundant_reduced_LKAO[, CV_type := "LKAO"]
peptide_only_ranker_redundant_reduced_LOLO <- fread("./res_all/LOLO_rank:pairwise_Rosetta_peptide_only_reduced_HA_rmsd_redundant_res.csv")
peptide_only_ranker_redundant_reduced_LOLO <- merge(peptide_only_ranker_redundant_reduced_LOLO, rmsds, by = c("pdb_code", "filename"))
peptide_only_ranker_redundant_reduced_LOLO[, `:=`(gradock = NULL, vina = NULL, vinardo = NULL, molPDF = NULL, `3pHLA` = NULL, `RF_pred` = NULL)]
peptide_only_ranker_redundant_reduced_LOLO <- melt.data.table(peptide_only_ranker_redundant_reduced_LOLO, measure.vars = c("Score", "HA_rmsd", "ref2015"))
peptide_only_ranker_redundant_reduced_LOLO[, variable := ifelse(variable == "Score", "Ranker_reduced", ifelse(variable == "HA_rmsd", "Best", "ref2015"))]
peptide_only_ranker_redundant_reduced_LOLO[, CV_type := "LOLO"]

peptide_only_ranker_redundant_LKPO <- fread("./res_all/LKPO_rank:pairwise_Rosetta_peptide_only_HA_rmsd_redundant_res.csv")
peptide_only_ranker_redundant_LKPO <- merge(peptide_only_ranker_redundant_LKPO, rmsds, by = c("pdb_code", "filename"))
peptide_only_ranker_redundant_LKPO[, `:=`(gradock = NULL, vina = NULL, vinardo = NULL, molPDF = NULL, `3pHLA` = NULL, `RF_pred` = NULL)]
peptide_only_ranker_redundant_LKPO <- melt.data.table(peptide_only_ranker_redundant_LKPO, measure.vars = c("Score", "HA_rmsd", "ref2015"))
peptide_only_ranker_redundant_LKPO[, variable := ifelse(variable == "Score", "Ranker", ifelse(variable == "HA_rmsd", "Best", "ref2015"))]
peptide_only_ranker_redundant_LKPO[, CV_type := "LKPO"]
peptide_only_ranker_redundant_LKAO <- fread("./res_all/LKAO_rank:pairwise_Rosetta_peptide_only_HA_rmsd_redundant_res.csv")
peptide_only_ranker_redundant_LKAO <- merge(peptide_only_ranker_redundant_LKAO, rmsds, by = c("pdb_code", "filename"))
peptide_only_ranker_redundant_LKAO[, `:=`(gradock = NULL, vina = NULL, vinardo = NULL, molPDF = NULL, `3pHLA` = NULL, `RF_pred` = NULL)]
peptide_only_ranker_redundant_LKAO <- melt.data.table(peptide_only_ranker_redundant_LKAO, measure.vars = c("Score", "HA_rmsd", "ref2015"))
peptide_only_ranker_redundant_LKAO[, variable := ifelse(variable == "Score", "Ranker", ifelse(variable == "HA_rmsd", "Best", "ref2015"))]
peptide_only_ranker_redundant_LKAO[, CV_type := "LKAO"]
peptide_only_ranker_redundant_LOLO <- fread("./res_all/LOLO_rank:pairwise_Rosetta_peptide_only_HA_rmsd_redundant_res.csv")
peptide_only_ranker_redundant_LOLO <- merge(peptide_only_ranker_redundant_LOLO, rmsds, by = c("pdb_code", "filename"))
peptide_only_ranker_redundant_LOLO[, `:=`(gradock = NULL, vina = NULL, vinardo = NULL, molPDF = NULL, `3pHLA` = NULL, `RF_pred` = NULL)]
peptide_only_ranker_redundant_LOLO <- melt.data.table(peptide_only_ranker_redundant_LOLO, measure.vars = c("Score", "HA_rmsd", "ref2015"))
peptide_only_ranker_redundant_LOLO[, variable := ifelse(variable == "Score", "Ranker", ifelse(variable == "HA_rmsd", "Best", "ref2015"))]
peptide_only_ranker_redundant_LOLO[, CV_type := "LOLO"]

peptide_only_ranker_redundant_pad_LKPO <- fread("./res_all/LKPO_rank:pairwise_Rosetta_peptide_only_pad_HA_rmsd_redundant_res.csv")
peptide_only_ranker_redundant_pad_LKPO <- merge(peptide_only_ranker_redundant_pad_LKPO, rmsds, by = c("pdb_code", "filename"))
peptide_only_ranker_redundant_pad_LKPO[, `:=`(gradock = NULL, vina = NULL, vinardo = NULL, molPDF = NULL, `3pHLA` = NULL, `RF_pred` = NULL)]
peptide_only_ranker_redundant_pad_LKPO <- melt.data.table(peptide_only_ranker_redundant_pad_LKPO, measure.vars = c("Score", "HA_rmsd", "ref2015"))
peptide_only_ranker_redundant_pad_LKPO[, variable := ifelse(variable == "Score", "Ranker_pad", ifelse(variable == "HA_rmsd", "Best", "ref2015"))]
peptide_only_ranker_redundant_pad_LKPO[, CV_type := "LKPO"]
peptide_only_ranker_redundant_pad_LKAO <- fread("./res_all/LKAO_rank:pairwise_Rosetta_peptide_only_pad_HA_rmsd_redundant_res.csv")
peptide_only_ranker_redundant_pad_LKAO <- merge(peptide_only_ranker_redundant_pad_LKAO, rmsds, by = c("pdb_code", "filename"))
peptide_only_ranker_redundant_pad_LKAO[, `:=`(gradock = NULL, vina = NULL, vinardo = NULL, molPDF = NULL, `3pHLA` = NULL, `RF_pred` = NULL)]
peptide_only_ranker_redundant_pad_LKAO <- melt.data.table(peptide_only_ranker_redundant_pad_LKAO, measure.vars = c("Score", "HA_rmsd", "ref2015"))
peptide_only_ranker_redundant_pad_LKAO[, variable := ifelse(variable == "Score", "Ranker_pad", ifelse(variable == "HA_rmsd", "Best", "ref2015"))]
peptide_only_ranker_redundant_pad_LKAO[, CV_type := "LKAO"]
peptide_only_ranker_redundant_pad_LOLO <- fread("./res_all/LOLO_rank:pairwise_Rosetta_peptide_only_pad_HA_rmsd_redundant_res.csv")
peptide_only_ranker_redundant_pad_LOLO <- merge(peptide_only_ranker_redundant_pad_LOLO, rmsds, by = c("pdb_code", "filename"))
peptide_only_ranker_redundant_pad_LOLO[, `:=`(gradock = NULL, vina = NULL, vinardo = NULL, molPDF = NULL, `3pHLA` = NULL, `RF_pred` = NULL)]
peptide_only_ranker_redundant_pad_LOLO <- melt.data.table(peptide_only_ranker_redundant_pad_LOLO, measure.vars = c("Score", "HA_rmsd", "ref2015"))
peptide_only_ranker_redundant_pad_LOLO[, variable := ifelse(variable == "Score", "Ranker_pad", ifelse(variable == "HA_rmsd", "Best", "ref2015"))]
peptide_only_ranker_redundant_pad_LOLO[, CV_type := "LOLO"]

peptide_HLA_ranker_redundant_LKPO <- fread("./res_all/LKPO_rank:pairwise_Rosetta_peptide_HLA_HA_rmsd_redundant_res.csv")
peptide_HLA_ranker_redundant_LKPO <- merge(peptide_HLA_ranker_redundant_LKPO, rmsds, by = c("pdb_code", "filename"))
peptide_HLA_ranker_redundant_LKPO[, `:=`(gradock = NULL, vina = NULL, vinardo = NULL, molPDF = NULL, `3pHLA` = NULL, `RF_pred` = NULL)]
peptide_HLA_ranker_redundant_LKPO <- melt.data.table(peptide_HLA_ranker_redundant_LKPO, measure.vars = c("Score", "HA_rmsd", "ref2015"))
peptide_HLA_ranker_redundant_LKPO[, variable := ifelse(variable == "Score", "Ranker_HLA", ifelse(variable == "HA_rmsd", "Best", "ref2015"))]
peptide_HLA_ranker_redundant_LKPO[, CV_type := "LKPO"]
peptide_HLA_ranker_redundant_LKAO <- fread("./res_all/LKAO_rank:pairwise_Rosetta_peptide_HLA_HA_rmsd_redundant_res.csv")
peptide_HLA_ranker_redundant_LKAO <- merge(peptide_HLA_ranker_redundant_LKAO, rmsds, by = c("pdb_code", "filename"))
peptide_HLA_ranker_redundant_LKAO[, `:=`(gradock = NULL, vina = NULL, vinardo = NULL, molPDF = NULL, `3pHLA` = NULL, `RF_pred` = NULL)]
peptide_HLA_ranker_redundant_LKAO <- melt.data.table(peptide_HLA_ranker_redundant_LKAO, measure.vars = c("Score", "HA_rmsd", "ref2015"))
peptide_HLA_ranker_redundant_LKAO[, variable := ifelse(variable == "Score", "Ranker_HLA", ifelse(variable == "HA_rmsd", "Best", "ref2015"))]
peptide_HLA_ranker_redundant_LKAO[, CV_type := "LKAO"]
peptide_HLA_ranker_redundant_LOLO <- fread("./res_all/LOLO_rank:pairwise_Rosetta_peptide_HLA_HA_rmsd_redundant_res.csv")
peptide_HLA_ranker_redundant_LOLO <- merge(peptide_HLA_ranker_redundant_LOLO, rmsds, by = c("pdb_code", "filename"))
peptide_HLA_ranker_redundant_LOLO[, `:=`(gradock = NULL, vina = NULL, vinardo = NULL, molPDF = NULL, `3pHLA` = NULL, `RF_pred` = NULL)]
peptide_HLA_ranker_redundant_LOLO <- melt.data.table(peptide_HLA_ranker_redundant_LOLO, measure.vars = c("Score", "HA_rmsd", "ref2015"))
peptide_HLA_ranker_redundant_LOLO[, variable := ifelse(variable == "Score", "Ranker_HLA", ifelse(variable == "HA_rmsd", "Best", "ref2015"))]
peptide_HLA_ranker_redundant_LOLO[, CV_type := "LOLO"]

pairwise_ranker_redundant_LKPO <- fread("./res_all/LKPO_rank:pairwise_Rosetta_pairwise_HA_rmsd_redundant_res.csv")
pairwise_ranker_redundant_LKPO <- merge(pairwise_ranker_redundant_LKPO, rmsds, by = c("pdb_code", "filename"))
pairwise_ranker_redundant_LKPO[, `:=`(gradock = NULL, vina = NULL, vinardo = NULL, molPDF = NULL, `3pHLA` = NULL, `RF_pred` = NULL)]
pairwise_ranker_redundant_LKPO <- melt.data.table(pairwise_ranker_redundant_LKPO, measure.vars = c("Score", "HA_rmsd", "ref2015"))
pairwise_ranker_redundant_LKPO[, variable := ifelse(variable == "Score", "Ranker_pairwise", ifelse(variable == "HA_rmsd", "Best", "ref2015"))]
pairwise_ranker_redundant_LKPO[, CV_type := "LKPO"]
pairwise_ranker_redundant_LKAO <- fread("./res_all/LKAO_rank:pairwise_Rosetta_pairwise_HA_rmsd_redundant_res.csv")
pairwise_ranker_redundant_LKAO <- merge(pairwise_ranker_redundant_LKAO, rmsds, by = c("pdb_code", "filename"))
pairwise_ranker_redundant_LKAO[, `:=`(gradock = NULL, vina = NULL, vinardo = NULL, molPDF = NULL, `3pHLA` = NULL, `RF_pred` = NULL)]
pairwise_ranker_redundant_LKAO <- melt.data.table(pairwise_ranker_redundant_LKAO, measure.vars = c("Score", "HA_rmsd", "ref2015"))
pairwise_ranker_redundant_LKAO[, variable := ifelse(variable == "Score", "Ranker_pairwise", ifelse(variable == "HA_rmsd", "Best", "ref2015"))]
pairwise_ranker_redundant_LKAO[, CV_type := "LKAO"]
pairwise_ranker_redundant_LOLO <- fread("./res_all/LOLO_rank:pairwise_Rosetta_pairwise_HA_rmsd_redundant_res.csv")
pairwise_ranker_redundant_LOLO <- merge(pairwise_ranker_redundant_LOLO, rmsds, by = c("pdb_code", "filename"))
pairwise_ranker_redundant_LOLO[, `:=`(gradock = NULL, vina = NULL, vinardo = NULL, molPDF = NULL, `3pHLA` = NULL, `RF_pred` = NULL)]
pairwise_ranker_redundant_LOLO <- melt.data.table(pairwise_ranker_redundant_LOLO, measure.vars = c("Score", "HA_rmsd", "ref2015"))
pairwise_ranker_redundant_LOLO[, variable := ifelse(variable == "Score", "Ranker_pairwise", ifelse(variable == "HA_rmsd", "Best", "ref2015"))]
pairwise_ranker_redundant_LOLO[, CV_type := "LOLO"]

df <- rbindlist(list(peptide_only_ranker_LKPO, peptide_only_ranker_LKAO, peptide_only_ranker_LOLO,
                     peptide_only_ranker_redundant_reduced_LKPO, peptide_only_ranker_redundant_reduced_LKAO, peptide_only_ranker_redundant_reduced_LOLO,
                     peptide_only_ranker_redundant_LKPO, peptide_only_ranker_redundant_LKAO, peptide_only_ranker_redundant_LOLO,
                     peptide_only_ranker_redundant_pad_LKPO, peptide_only_ranker_redundant_pad_LKAO, peptide_only_ranker_redundant_pad_LOLO,
                     peptide_HLA_ranker_redundant_LKPO, peptide_HLA_ranker_redundant_LKAO, peptide_HLA_ranker_redundant_LOLO,
                     pairwise_ranker_redundant_LKPO, pairwise_ranker_redundant_LKAO, pairwise_ranker_redundant_LOLO))
df <- merge(df, rmsds, by = c("pdb_code", "filename"))
df <- unique(df)
df[, variable := factor(variable, levels = c("ref2015", "Ranker", "Ranker_reduced", "Ranker_nonredundant", "Ranker_pad", "Ranker_HLA", "Ranker_pairwise", "Best"), ordered = T)]
df[, CV_type := factor(CV_type, levels = c("LKPO", "LKAO", "LOLO"), ordered = T)]
df[, value := ifelse(variable == "Ranker", -value, value)]
df[, value := ifelse(variable == "Ranker_reduced", -value, value)]
df[, value := ifelse(variable == "Ranker_nonredundant", -value, value)]
df[, value := ifelse(variable == "Ranker_pad", -value, value)]
df[, value := ifelse(variable == "Ranker_HLA", -value, value)]
df[, value := ifelse(variable == "Ranker_pairwise", -value, value)]

N <- df[variable == "Ranker" & CV_type == "LKPO", .N, by = pdb_code][, N]
df_corr <- df[, cor(value, HA_rmsd, method = "spearman"), by = .(pdb_code, qid, peptide, allele, variable, CV_type)]
df_corr <- df_corr[variable != "Best"]
df_corr[, metric := "Spearman"]
df_corr_ref2015_LKPO <- df_corr[variable == "ref2015" & CV_type == "LKPO", V1]
df_corr_ref2015_LKPO <- MeanR(df_corr_ref2015_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_ref2015_LKAO <- df_corr[variable == "ref2015" & CV_type == "LKAO", V1]
df_corr_ref2015_LKAO <- MeanR(df_corr_ref2015_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_ref2015_LOLO <- df_corr[variable == "ref2015" & CV_type == "LOLO", V1]
df_corr_ref2015_LOLO <- MeanR(df_corr_ref2015_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_Ranker_LKPO <- df_corr[variable == "Ranker" & CV_type == "LKPO", V1]
df_corr_Ranker_LKPO <- MeanR(df_corr_Ranker_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_Ranker_LKAO <- df_corr[variable == "Ranker" & CV_type == "LKAO", V1]
df_corr_Ranker_LKAO <- MeanR(df_corr_Ranker_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_Ranker_LOLO <- df_corr[variable == "Ranker" & CV_type == "LOLO", V1]
df_corr_Ranker_LOLO <- MeanR(df_corr_Ranker_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_Reduced_LKPO <- df_corr[variable == "Ranker_reduced" & CV_type == "LKPO", V1]
df_corr_Reduced_LKPO <- MeanR(df_corr_Reduced_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_Reduced_LKAO <- df_corr[variable == "Ranker_reduced" & CV_type == "LKAO", V1]
df_corr_Reduced_LKAO <- MeanR(df_corr_Reduced_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_Reduced_LOLO <- df_corr[variable == "Ranker_reduced" & CV_type == "LOLO", V1]
df_corr_Reduced_LOLO <- MeanR(df_corr_Reduced_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_NonRedundant_LKPO <- df_corr[variable == "Ranker_nonredundant" & CV_type == "LKPO", V1]
df_corr_NonRedundant_LKPO <- MeanR(df_corr_NonRedundant_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_NonRedundant_LKAO <- df_corr[variable == "Ranker_nonredundant" & CV_type == "LKAO", V1]
df_corr_NonRedundant_LKAO <- MeanR(df_corr_NonRedundant_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_NonRedundant_LOLO <- df_corr[variable == "Ranker_nonredundant" & CV_type == "LOLO", V1]
df_corr_NonRedundant_LOLO <- MeanR(df_corr_NonRedundant_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_pad_LKPO <- df_corr[variable == "Ranker_pad" & CV_type == "LKPO", V1]
df_corr_pad_LKPO <- MeanR(df_corr_pad_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_pad_LKAO <- df_corr[variable == "Ranker_pad" & CV_type == "LKAO", V1]
df_corr_pad_LKAO <- MeanR(df_corr_pad_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_pad_LOLO <- df_corr[variable == "Ranker_pad" & CV_type == "LOLO", V1]
df_corr_pad_LOLO <- MeanR(df_corr_pad_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_HLA_LKPO <- df_corr[variable == "Ranker_HLA" & CV_type == "LKPO", V1]
df_corr_HLA_LKPO <- MeanR(df_corr_HLA_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_HLA_LKAO <- df_corr[variable == "Ranker_HLA" & CV_type == "LKAO", V1]
df_corr_HLA_LKAO <- MeanR(df_corr_HLA_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_HLA_LOLO <- df_corr[variable == "Ranker_HLA" & CV_type == "LOLO", V1]
df_corr_HLA_LOLO <- MeanR(df_corr_HLA_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_pairwise_LKPO <- df_corr[variable == "Ranker_pairwise" & CV_type == "LKPO", V1]
df_corr_pairwise_LKPO <- MeanR(df_corr_pairwise_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_pairwise_LKAO <- df_corr[variable == "Ranker_pairwise" & CV_type == "LKAO", V1]
df_corr_pairwise_LKAO <- MeanR(df_corr_pairwise_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_pairwise_LOLO <- df_corr[variable == "Ranker_pairwise" & CV_type == "LOLO", V1]
df_corr_pairwise_LOLO <- MeanR(df_corr_pairwise_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_mean <- data.table(variable=rep(c("ref2015", "Ranker", "Ranker_reduced", "Ranker_nonredundant", "Ranker_pad", "Ranker_HLA", "Ranker_pairwise"), 3),
                           CV_type=c(rep("LKPO", 7), rep("LKAO", 7), rep("LOLO", 7)),
                           metric=rep("ρ", 21),
                           value=c(df_corr_ref2015_LKPO, df_corr_Ranker_LKPO, df_corr_Reduced_LKPO, df_corr_NonRedundant_LKPO, df_corr_pad_LKPO, df_corr_HLA_LKPO, df_corr_pairwise_LKPO,
                                   df_corr_ref2015_LKAO, df_corr_Ranker_LKAO, df_corr_Reduced_LKAO, df_corr_NonRedundant_LKAO, df_corr_pad_LKAO, df_corr_HLA_LKAO, df_corr_pairwise_LKAO,
                                   df_corr_ref2015_LOLO, df_corr_Ranker_LOLO, df_corr_Reduced_LOLO, df_corr_NonRedundant_LOLO, df_corr_pad_LOLO, df_corr_HLA_LOLO, df_corr_pairwise_LOLO))
df_ndcg <- df[, ndcg(value, -HA_rmsd), by = .(pdb_code, qid, peptide, allele, variable, CV_type)]
df_ndcg <- df_ndcg[variable != "Best"]
df_ndcg[, metric := "NDCG"]
df_rr <- df[, reciprocal_rank(value, -HA_rmsd), by = .(pdb_code, qid, peptide, allele, variable, CV_type)]
df_rr <- df_rr[variable != "Best"]
df_rr[, metric := "MRR"]
df_p1 <- df[order(pdb_code, value), head(.SD, 1), by = .(pdb_code, variable, CV_type), .SDcols = names(df)]
df_p1 [, min_rmsd := min(HA_rmsd), by = .(pdb_code, variable, CV_type)]
df_p1[, `:=`(pdb_code = NULL, variable = NULL, CV_type = NULL, HA_rmsd = NULL)]
df_precision <- merge(rmsds, df_p1, by = c("pdb_code", "filename"))
df_precision[, `:=`(V1 = NULL, qid = NULL, value = NULL, min_rmsd = NULL, ref2015.x = NULL, ref2015.y = NULL)]
df_p1 <- df_p1[variable != "Best"]
df_p1[, metric := "LRMSD@1"]
df_p1[, `:=`(filename = NULL, V1 = min_rmsd, value = NULL, min_rmsd = NULL, ref2015 = NULL, gradock = NULL, molPDF = NULL, vina = NULL, vinardo = NULL, `3pHLA` = NULL, `RF_pred` = NULL)]
df_precision <- dcast.data.table(df_precision, ... ~ variable, value.var = c("HA_rmsd"))
df_precision <- df_precision[!is.na(Best)]
df_precision <- df_precision[, lapply(.SD, function(x) {sum(!is.na(x))/566}), by = CV_type, .SDcols = c("ref2015", "Ranker", "Ranker_reduced", "Ranker_nonredundant", "Ranker_pad", "Ranker_HLA", "Ranker_pairwise", "Best")]
df_precision[, Best := NULL]
df_precision <- melt.data.table(df_precision, id.vars = c("CV_type"))
df_precision[, metric := "P@1"]
df_p2 <- df[order(pdb_code, value), head(.SD, 2), by = .(pdb_code, variable, CV_type), .SDcols = names(df)]
df_p2 [, min_rmsd := min(HA_rmsd), by = .(pdb_code, variable, CV_type)]
df_p2[, `:=`(pdb_code = NULL, variable = NULL, CV_type = NULL, HA_rmsd = NULL)]
df_p2 <- unique(df_p2, by = c("pdb_code", "qid", "peptide", "allele", "variable", "CV_type"))
df_p2 <- df_p2[variable != "Best"]
df_p2[, metric := "Best_RMSD@2"]
df_p2[, `:=`(filename = NULL, V1 = min_rmsd, value = NULL, min_rmsd = NULL, ref2015 = NULL, gradock = NULL, molPDF = NULL, vina = NULL, vinardo = NULL, `3pHLA` = NULL, `RF_pred` = NULL)]
df_p3 <- df[order(pdb_code, value), head(.SD, 3), by = .(pdb_code, variable, CV_type), .SDcols = names(df)]
df_p3 [, min_rmsd := min(HA_rmsd), by = .(pdb_code, variable, CV_type)]
df_p3[, `:=`(pdb_code = NULL, variable = NULL, CV_type = NULL, HA_rmsd = NULL)]
df_p3 <- unique(df_p3, by = c("pdb_code", "qid", "peptide", "allele", "variable", "CV_type"))
df_p3 <- df_p3[variable != "Best"]
df_p3[, metric := "Best_RMSD@3"]
df_p3[, `:=`(filename = NULL, V1 = min_rmsd, value = NULL, min_rmsd = NULL, ref2015 = NULL, gradock = NULL, molPDF = NULL, vina = NULL, vinardo = NULL, `3pHLA` = NULL, `RF_pred` = NULL)]
df_metric <- rbindlist(list(df_ndcg, df_p1, df_rr), use.names=TRUE)
df_metric[, metric := factor(metric, levels = c("NDCG", "LRMSD@1", "MRR"), ordered = T)]
df_metric[, value := V1]
df_metric[, V1 := NULL]
mean_data <- df_metric[, mean(`value`), by = .(variable, CV_type, metric)]
mean_data[, value := round(V1, digits = 3)]
mean_data[, V1 := NULL]
mean_data <- rbindlist(list(mean_data, df_corr_mean, df_precision), use.names = TRUE)
mean_data <- mean_data[variable != "ref2015"][order(metric, CV_type, variable)]
diff_data <- mean_data %>%
  group_by(CV_type, metric) %>%
  mutate(DIFF = ifelse(row_number() == 1, NA, value - first(value))) %>% 
  ungroup() %>%
  as.data.table()
diff_data[, DIFF := ifelse(is.na(DIFF), 0, DIFF)]
diff_data[, DIFF := ifelse(metric == "LRMSD@1", -DIFF, DIFF)]
diff_data <- diff_data[variable != "Ranker"]

diff_data[, variable := as.character(variable)]
diff_data[, variable := ifelse(variable == "Ranker_reduced", "RankMHC - additional features", variable)]
diff_data[, variable := ifelse(variable == "Ranker_nonredundant", "RankMHC - redundant structures", variable)]
diff_data[, variable := ifelse(variable == "Ranker_pad", "RankMHC (padding)", variable)]
diff_data[, variable := ifelse(variable == "Ranker_HLA", "RankMHC + MHC features", variable)]
diff_data[, variable := ifelse(variable == "Ranker_pairwise", "RankMHC + pairwise features", variable)]
symmetric_limits <- function (x) 
{
  max <- max(abs(x))
  c(-max, max)
}

diff_data <- diff_data[variable %in% c("RankMHC - additional features", "RankMHC - redundant structures", "RankMHC (padding)", "RankMHC + MHC features", "RankMHC + pairwise features")]
diff_data[, variable := as.factor(variable)]
diff_data[, variable := factor(variable, levels = c("RankMHC - redundant structures", "RankMHC (padding)", "RankMHC - additional features", "RankMHC + MHC features", "RankMHC + pairwise features"), ordered = T)]
diff_data[metric == "P@1", metric := "P@1 (%)"]
diff_data[metric == "P@1 (%)", DIFF := DIFF*100]
diff_data[, metric := factor(metric, levels = c("ρ", "MRR", "P@1 (%)", "NDCG", "LRMSD@1"), ordered = T)]
ggplot(data = diff_data[CV_type == "LKPO"], aes(x = variable, y = DIFF)) +
  geom_bar(position="dodge", stat="identity", fill = '#287C8EFF') +
  # facet_grid(vars(metric), vars(CV_type)) + 
  facet_wrap(vars(metric), nrow = 1, scales="free_y") + 
  #  coord_flip() +
  geom_hline(yintercept = 0, 
             linetype = "dashed", 
             color = "black",
             linewidth = .3) +
  geom_blank(aes(y = 0, ymin = -0.05, ymax = +0.05)) +
  xlab("") +
  ylab("Difference from baseline RankMHC") +
  theme_bw() +
  theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12), 
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
        strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16),
        axis.title.y = element_text(size = rel(1.3)), axis.title.x = element_text(size = rel(1.3)), 
        axis.text.y = element_text(angle = 0, hjust = 1, size = 16), plot.margin = margin(l = 0 + 44))

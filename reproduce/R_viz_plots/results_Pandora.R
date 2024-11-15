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

rmsds <- fread("res_Pandora/Pandora_rmsd_scores.csv") 
names(rmsds) <- c("pdb_code", "filename", "CA_rmsd", "BB_rmsd", "HA_rmsd")
rmsds <- rmsds[, .(pdb_code, filename, HA_rmsd)]
rmsds <- rmsds[!is.na(rmsds[, HA_rmsd])]
rmsds <- rmsds[pdb_code != "2BVQ.pdb"]
rmsds <- unique(rmsds)

rmsds_pdbfixer <- fread("res_Pandora/Pandora_rmsd_scores_pdbfixer.csv") 
names(rmsds_pdbfixer) <- c("pdb_code", "filename", "CA_rmsd", "BB_rmsd", "HA_rmsd")
rmsds_pdbfixer <- rmsds_pdbfixer[, .(pdb_code, filename, HA_rmsd)]
rmsds_pdbfixer <- rmsds_pdbfixer[!is.na(rmsds_pdbfixer[, HA_rmsd])]
rmsds_pdbfixer <- rmsds_pdbfixer[pdb_code != "2BVQ.pdb"]
rmsds_pdbfixer <- unique(rmsds_pdbfixer)

rmsds_openmm <- fread("res_Pandora/Pandora_rmsd_scores_openmm.csv") 
names(rmsds_openmm) <- c("pdb_code", "filename", "CA_rmsd", "BB_rmsd", "HA_rmsd")
rmsds_openmm <- rmsds_openmm[, .(pdb_code, filename, HA_rmsd)]
rmsds_openmm <- rmsds_openmm[!is.na(rmsds_openmm[, HA_rmsd])]
rmsds_openmm <- rmsds_openmm[pdb_code != "2BVQ.pdb"]
rmsds_openmm <- unique(rmsds_openmm)

rosetta_score <- fread("res_Pandora/Pandora_peptide_only.csv")[, .(pdb_code, filename, global_total_score)]
names(rosetta_score) <- c("pdb_code", "filename", "ref2015")
rosetta_score <- unique(rosetta_score)
rmsds <- merge(rmsds, rosetta_score, by = c("pdb_code", "filename"), nomatch=0)

pandora_vvm <- fread("res_Pandora/pandora_vvm.csv")
rmsds <- merge(rmsds, pandora_vvm, by = c("pdb_code", "filename"), nomatch=0)

rosetta_score_pdbfixer <- fread("res_Pandora/Pandorapdbfixer_peptide_only.csv")[, .(pdb_code, filename, global_total_score)]
names(rosetta_score_pdbfixer) <- c("pdb_code", "filename", "ref2015")
rosetta_score_pdbfixer <- unique(rosetta_score_pdbfixer)
rmsds_pdbfixer <- merge(rmsds_pdbfixer, rosetta_score_pdbfixer, by = c("pdb_code", "filename"), nomatch=0)

pandora_vvm_pdbfixer <- fread("res_Pandora/pandora_pdbfixer_vvm.csv")
rmsds_pdbfixer <- merge(rmsds_pdbfixer, pandora_vvm_pdbfixer, by = c("pdb_code", "filename"), nomatch=0)

rosetta_score_openmm <- fread("res_Pandora/Pandoraopenmm_peptide_only.csv")[, .(pdb_code, filename, global_total_score)]
names(rosetta_score_openmm) <- c("pdb_code", "filename", "ref2015")
rosetta_score_openmm <- unique(rosetta_score_openmm)
rmsds_openmm <- merge(rmsds_openmm, rosetta_score_openmm, by = c("pdb_code", "filename"), nomatch=0)

pandora_vvm_openmm <- fread("res_Pandora/pandora_openmm_vvm.csv")
rmsds_openmm <- merge(rmsds_openmm, pandora_vvm_openmm, by = c("pdb_code", "filename"), nomatch=0)

pdb_code_list <- rmsds[, pdb_code] %>% unique
rmsds_pdbfixer <- rmsds_pdbfixer[pdb_code %in% pdb_code_list]
rmsds_openmm <- rmsds_openmm[pdb_code %in% pdb_code_list]

peptide_only_ndcg_LKPO <- fread("res_Pandora/LKPO_Pandora.csv")
peptide_only_ndcg_LKPO[, `V1` := NULL]
peptide_only_ndcg_LKPO <- peptide_only_ndcg_LKPO[pdb_code != "5GSD.pdb"]
peptide_only_ndcg_LKPO <- unique(peptide_only_ndcg_LKPO)
peptide_only_ndcg_LKPO <- peptide_only_ndcg_LKPO[complete.cases(peptide_only_ndcg_LKPO),]
peptide_only_ndcg_LKPO <- merge(rmsds, peptide_only_ndcg_LKPO, by = c("pdb_code", "filename"), nomatch=0)
peptide_only_ndcg_LKPO <- unique(peptide_only_ndcg_LKPO)
peptide_only_ndcg_LKPO <- peptide_only_ndcg_LKPO[complete.cases(peptide_only_ndcg_LKPO),]
peptide_only_ndcg_LKPO <- melt.data.table(peptide_only_ndcg_LKPO, measure.vars = c("Score", "HA_rmsd", "ref2015", "vina", "vinardo", "molPDF"))
peptide_only_ndcg_LKPO[variable == "Score", variable := "RankMHC"]
peptide_only_ndcg_LKPO[variable == "HA_rmsd", variable := "Best"]
peptide_only_ndcg_LKPO[, CV_type := "LKPO"]

peptide_only_ndcg_LKAO <- fread("res_Pandora/LKAO_Pandora.csv")
peptide_only_ndcg_LKAO[, `V1` := NULL]
peptide_only_ndcg_LKAO <- peptide_only_ndcg_LKAO[pdb_code != "5GSD.pdb"]
peptide_only_ndcg_LKAO <- unique(peptide_only_ndcg_LKAO)
peptide_only_ndcg_LKAO <- peptide_only_ndcg_LKAO[complete.cases(peptide_only_ndcg_LKAO),]
peptide_only_ndcg_LKAO <- merge(rmsds, peptide_only_ndcg_LKAO, by = c("pdb_code", "filename"), nomatch=0)
peptide_only_ndcg_LKAO <- unique(peptide_only_ndcg_LKAO)
peptide_only_ndcg_LKAO <- peptide_only_ndcg_LKAO[complete.cases(peptide_only_ndcg_LKAO),]
peptide_only_ndcg_LKAO <- melt.data.table(peptide_only_ndcg_LKAO, measure.vars = c("Score", "HA_rmsd", "ref2015", "vina", "vinardo", "molPDF"))
peptide_only_ndcg_LKAO[variable == "Score", variable := "RankMHC"]
peptide_only_ndcg_LKAO[variable == "HA_rmsd", variable := "Best"]
peptide_only_ndcg_LKAO[, CV_type := "LKAO"]

peptide_only_ndcg_LOLO <- fread("res_Pandora/LOLO_Pandora.csv")
peptide_only_ndcg_LOLO[, `V1` := NULL]
peptide_only_ndcg_LOLO <- peptide_only_ndcg_LOLO[pdb_code != "5GSD.pdb"]
peptide_only_ndcg_LOLO <- unique(peptide_only_ndcg_LOLO)
peptide_only_ndcg_LOLO <- peptide_only_ndcg_LOLO[complete.cases(peptide_only_ndcg_LOLO),]
peptide_only_ndcg_LOLO <- merge(rmsds, peptide_only_ndcg_LOLO, by = c("pdb_code", "filename"), nomatch=0)
peptide_only_ndcg_LOLO <- unique(peptide_only_ndcg_LOLO)
peptide_only_ndcg_LOLO <- peptide_only_ndcg_LOLO[complete.cases(peptide_only_ndcg_LOLO),]
peptide_only_ndcg_LOLO <- melt.data.table(peptide_only_ndcg_LOLO, measure.vars = c("Score", "HA_rmsd", "ref2015", "vina", "vinardo", "molPDF"))
peptide_only_ndcg_LOLO[variable == "Score", variable := "RankMHC"]
peptide_only_ndcg_LOLO[variable == "HA_rmsd", variable := "Best"]
peptide_only_ndcg_LOLO[, CV_type := "LOLO"]

df <- rbindlist(list(peptide_only_ndcg_LKPO, peptide_only_ndcg_LKAO, peptide_only_ndcg_LOLO))
df <- merge(df, rmsds, by = c("pdb_code", "filename"))
df <- unique(df)
df[, variable := factor(variable, levels = c("ref2015", "vina", "vinardo", "molPDF", "RankMHC", "Best"), ordered = T)]
df[, CV_type := factor(CV_type, levels = c("LKPO", "LKAO", "LOLO"), ordered = T)]
df[, value := ifelse(variable == "RankMHC", -value, value)]

N <- df[variable == "RankMHC" & CV_type == "LOLO", .N, by = pdb_code][, N]
df_corr <- df[, cor(value, HA_rmsd, method = "spearman"), by = .(pdb_code, qid, peptide, allele, variable, CV_type)]
df_corr <- df_corr[variable != "Best"]
df_corr[, metric := "Spearman"]
df_corr_ref2015_LKPO <- df_corr[variable == "ref2015" & CV_type == "LKPO", V1]
df_corr_ref2015_LKPO[is.na(df_corr_ref2015_LKPO)] <- 0
df_corr_ref2015_LKPO <- MeanR(df_corr_ref2015_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_ref2015_LKAO <- df_corr[variable == "ref2015" & CV_type == "LKAO", V1]
df_corr_ref2015_LKAO[is.na(df_corr_ref2015_LKAO)] <- 0
df_corr_ref2015_LKAO <- MeanR(df_corr_ref2015_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_ref2015_LOLO <- df_corr[variable == "ref2015" & CV_type == "LOLO", V1]
df_corr_ref2015_LOLO[is.na(df_corr_ref2015_LOLO)] <- 0
df_corr_ref2015_LOLO <- MeanR(df_corr_ref2015_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_vina_LKPO <- df_corr[variable == "vina" & CV_type == "LKPO", V1]
df_corr_vina_LKPO[is.na(df_corr_vina_LKPO)] <- 0
df_corr_vina_LKPO <- MeanR(df_corr_vina_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_vina_LKAO <- df_corr[variable == "vina" & CV_type == "LKAO", V1]
df_corr_vina_LKAO[is.na(df_corr_vina_LKAO)] <- 0
df_corr_vina_LKAO <- MeanR(df_corr_vina_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_vina_LOLO <- df_corr[variable == "vina" & CV_type == "LOLO", V1]
df_corr_vina_LOLO[is.na(df_corr_vina_LOLO)] <- 0
df_corr_vina_LOLO <- MeanR(df_corr_vina_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_vinardo_LKPO <- df_corr[variable == "vinardo" & CV_type == "LKPO", V1]
df_corr_vinardo_LKPO[is.na(df_corr_vinardo_LKPO)] <- 0
df_corr_vinardo_LKPO <- MeanR(df_corr_vinardo_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_vinardo_LKAO <- df_corr[variable == "vinardo" & CV_type == "LKAO", V1]
df_corr_vinardo_LKAO[is.na(df_corr_vinardo_LKAO)] <- 0
df_corr_vinardo_LKAO <- MeanR(df_corr_vinardo_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_vinardo_LOLO <- df_corr[variable == "vinardo" & CV_type == "LOLO", V1]
df_corr_vinardo_LOLO[is.na(df_corr_vinardo_LOLO)] <- 0
df_corr_vinardo_LOLO <- MeanR(df_corr_vinardo_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_molPDF_LKPO <- df_corr[variable == "molPDF" & CV_type == "LKPO", V1]
df_corr_molPDF_LKPO[is.na(df_corr_molPDF_LKPO)] <- 0
df_corr_molPDF_LKPO <- MeanR(df_corr_molPDF_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_molPDF_LKAO <- df_corr[variable == "molPDF" & CV_type == "LKAO", V1]
df_corr_molPDF_LKAO[is.na(df_corr_molPDF_LKAO)] <- 0
df_corr_molPDF_LKAO <- MeanR(df_corr_molPDF_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_molPDF_LOLO <- df_corr[variable == "molPDF" & CV_type == "LOLO", V1]
df_corr_molPDF_LOLO[is.na(df_corr_molPDF_LOLO)] <- 0
df_corr_molPDF_LOLO <- MeanR(df_corr_molPDF_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_XGBndcg_LKPO <- df_corr[variable == "RankMHC" & CV_type == "LKPO", V1]
df_corr_XGBndcg_LKPO[is.na(df_corr_XGBndcg_LKPO)] <- 0
df_corr_XGBndcg_LKPO <- MeanR(df_corr_XGBndcg_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_XGBndcg_LKAO <- df_corr[variable == "RankMHC" & CV_type == "LKAO", V1]
df_corr_XGBndcg_LKAO[is.na(df_corr_XGBndcg_LKAO)] <- 0
df_corr_XGBndcg_LKAO <- MeanR(df_corr_XGBndcg_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_XGBndcg_LOLO <- df_corr[variable == "RankMHC" & CV_type == "LOLO", V1]
df_corr_XGBndcg_LOLO[is.na(df_corr_XGBndcg_LOLO)] <- 0
df_corr_XGBndcg_LOLO <- MeanR(df_corr_XGBndcg_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_mean <- data.table(variable=rep(c("ref2015", "vina", "vinardo", "molPDF", "RankMHC"), 3),
                           CV_type=c(rep("LKPO", 5), rep("LKAO", 5), rep("LOLO", 5)),
                           metric=rep("ρ", 15),
                           value=c(df_corr_ref2015_LKPO, df_corr_vina_LKPO, df_corr_vinardo_LKPO, df_corr_molPDF_LKPO, df_corr_XGBndcg_LKPO,
                                   df_corr_ref2015_LKAO, df_corr_vina_LKAO, df_corr_vinardo_LKAO, df_corr_molPDF_LKAO, df_corr_XGBndcg_LKAO,
                                   df_corr_ref2015_LOLO, df_corr_vina_LOLO, df_corr_vinardo_LOLO, df_corr_molPDF_LOLO, df_corr_XGBndcg_LOLO))
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
df_precision[, `:=`(qid = NULL, value = NULL, min_rmsd = NULL, ref2015.x = NULL, ref2015.y = NULL)]
df_p1 <- df_p1[variable != "Best"]
df_p1[, metric := "LRMSD@1"]
df_p1[, `:=`(filename = NULL, V1 = min_rmsd, value = NULL, min_rmsd = NULL, ref2015 = NULL, vina = NULL, vinardo = NULL, molPDF = NULL)]
df_precision <- dcast.data.table(df_precision, ... ~ variable, value.var = c("HA_rmsd"))
df_precision <- df_precision[!is.na(Best)]
len_p <- df[, pdb_code] %>% unique %>% length
df_precision <- df_precision[, lapply(.SD, function(x) {sum(!is.na(x))/len_p}), by = CV_type, .SDcols = c("ref2015", "vina", "vinardo", "molPDF", "RankMHC", "Best")]
df_precision[, Best := NULL]
df_precision <- melt.data.table(df_precision, id.vars = c("CV_type"))
df_precision[, metric := "P@1"]
df_p2 <- df[order(pdb_code, value), head(.SD, 2), by = .(pdb_code, variable, CV_type), .SDcols = names(df)]
df_p2 [, min_rmsd := min(HA_rmsd), by = .(pdb_code, variable, CV_type)]
df_p2[, `:=`(pdb_code = NULL, variable = NULL, CV_type = NULL, HA_rmsd = NULL)]
df_p2 <- unique(df_p2, by = c("pdb_code", "qid", "peptide", "allele", "variable", "CV_type"))
df_p2 <- df_p2[variable != "Best"]
df_p2[, metric := "Best_RMSD@2"]
df_p2[, `:=`(filename = NULL, V1 = min_rmsd, value = NULL, min_rmsd = NULL, ref2015 = NULL, vina = NULL, vinardo = NULL, molPDF = NULL)]
df_p3 <- df[order(pdb_code, value), head(.SD, 3), by = .(pdb_code, variable, CV_type), .SDcols = names(df)]
df_p3 [, min_rmsd := min(HA_rmsd), by = .(pdb_code, variable, CV_type)]
df_p3[, `:=`(pdb_code = NULL, variable = NULL, CV_type = NULL, HA_rmsd = NULL)]
df_p3 <- unique(df_p3, by = c("pdb_code", "qid", "peptide", "allele", "variable", "CV_type"))
df_p3 <- df_p3[variable != "Best"]
df_p3[, metric := "Best_RMSD@3"]
df_p3[, `:=`(filename = NULL, V1 = min_rmsd, value = NULL, min_rmsd = NULL, ref2015 = NULL, vina = NULL, vinardo = NULL, molPDF = NULL)]
df_metric <- rbindlist(list(df_ndcg, df_p1, df_rr), use.names=TRUE)
df_metric[, metric := factor(metric, levels = c("NDCG", "LRMSD@1", "MRR"), ordered = T)]
df_metric[, value := V1]
df_metric[, V1 := NULL]
mean_data <- df_metric[, mean(`value`), by = .(variable, CV_type, metric)]
mean_data[, value := round(V1, digits = 3)]
mean_data[, V1 := NULL]
mean_data <- rbindlist(list(mean_data, df_corr_mean, df_precision), use.names = TRUE)
mean_data[, type := "PANDORA"]

peptide_only_ndcg_LKPO <- fread("res_Pandora/LKPO_Pandora_pdbfixer.csv")
peptide_only_ndcg_LKPO[, `V1` := NULL]
peptide_only_ndcg_LKPO <- peptide_only_ndcg_LKPO[pdb_code != "5GSD.pdb"]
peptide_only_ndcg_LKPO <- unique(peptide_only_ndcg_LKPO)
peptide_only_ndcg_LKPO <- peptide_only_ndcg_LKPO[complete.cases(peptide_only_ndcg_LKPO),]
peptide_only_ndcg_LKPO <- merge(rmsds_pdbfixer, peptide_only_ndcg_LKPO, by = c("pdb_code", "filename"), nomatch=0)
peptide_only_ndcg_LKPO <- unique(peptide_only_ndcg_LKPO)
peptide_only_ndcg_LKPO <- peptide_only_ndcg_LKPO[complete.cases(peptide_only_ndcg_LKPO),]
peptide_only_ndcg_LKPO <- melt.data.table(peptide_only_ndcg_LKPO, measure.vars = c("Score", "HA_rmsd", "ref2015", "vina", "vinardo", "molPDF"))
peptide_only_ndcg_LKPO[variable == "Score", variable := "RankMHC"]
peptide_only_ndcg_LKPO[variable == "HA_rmsd", variable := "Best"]
peptide_only_ndcg_LKPO[, CV_type := "LKPO"]

peptide_only_ndcg_LKAO <- fread("res_Pandora/LKAO_Pandora_pdbfixer.csv")
peptide_only_ndcg_LKAO[, `V1` := NULL]
peptide_only_ndcg_LKAO <- peptide_only_ndcg_LKAO[pdb_code != "5GSD.pdb"]
peptide_only_ndcg_LKAO <- unique(peptide_only_ndcg_LKAO)
peptide_only_ndcg_LKAO <- peptide_only_ndcg_LKAO[complete.cases(peptide_only_ndcg_LKAO),]
peptide_only_ndcg_LKAO <- merge(rmsds_pdbfixer, peptide_only_ndcg_LKAO, by = c("pdb_code", "filename"), nomatch=0)
peptide_only_ndcg_LKAO <- unique(peptide_only_ndcg_LKAO)
peptide_only_ndcg_LKAO <- peptide_only_ndcg_LKAO[complete.cases(peptide_only_ndcg_LKAO),]
peptide_only_ndcg_LKAO <- melt.data.table(peptide_only_ndcg_LKAO, measure.vars = c("Score", "HA_rmsd", "ref2015", "vina", "vinardo", "molPDF"))
peptide_only_ndcg_LKAO[variable == "Score", variable := "RankMHC"]
peptide_only_ndcg_LKAO[variable == "HA_rmsd", variable := "Best"]
peptide_only_ndcg_LKAO[, CV_type := "LKAO"]

peptide_only_ndcg_LOLO <- fread("res_Pandora/LOLO_Pandora_pdbfixer.csv")
peptide_only_ndcg_LOLO[, `V1` := NULL]
peptide_only_ndcg_LOLO <- peptide_only_ndcg_LOLO[pdb_code != "5GSD.pdb"]
peptide_only_ndcg_LOLO <- unique(peptide_only_ndcg_LOLO)
peptide_only_ndcg_LOLO <- peptide_only_ndcg_LOLO[complete.cases(peptide_only_ndcg_LOLO),]
peptide_only_ndcg_LOLO <- merge(rmsds_pdbfixer, peptide_only_ndcg_LOLO, by = c("pdb_code", "filename"), nomatch=0)
peptide_only_ndcg_LOLO <- unique(peptide_only_ndcg_LOLO)
peptide_only_ndcg_LOLO <- peptide_only_ndcg_LOLO[complete.cases(peptide_only_ndcg_LOLO),]
peptide_only_ndcg_LOLO <- melt.data.table(peptide_only_ndcg_LOLO, measure.vars = c("Score", "HA_rmsd", "ref2015", "vina", "vinardo", "molPDF"))
peptide_only_ndcg_LOLO[variable == "Score", variable := "RankMHC"]
peptide_only_ndcg_LOLO[variable == "HA_rmsd", variable := "Best"]
peptide_only_ndcg_LOLO[, CV_type := "LOLO"]

df <- rbindlist(list(peptide_only_ndcg_LKPO, peptide_only_ndcg_LKAO, peptide_only_ndcg_LOLO))
df <- merge(df, rmsds, by = c("pdb_code", "filename"))
df <- unique(df)
df[, variable := factor(variable, levels = c("ref2015", "vina", "vinardo", "molPDF", "RankMHC", "Best"), ordered = T)]
df[, CV_type := factor(CV_type, levels = c("LKPO", "LKAO", "LOLO"), ordered = T)]
df[, value := ifelse(variable == "RankMHC", -value, value)]

N <- df[variable == "RankMHC" & CV_type == "LOLO", .N, by = pdb_code][, N]
df_corr <- df[, cor(value, HA_rmsd, method = "spearman"), by = .(pdb_code, qid, peptide, allele, variable, CV_type)]
df_corr <- df_corr[variable != "Best"]
df_corr[, metric := "Spearman"]
df_corr_ref2015_LKPO <- df_corr[variable == "ref2015" & CV_type == "LKPO", V1]
df_corr_ref2015_LKPO[is.na(df_corr_ref2015_LKPO)] <- 0
df_corr_ref2015_LKPO <- MeanR(df_corr_ref2015_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_ref2015_LKAO <- df_corr[variable == "ref2015" & CV_type == "LKAO", V1]
df_corr_ref2015_LKAO[is.na(df_corr_ref2015_LKAO)] <- 0
df_corr_ref2015_LKAO <- MeanR(df_corr_ref2015_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_ref2015_LOLO <- df_corr[variable == "ref2015" & CV_type == "LOLO", V1]
df_corr_ref2015_LOLO[is.na(df_corr_ref2015_LOLO)] <- 0
df_corr_ref2015_LOLO <- MeanR(df_corr_ref2015_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_vina_LKPO <- df_corr[variable == "vina" & CV_type == "LKPO", V1]
df_corr_vina_LKPO[is.na(df_corr_vina_LKPO)] <- 0
df_corr_vina_LKPO <- MeanR(df_corr_vina_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_vina_LKAO <- df_corr[variable == "vina" & CV_type == "LKAO", V1]
df_corr_vina_LKAO[is.na(df_corr_vina_LKAO)] <- 0
df_corr_vina_LKAO <- MeanR(df_corr_vina_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_vina_LOLO <- df_corr[variable == "vina" & CV_type == "LOLO", V1]
df_corr_vina_LOLO[is.na(df_corr_vina_LOLO)] <- 0
df_corr_vina_LOLO <- MeanR(df_corr_vina_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_vinardo_LKPO <- df_corr[variable == "vinardo" & CV_type == "LKPO", V1]
df_corr_vinardo_LKPO[is.na(df_corr_vinardo_LKPO)] <- 0
df_corr_vinardo_LKPO <- MeanR(df_corr_vinardo_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_vinardo_LKAO <- df_corr[variable == "vinardo" & CV_type == "LKAO", V1]
df_corr_vinardo_LKAO[is.na(df_corr_vinardo_LKAO)] <- 0
df_corr_vinardo_LKAO <- MeanR(df_corr_vinardo_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_vinardo_LOLO <- df_corr[variable == "vinardo" & CV_type == "LOLO", V1]
df_corr_vinardo_LOLO[is.na(df_corr_vinardo_LOLO)] <- 0
df_corr_vinardo_LOLO <- MeanR(df_corr_vinardo_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_molPDF_LKPO <- df_corr[variable == "molPDF" & CV_type == "LKPO", V1]
df_corr_molPDF_LKPO[is.na(df_corr_molPDF_LKPO)] <- 0
df_corr_molPDF_LKPO <- MeanR(df_corr_molPDF_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_molPDF_LKAO <- df_corr[variable == "molPDF" & CV_type == "LKAO", V1]
df_corr_molPDF_LKAO[is.na(df_corr_molPDF_LKAO)] <- 0
df_corr_molPDF_LKAO <- MeanR(df_corr_molPDF_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_molPDF_LOLO <- df_corr[variable == "molPDF" & CV_type == "LOLO", V1]
df_corr_molPDF_LOLO[is.na(df_corr_molPDF_LOLO)] <- 0
df_corr_molPDF_LOLO <- MeanR(df_corr_molPDF_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_XGBndcg_LKPO <- df_corr[variable == "RankMHC" & CV_type == "LKPO", V1]
df_corr_XGBndcg_LKPO[is.na(df_corr_XGBndcg_LKPO)] <- 0
df_corr_XGBndcg_LKPO <- MeanR(df_corr_XGBndcg_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_XGBndcg_LKAO <- df_corr[variable == "RankMHC" & CV_type == "LKAO", V1]
df_corr_XGBndcg_LKAO[is.na(df_corr_XGBndcg_LKAO)] <- 0
df_corr_XGBndcg_LKAO <- MeanR(df_corr_XGBndcg_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_XGBndcg_LOLO <- df_corr[variable == "RankMHC" & CV_type == "LOLO", V1]
df_corr_XGBndcg_LOLO[is.na(df_corr_XGBndcg_LOLO)] <- 0
df_corr_XGBndcg_LOLO <- MeanR(df_corr_XGBndcg_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_mean <- data.table(variable=rep(c("ref2015", "vina", "vinardo", "molPDF", "RankMHC"), 3),
                           CV_type=c(rep("LKPO", 5), rep("LKAO", 5), rep("LOLO", 5)),
                           metric=rep("ρ", 15),
                           value=c(df_corr_ref2015_LKPO, df_corr_vina_LKPO, df_corr_vinardo_LKPO, df_corr_molPDF_LKPO, df_corr_XGBndcg_LKPO,
                                   df_corr_ref2015_LKAO, df_corr_vina_LKAO, df_corr_vinardo_LKAO, df_corr_molPDF_LKAO, df_corr_XGBndcg_LKAO,
                                   df_corr_ref2015_LOLO, df_corr_vina_LOLO, df_corr_vinardo_LOLO, df_corr_molPDF_LOLO, df_corr_XGBndcg_LOLO))
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
df_precision[, `:=`(qid = NULL, value = NULL, min_rmsd = NULL, ref2015.x = NULL, ref2015.y = NULL)]
df_p1 <- df_p1[variable != "Best"]
df_p1[, metric := "LRMSD@1"]
df_p1[, `:=`(filename = NULL, V1 = min_rmsd, value = NULL, min_rmsd = NULL, ref2015 = NULL, vina = NULL, vinardo = NULL, molPDF = NULL)]
df_precision <- dcast.data.table(df_precision, ... ~ variable, value.var = c("HA_rmsd"))
df_precision <- df_precision[!is.na(Best)]
len_p <- df[, pdb_code] %>% unique %>% length
df_precision <- df_precision[, lapply(.SD, function(x) {sum(!is.na(x))/len_p}), by = CV_type, .SDcols = c("ref2015", "vina", "vinardo", "molPDF", "RankMHC", "Best")]
df_precision[, Best := NULL]
df_precision <- melt.data.table(df_precision, id.vars = c("CV_type"))
df_precision[, metric := "P@1"]
df_p2 <- df[order(pdb_code, value), head(.SD, 2), by = .(pdb_code, variable, CV_type), .SDcols = names(df)]
df_p2 [, min_rmsd := min(HA_rmsd), by = .(pdb_code, variable, CV_type)]
df_p2[, `:=`(pdb_code = NULL, variable = NULL, CV_type = NULL, HA_rmsd = NULL)]
df_p2 <- unique(df_p2, by = c("pdb_code", "qid", "peptide", "allele", "variable", "CV_type"))
df_p2 <- df_p2[variable != "Best"]
df_p2[, metric := "Best_RMSD@2"]
df_p2[, `:=`(filename = NULL, V1 = min_rmsd, value = NULL, min_rmsd = NULL, ref2015 = NULL, vina = NULL, vinardo = NULL, molPDF = NULL)]
df_p3 <- df[order(pdb_code, value), head(.SD, 3), by = .(pdb_code, variable, CV_type), .SDcols = names(df)]
df_p3 [, min_rmsd := min(HA_rmsd), by = .(pdb_code, variable, CV_type)]
df_p3[, `:=`(pdb_code = NULL, variable = NULL, CV_type = NULL, HA_rmsd = NULL)]
df_p3 <- unique(df_p3, by = c("pdb_code", "qid", "peptide", "allele", "variable", "CV_type"))
df_p3 <- df_p3[variable != "Best"]
df_p3[, metric := "Best_RMSD@3"]
df_p3[, `:=`(filename = NULL, V1 = min_rmsd, value = NULL, min_rmsd = NULL, ref2015 = NULL, vina = NULL, vinardo = NULL, molPDF = NULL)]
df_metric <- rbindlist(list(df_ndcg, df_p1, df_rr), use.names=TRUE)
df_metric[, metric := factor(metric, levels = c("NDCG", "LRMSD@1", "MRR"), ordered = T)]
df_metric[, value := V1]
df_metric[, V1 := NULL]
mean_data_pdbfixer <- df_metric[, mean(`value`), by = .(variable, CV_type, metric)]
mean_data_pdbfixer[, value := round(V1, digits = 3)]
mean_data_pdbfixer[, V1 := NULL]
mean_data_pdbfixer <- rbindlist(list(mean_data_pdbfixer, df_corr_mean, df_precision), use.names = TRUE)
mean_data_pdbfixer[, type := "PANDORA + HA"]

peptide_only_ndcg_LKPO <- fread("res_Pandora/LKPO_Pandora_openmm.csv")
peptide_only_ndcg_LKPO[, `V1` := NULL]
peptide_only_ndcg_LKPO <- peptide_only_ndcg_LKPO[pdb_code != "5GSD.pdb"]
peptide_only_ndcg_LKPO <- unique(peptide_only_ndcg_LKPO)
peptide_only_ndcg_LKPO <- peptide_only_ndcg_LKPO[complete.cases(peptide_only_ndcg_LKPO),]
peptide_only_ndcg_LKPO <- merge(rmsds_openmm, peptide_only_ndcg_LKPO, by = c("pdb_code", "filename"), nomatch=0)
peptide_only_ndcg_LKPO <- unique(peptide_only_ndcg_LKPO)
peptide_only_ndcg_LKPO <- peptide_only_ndcg_LKPO[complete.cases(peptide_only_ndcg_LKPO),]
peptide_only_ndcg_LKPO <- melt.data.table(peptide_only_ndcg_LKPO, measure.vars = c("Score", "HA_rmsd", "ref2015", "vina", "vinardo", "molPDF"))
peptide_only_ndcg_LKPO[variable == "Score", variable := "RankMHC"]
peptide_only_ndcg_LKPO[variable == "HA_rmsd", variable := "Best"]
peptide_only_ndcg_LKPO[, CV_type := "LKPO"]

peptide_only_ndcg_LKAO <- fread("res_Pandora/LKAO_Pandora_openmm.csv")
peptide_only_ndcg_LKAO[, `V1` := NULL]
peptide_only_ndcg_LKAO <- peptide_only_ndcg_LKAO[pdb_code != "5GSD.pdb"]
peptide_only_ndcg_LKAO <- unique(peptide_only_ndcg_LKAO)
peptide_only_ndcg_LKAO <- peptide_only_ndcg_LKAO[complete.cases(peptide_only_ndcg_LKAO),]
peptide_only_ndcg_LKAO <- merge(rmsds_openmm, peptide_only_ndcg_LKAO, by = c("pdb_code", "filename"), nomatch=0)
peptide_only_ndcg_LKAO <- unique(peptide_only_ndcg_LKAO)
peptide_only_ndcg_LKAO <- peptide_only_ndcg_LKAO[complete.cases(peptide_only_ndcg_LKAO),]
peptide_only_ndcg_LKAO <- melt.data.table(peptide_only_ndcg_LKAO, measure.vars = c("Score", "HA_rmsd", "ref2015", "vina", "vinardo", "molPDF"))
peptide_only_ndcg_LKAO[variable == "Score", variable := "RankMHC"]
peptide_only_ndcg_LKAO[variable == "HA_rmsd", variable := "Best"]
peptide_only_ndcg_LKAO[, CV_type := "LKAO"]

peptide_only_ndcg_LOLO <- fread("res_Pandora/LOLO_Pandora_openmm.csv")
peptide_only_ndcg_LOLO[, `V1` := NULL]
peptide_only_ndcg_LOLO <- peptide_only_ndcg_LOLO[pdb_code != "5GSD.pdb"]
peptide_only_ndcg_LOLO <- unique(peptide_only_ndcg_LOLO)
peptide_only_ndcg_LOLO <- peptide_only_ndcg_LOLO[complete.cases(peptide_only_ndcg_LOLO),]
peptide_only_ndcg_LOLO <- merge(rmsds_openmm, peptide_only_ndcg_LOLO, by = c("pdb_code", "filename"), nomatch=0)
peptide_only_ndcg_LOLO <- unique(peptide_only_ndcg_LOLO)
peptide_only_ndcg_LOLO <- peptide_only_ndcg_LOLO[complete.cases(peptide_only_ndcg_LOLO),]
peptide_only_ndcg_LOLO <- melt.data.table(peptide_only_ndcg_LOLO, measure.vars = c("Score", "HA_rmsd", "ref2015", "vina", "vinardo", "molPDF"))
peptide_only_ndcg_LOLO[variable == "Score", variable := "RankMHC"]
peptide_only_ndcg_LOLO[variable == "HA_rmsd", variable := "Best"]
peptide_only_ndcg_LOLO[, CV_type := "LOLO"]

df <- rbindlist(list(peptide_only_ndcg_LKPO, peptide_only_ndcg_LKAO, peptide_only_ndcg_LOLO))
df <- merge(df, rmsds, by = c("pdb_code", "filename"))
df <- unique(df)
df[, variable := factor(variable, levels = c("ref2015", "vina", "vinardo", "molPDF", "RankMHC", "Best"), ordered = T)]
df[, CV_type := factor(CV_type, levels = c("LKPO", "LKAO", "LOLO"), ordered = T)]
df[, value := ifelse(variable == "RankMHC", -value, value)]

N <- df[variable == "RankMHC" & CV_type == "LOLO", .N, by = pdb_code][, N]
df_corr <- df[, cor(value, HA_rmsd, method = "spearman"), by = .(pdb_code, qid, peptide, allele, variable, CV_type)]
df_corr <- df_corr[variable != "Best"]
df_corr[, metric := "Spearman"]
df_corr_ref2015_LKPO <- df_corr[variable == "ref2015" & CV_type == "LKPO", V1]
df_corr_ref2015_LKPO[is.na(df_corr_ref2015_LKPO)] <- 0
df_corr_ref2015_LKPO <- MeanR(df_corr_ref2015_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_ref2015_LKAO <- df_corr[variable == "ref2015" & CV_type == "LKAO", V1]
df_corr_ref2015_LKAO[is.na(df_corr_ref2015_LKAO)] <- 0
df_corr_ref2015_LKAO <- MeanR(df_corr_ref2015_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_ref2015_LOLO <- df_corr[variable == "ref2015" & CV_type == "LOLO", V1]
df_corr_ref2015_LOLO[is.na(df_corr_ref2015_LOLO)] <- 0
df_corr_ref2015_LOLO <- MeanR(df_corr_ref2015_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_vina_LKPO <- df_corr[variable == "vina" & CV_type == "LKPO", V1]
df_corr_vina_LKPO[is.na(df_corr_vina_LKPO)] <- 0
df_corr_vina_LKPO <- MeanR(df_corr_vina_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_vina_LKAO <- df_corr[variable == "vina" & CV_type == "LKAO", V1]
df_corr_vina_LKAO[is.na(df_corr_vina_LKAO)] <- 0
df_corr_vina_LKAO <- MeanR(df_corr_vina_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_vina_LOLO <- df_corr[variable == "vina" & CV_type == "LOLO", V1]
df_corr_vina_LOLO[is.na(df_corr_vina_LOLO)] <- 0
df_corr_vina_LOLO <- MeanR(df_corr_vina_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_vinardo_LKPO <- df_corr[variable == "vinardo" & CV_type == "LKPO", V1]
df_corr_vinardo_LKPO[is.na(df_corr_vinardo_LKPO)] <- 0
df_corr_vinardo_LKPO <- MeanR(df_corr_vinardo_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_vinardo_LKAO <- df_corr[variable == "vinardo" & CV_type == "LKAO", V1]
df_corr_vinardo_LKAO[is.na(df_corr_vinardo_LKAO)] <- 0
df_corr_vinardo_LKAO <- MeanR(df_corr_vinardo_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_vinardo_LOLO <- df_corr[variable == "vinardo" & CV_type == "LOLO", V1]
df_corr_vinardo_LOLO[is.na(df_corr_vinardo_LOLO)] <- 0
df_corr_vinardo_LOLO <- MeanR(df_corr_vinardo_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_molPDF_LKPO <- df_corr[variable == "molPDF" & CV_type == "LKPO", V1]
df_corr_molPDF_LKPO[is.na(df_corr_molPDF_LKPO)] <- 0
df_corr_molPDF_LKPO <- MeanR(df_corr_molPDF_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_molPDF_LKAO <- df_corr[variable == "molPDF" & CV_type == "LKAO", V1]
df_corr_molPDF_LKAO[is.na(df_corr_molPDF_LKAO)] <- 0
df_corr_molPDF_LKAO <- MeanR(df_corr_molPDF_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_molPDF_LOLO <- df_corr[variable == "molPDF" & CV_type == "LOLO", V1]
df_corr_molPDF_LOLO[is.na(df_corr_molPDF_LOLO)] <- 0
df_corr_molPDF_LOLO <- MeanR(df_corr_molPDF_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_XGBndcg_LKPO <- df_corr[variable == "RankMHC" & CV_type == "LKPO", V1]
df_corr_XGBndcg_LKPO[is.na(df_corr_XGBndcg_LKPO)] <- 0
df_corr_XGBndcg_LKPO <- MeanR(df_corr_XGBndcg_LKPO[N > 4], N[N > 4], "MinVar")
df_corr_XGBndcg_LKAO <- df_corr[variable == "RankMHC" & CV_type == "LKAO", V1]
df_corr_XGBndcg_LKAO[is.na(df_corr_XGBndcg_LKAO)] <- 0
df_corr_XGBndcg_LKAO <- MeanR(df_corr_XGBndcg_LKAO[N > 4], N[N > 4], "MinVar")
df_corr_XGBndcg_LOLO <- df_corr[variable == "RankMHC" & CV_type == "LOLO", V1]
df_corr_XGBndcg_LOLO[is.na(df_corr_XGBndcg_LOLO)] <- 0
df_corr_XGBndcg_LOLO <- MeanR(df_corr_XGBndcg_LOLO[N > 4], N[N > 4], "MinVar")
df_corr_mean <- data.table(variable=rep(c("ref2015", "vina", "vinardo", "molPDF", "RankMHC"), 3),
                           CV_type=c(rep("LKPO", 5), rep("LKAO", 5), rep("LOLO", 5)),
                           metric=rep("ρ", 15),
                           value=c(df_corr_ref2015_LKPO, df_corr_vina_LKPO, df_corr_vinardo_LKPO, df_corr_molPDF_LKPO, df_corr_XGBndcg_LKPO,
                                   df_corr_ref2015_LKAO, df_corr_vina_LKAO, df_corr_vinardo_LKAO, df_corr_molPDF_LKAO, df_corr_XGBndcg_LKAO,
                                   df_corr_ref2015_LOLO, df_corr_vina_LOLO, df_corr_vinardo_LOLO, df_corr_molPDF_LOLO, df_corr_XGBndcg_LOLO))
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
df_precision[, `:=`(qid = NULL, value = NULL, min_rmsd = NULL, ref2015.x = NULL, ref2015.y = NULL)]
df_p1 <- df_p1[variable != "Best"]
df_p1[, metric := "LRMSD@1"]
df_p1[, `:=`(filename = NULL, V1 = min_rmsd, value = NULL, min_rmsd = NULL, ref2015 = NULL, vina = NULL, vinardo = NULL, molPDF = NULL)]
df_precision <- dcast.data.table(df_precision, ... ~ variable, value.var = c("HA_rmsd"))
df_precision <- df_precision[!is.na(Best)]
len_p <- df[, pdb_code] %>% unique %>% length
df_precision <- df_precision[, lapply(.SD, function(x) {sum(!is.na(x))/len_p}), by = CV_type, .SDcols = c("ref2015", "vina", "vinardo", "molPDF", "RankMHC", "Best")]
df_precision[, Best := NULL]
df_precision <- melt.data.table(df_precision, id.vars = c("CV_type"))
df_precision[, metric := "P@1"]
df_p2 <- df[order(pdb_code, value), head(.SD, 2), by = .(pdb_code, variable, CV_type), .SDcols = names(df)]
df_p2 [, min_rmsd := min(HA_rmsd), by = .(pdb_code, variable, CV_type)]
df_p2[, `:=`(pdb_code = NULL, variable = NULL, CV_type = NULL, HA_rmsd = NULL)]
df_p2 <- unique(df_p2, by = c("pdb_code", "qid", "peptide", "allele", "variable", "CV_type"))
df_p2 <- df_p2[variable != "Best"]
df_p2[, metric := "Best_RMSD@2"]
df_p2[, `:=`(filename = NULL, V1 = min_rmsd, value = NULL, min_rmsd = NULL, ref2015 = NULL, vina = NULL, vinardo = NULL, molPDF = NULL)]
df_p3 <- df[order(pdb_code, value), head(.SD, 3), by = .(pdb_code, variable, CV_type), .SDcols = names(df)]
df_p3 [, min_rmsd := min(HA_rmsd), by = .(pdb_code, variable, CV_type)]
df_p3[, `:=`(pdb_code = NULL, variable = NULL, CV_type = NULL, HA_rmsd = NULL)]
df_p3 <- unique(df_p3, by = c("pdb_code", "qid", "peptide", "allele", "variable", "CV_type"))
df_p3 <- df_p3[variable != "Best"]
df_p3[, metric := "Best_RMSD@3"]
df_p3[, `:=`(filename = NULL, V1 = min_rmsd, value = NULL, min_rmsd = NULL, ref2015 = NULL, vina = NULL, vinardo = NULL, molPDF = NULL)]
df_metric <- rbindlist(list(df_ndcg, df_p1, df_rr), use.names=TRUE)
df_metric[, metric := factor(metric, levels = c("NDCG", "LRMSD@1", "MRR"), ordered = T)]
df_metric[, value := V1]
df_metric[, V1 := NULL]
mean_data_openmm <- df_metric[, mean(`value`), by = .(variable, CV_type, metric)]
mean_data_openmm[, value := round(V1, digits = 3)]
mean_data_openmm[, V1 := NULL]
mean_data_openmm <- rbindlist(list(mean_data_openmm, df_corr_mean, df_precision), use.names = TRUE)
mean_data_openmm[, type := "PANDORA + HA + min"]

final_data <- rbindlist(list(mean_data[CV_type == "LKPO"], mean_data_pdbfixer[CV_type == "LKPO"], mean_data_openmm[CV_type == "LKPO"]))
final_data[, variable := factor(variable, levels = c("ref2015", "vina", "vinardo", "molPDF", "RankMHC"), ordered = T)]
final_data[, metric := factor(metric, levels = c("ρ", "MRR", "P@1", "NDCG", "LRMSD@1"), ordered = T)]
final_data[metric == "P@1", value := value*100]
ggplot(data = final_data[type != "PANDORA + HA"], aes(x = variable, y = value, fill = variable)) +
  geom_bar(position="dodge", stat="identity", width = 0.5) +
  scale_fill_viridis(discrete = TRUE) +
  ggh4x::facet_grid2(vars(type), vars(metric), scales="free_y", independent = "y") +
  xlab("Model") +
  ylab("Value") +
  theme_bw() +
  theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12), plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
        strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16),
        axis.title.y = element_text(size = rel(1.6)), axis.title.x = element_text(size = rel(1.6)), 
        axis.text.y = element_text(angle = 0, hjust = 1, size = 16))

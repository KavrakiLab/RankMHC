import numpy as np
import pandas as pd

import xgboost as xgb

from copy import deepcopy
import pickle as pkl
from tqdm import tqdm
from collections import defaultdict
import sys

from scipy.stats import spearmanr
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform

from mpire import WorkerPool

def permutation_importance(predictor, n_repeats, features_list, model_list, reference_spearman, reference_p1):
    df_test_list=[]
    spearman_repeats_list = []
    best_lrmsd_repeats_list = []
    for j in range(1, n_repeats + 1):
        for i in range(1, 7):
            ranker = model_list[i - 1]
            with open('./pickle_splits/LKPO_Outer_' + str(i) +'.pkl', 'rb') as f:
                outer_split_pdb_list = pkl.load(f)
            df_test = df[~df["pdb_code"].isin(outer_split_pdb_list)]
            df_test = pd.merge(df_test, non_redundand)
            y_test = df_test["HA_rmsd"].values
            df_test_test_ranker = df_test.drop(columns=column_drop)
            df_test_test_ranker = df_test_test_ranker[features_list]
            df_test_test_ranker_copy = deepcopy(df_test_test_ranker)
            df_test_test_ranker_copy[predictor] = df_test_test_ranker[predictor].sample(frac=1).values
            scores = ranker.predict(df_test_test_ranker_copy)
            df_test = df_test.assign(Score = scores)
            df_test_list.append(df_test)
        final_df = pd.concat(df_test_list)[['qid', 'pdb_code', 'filename', 'peptide', 'allele', 'HA_rmsd', 'Score']]
        spearman_list = []
        best_lrmsd_list = []
        split = [y for x, y in final_df.groupby('qid')]
        for split_df in split:
            spearman_list.append(spearmanr(split_df['HA_rmsd'], -split_df['Score']).correlation)
            best_lrmsd = split_df['HA_rmsd'].min()
            best_score_lrmsd = split_df[split_df['Score']==split_df['Score'].max()]['HA_rmsd'].values[0]
            best_lrmsd_list.append(best_lrmsd == best_score_lrmsd)
        spearman_repeats_list.append(np.mean(spearman_list))
        best_lrmsd_repeats_list.append(np.mean(best_lrmsd_list))
    return {'pred': predictor, 'spearman': reference_spearman - np.mean(spearman_repeats_list), 'p1': np.mean(best_lrmsd_repeats_list) - reference_p1}

# Clustering Threshold and repeats determination
threshold = float(sys.argv[1])
n_repeats = int(sys.argv[2])

# Load necessary files
rmsd_scores = pd.read_csv("./rmsd_scores.csv")
df = pd.read_csv("../Scoring_function_XGBoost/features/Rosetta_peptide_only.csv")
df = df.merge(rmsd_scores, how='inner', on=['pdb_code', 'filename'])
non_redundand = pd.read_csv("./Non_rendundant_structs.csv")

#Convert pdb_code to qid for ranking
df["qid"] = deepcopy(df["pdb_code"])
df["qid"] = pd.Categorical(df["qid"], ordered=True)
df["qid"] = df["qid"].cat.codes
cols = df.columns.tolist()
cols = cols[-1:] + cols[:-1]
df = df[cols]

# Get list of features
df.insert(5, "pep_len", df['peptide'].str.len())
column_drop = ['pdb_code', 'filename', 'allele', 'peptide', 'pep_len', 'CA_rmsd', 'BB_rmsd', 'HA_rmsd', 'dscore']

# Calculate correlation and select features based on a hierarchical clustering threshold
print("Feauture selection step")
features_list = []
for i in tqdm(range(1, 7)):
    with open('./pickle_splits/LKPO_Outer_' + str(i) +'.pkl', 'rb') as f:
        outer_split_pdb_list = pkl.load(f)
    df_train = df[df["pdb_code"].isin(outer_split_pdb_list)]
    df_train = df_train.drop(columns=column_drop)
    cols = df_train.columns
    corr = spearmanr(df_train.drop(columns=['qid'])).correlation

    # Ensure the correlation matrix is symmetric
    corr = (corr + corr.T) / 2
    np.fill_diagonal(corr, 1)
    indexes = np.argwhere(corr > 0.8).tolist()

    # Cluster and remove features
    distance_matrix = 1 - np.abs(corr)
    distance_matrix = np.nan_to_num(distance_matrix)
    dist_linkage = hierarchy.ward(squareform(distance_matrix, checks=False))
    cluster_ids = hierarchy.fcluster(dist_linkage, threshold, criterion="distance")
    cluster_id_to_feature_ids = defaultdict(list)
    for idx, cluster_id in enumerate(cluster_ids):
        cluster_id_to_feature_ids[cluster_id].append(idx)
    selected_features = [v[0] for v in cluster_id_to_feature_ids.values()]
    selected_features_names = df_train.drop(columns=['qid']).columns[selected_features]    
    features_list.append(selected_features_names.tolist())
features_list = [item for sublist in features_list for item in sublist]
features_list = list(set(features_list))
features_list.append('qid')

# Train model with the reduced feature set
print("Train models")
model_list = []
for i in tqdm(range(1, 7)):
    with open('./pickle_splits/LKPO_Outer_' + str(i) +'.pkl', 'rb') as f:
        outer_split_pdb_list = pkl.load(f)
    df_train = df[df["pdb_code"].isin(outer_split_pdb_list)]
    y_train = df_train["HA_rmsd"].values
    df_train = df_train.drop(columns=column_drop)
    df_train = df_train[features_list]
    ranker = xgb.XGBRanker(tree_method="hist", eta=0.01, gamma=1, max_depth=8, min_child_weight=1.0, max_delta_step=0, subsample=0.75, 
                           reg_lambda=1, alpha=1, lambdarank_num_pair_per_sample=10, objective='rank:pairwise',
                           lambdarank_pair_method="mean", nthread=8)
    ranker.fit(df_train, -y_train)
    model_list.append(ranker)

# Calculate reference scores
print("Calculating reference scores")
df_test_list = []
for i in tqdm(range(1, 7)):
    ranker = model_list[i - 1]
    with open('./pickle_splits/LKPO_Outer_' + str(i) +'.pkl', 'rb') as f:
        outer_split_pdb_list = pkl.load(f)
    df_test = df[~df["pdb_code"].isin(outer_split_pdb_list)]
    df_test = pd.merge(df_test, non_redundand)
    y_test = df_test["HA_rmsd"].values
    df_test_test_ranker = df_test.drop(columns=column_drop)
    df_test_test_ranker = df_test_test_ranker[features_list]
    scores = ranker.predict(df_test_test_ranker)
    df_test = df_test.assign(Score = scores)
    df_test_list.append(df_test)
final_df = pd.concat(df_test_list)[['qid', 'pdb_code', 'filename', 'peptide', 'allele', 'HA_rmsd', 'Score']]
spearman_list = []
best_lrmsd_list = []
split = [y for x, y in final_df.groupby('qid')]
for split_df in split:
    spearman_list.append(spearmanr(split_df['HA_rmsd'], -split_df['Score']).correlation)
    best_lrmsd = split_df['HA_rmsd'].min()
    best_score_lrmsd = split_df[split_df['Score']==split_df['Score'].max()]['HA_rmsd'].values[0]
    best_lrmsd_list.append(best_lrmsd == best_score_lrmsd)
reference_spearman = np.mean(spearman_list)
reference_p1 = np.mean(best_lrmsd_list)

arg_list = []
predictors = [item for item in features_list if item not in ['qid']] 
for predictor in predictors:
    arg_list.append((predictor, n_repeats, features_list, model_list, reference_spearman, reference_p1))

print("Calculating Permutation Importance")
results = []
for arg in tqdm(arg_list):
    results.append(permutation_importance(arg[0], arg[1], arg[2], arg[3], arg[4], arg[5]))
    
resultsdf = pd.DataFrame(results)
resultsdf.to_csv("models/permutation_importance_v2_" + str(threshold) + "_" + str(n_repeats) + ".csv", index = True)
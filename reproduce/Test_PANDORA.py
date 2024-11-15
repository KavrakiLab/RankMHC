import argparse
from tqdm import tqdm
import pickle as pkl
from copy import deepcopy
import itertools
import os

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import MinMaxScaler

import xgboost as xgb

def MHCscore_parser():

    parser = argparse.ArgumentParser(description="Hyperparameter Tuning parser", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--features', type=str, default='peptide_only', choices=['peptide_only', 'peptide_HLA', 'pairwise'], help='Feature file that includes different features to use')
    parser.add_argument('--padding', type=str, default='', choices=['x', 'oct', 'pad'], help='Feature file that includes different features to use')
    parser.add_argument("--labels", type=str, default='HA_rmsd', choices=['HA_rmsd', 'dscore'], help='The training label (RMSDs, D-score)')
    parser.add_argument("--reduced", type=str, default='', choices=['x', 'reduced'], help='Either use the full set (almost Graddock) of features or the ref2015 set (reduced)')
    parser.add_argument("--redundancy", type=str, default='', choices=['x', 'redundant'], help='Whether to include redundant structures for data augmentation or not.')
    parser.add_argument("--objective", type=str, default='rank:pairwise', choices=['rank:ndcg', 'rank:pairwise', 'regr'], help='The training objective (ranking, regression)')
    parser.add_argument("--num_cores", type=int, default=10, help='Number of cores to use for training/evaluation')
    parser.add_argument("--learning_rate_LKPO", type=float, default=0.2, help='Learning rate')
    parser.add_argument("--gamma_LKPO", type=float, default=0.0, help='Minimum loss reduction required to make a further partition on a leaf node of a tree')
    parser.add_argument("--max_depth_LKPO", type=int, default=8, help='Maximum depth of a tree. Increasing this value will make the model more complex and more likely to overfit.')
    parser.add_argument("--min_child_weight_LKPO", type=float, default=1, help='Higher values prevent a model from learning relations which might be highly specific to the particular sample selected for a tree.')
    parser.add_argument("--max_delta_step_LKPO", type=float, default=0, help='If this is set to a positive value, it can help making the update step more conservative')
    parser.add_argument("--subsample_LKPO", type=float, default=1, help='Setting it to 0.5 means that XGBoost would randomly sample half of the training data prior to growing trees')
    parser.add_argument("--lambda_l2_LKPO", type=float, default=1.0, help='L2 norm')
    parser.add_argument("--alpha_l1_LKPO", type=float, default=0, help='L1 norm')
    parser.add_argument("--lambda_num_pair_per_sample_LKPO", type=int, default=5, help='Ranking-only parameter')
    parser.add_argument("--learning_rate_LKAO", type=float, default=0.2, help='Learning rate')
    parser.add_argument("--gamma_LKAO", type=float, default=0.0, help='Minimum loss reduction required to make a further partition on a leaf node of a tree')
    parser.add_argument("--max_depth_LKAO", type=int, default=8, help='Maximum depth of a tree. Increasing this value will make the model more complex and more likely to overfit.')
    parser.add_argument("--min_child_weight_LKAO", type=float, default=1, help='Higher values prevent a model from learning relations which might be highly specific to the particular sample selected for a tree.')
    parser.add_argument("--max_delta_step_LKAO", type=float, default=0, help='If this is set to a positive value, it can help making the update step more conservative')
    parser.add_argument("--subsample_LKAO", type=float, default=1, help='Setting it to 0.5 means that XGBoost would randomly sample half of the training data prior to growing trees')
    parser.add_argument("--lambda_l2_LKAO", type=float, default=1.0, help='L2 norm')
    parser.add_argument("--alpha_l1_LKAO", type=float, default=0, help='L1 norm')
    parser.add_argument("--lambda_num_pair_per_sample_LKAO", type=int, default=5, help='Ranking-only parameter')
    parser.add_argument("--learning_rate_LOLO", type=float, default=0.2, help='Learning rate')
    parser.add_argument("--gamma_LOLO", type=float, default=0.0, help='Minimum loss reduction required to make a further partition on a leaf node of a tree')
    parser.add_argument("--max_depth_LOLO", type=int, default=8, help='Maximum depth of a tree. Increasing this value will make the model more complex and more likely to overfit.')
    parser.add_argument("--min_child_weight_LOLO", type=float, default=1, help='Higher values prevent a model from learning relations which might be highly specific to the particular sample selected for a tree.')
    parser.add_argument("--max_delta_step_LOLO", type=float, default=0, help='If this is set to a positive value, it can help making the update step more conservative')
    parser.add_argument("--subsample_LOLO", type=float, default=1, help='Setting it to 0.5 means that XGBoost would randomly sample half of the training data prior to growing trees')
    parser.add_argument("--lambda_l2_LOLO", type=float, default=1.0, help='L2 norm')
    parser.add_argument("--alpha_l1_LOLO", type=float, default=0, help='L1 norm')
    parser.add_argument("--lambda_num_pair_per_sample_LOLO", type=int, default=5, help='Ranking-only parameter')
    return parser

parser = MHCscore_parser()
args = parser.parse_args()
features = args.features
padding = args.padding
RMSD_type = args.labels
reduced_features = args.reduced
redundancy = args.redundancy
objective = args.objective
num_cores = args.num_cores
learning_rate_LKPO = args.learning_rate_LKPO
gamma_LKPO = args.gamma_LKPO
max_depth_LKPO = args.max_depth_LKPO
min_child_weight_LKPO = args.min_child_weight_LKPO
max_delta_step_LKPO = args.max_delta_step_LKPO
subsample_LKPO = args.subsample_LKPO
lambda_l2_LKPO = args.lambda_l2_LKPO
alpha_l1_LKPO = args.alpha_l1_LKPO
lambda_num_pair_per_sample_LKPO = args.lambda_num_pair_per_sample_LKPO
learning_rate_LKAO = args.learning_rate_LKAO
gamma_LKAO = args.gamma_LKAO
max_depth_LKAO = args.max_depth_LKAO
min_child_weight_LKAO = args.min_child_weight_LKAO
max_delta_step_LKAO = args.max_delta_step_LKAO
subsample_LKAO = args.subsample_LKAO
lambda_l2_LKAO = args.lambda_l2_LKAO
alpha_l1_LKAO = args.alpha_l1_LKAO
lambda_num_pair_per_sample_LKAO = args.lambda_num_pair_per_sample_LKAO
learning_rate_LOLO = args.learning_rate_LOLO
gamma_LOLO = args.gamma_LOLO
max_depth_LOLO = args.max_depth_LOLO
min_child_weight_LOLO = args.min_child_weight_LOLO
max_delta_step_LOLO = args.max_delta_step_LOLO
subsample_LOLO = args.subsample_LOLO
lambda_l2_LOLO = args.lambda_l2_LOLO
alpha_l1_LOLO = args.alpha_l1_LOLO
lambda_num_pair_per_sample_LOLO = args.lambda_num_pair_per_sample_LOLO
effort = features + '_' + RMSD_type + '_' + reduced_features + '_' + padding + '_' + redundancy + '_' + objective

# Define input/output files
feature_file = './features/Rosetta_'
feature_file += features
if padding == 'pad':
    feature_file += '_' + padding
if reduced_features == 'reduced':
    feature_file += '_' + reduced_features
log_file = feature_file + '.log'
feature_file += '.csv'
res_file = feature_file.split("/")[2].split(".")[0] + "_" + RMSD_type
if redundancy == 'redundant':
    res_file += '_redundant'
res_file += '_res.csv'

rmsd_scores = pd.read_csv("./rmsd_scores.csv")
Pandora_rmsd_scores = pd.read_csv("./Pandora_features/Pandora_rmsd_scores_openmm.csv")

# Load the feature file
df = pd.read_csv(feature_file)
df = df.merge(rmsd_scores, how='inner', on=['pdb_code', 'filename'])

Pandora_df = pd.read_csv("./Pandora_features/Pandoraopenmm_peptide_only.csv")
Pandora_df = Pandora_df.merge(Pandora_rmsd_scores, how='inner', on=['pdb_code', 'filename'])

# Filter out if you want to include redundand structures:
non_redundand = pd.read_csv("./Non_rendundant_structs.csv")
if redundancy != 'redundant':
    df = pd.merge(df, non_redundand)

#Convert pdb_code to qid for ranking
df["qid"] = deepcopy(df["pdb_code"])
df["qid"] = pd.Categorical(df["qid"], ordered=True)
df["qid"] = df["qid"].cat.codes
cols = df.columns.tolist()
cols = cols[-1:] + cols[:-1]
df = df[cols]

Pandora_df["qid"] = deepcopy(Pandora_df["pdb_code"])
Pandora_df["qid"] = pd.Categorical(Pandora_df["qid"], ordered=True)
Pandora_df["qid"] = Pandora_df["qid"].cat.codes
cols = Pandora_df.columns.tolist()
cols = cols[-1:] + cols[:-1]
Pandora_df = Pandora_df[cols]


# Add peptide length
df.insert(5, "pep_len", df['peptide'].str.len())

column_drop = ['pdb_code', 'filename', 'allele', 'peptide', 'pep_len', 'CA_rmsd', 'BB_rmsd', 'HA_rmsd', 'dscore']

Pandora_df.insert(5, "pep_len", Pandora_df['peptide'].str.len())
print(Pandora_df)
Pandora_column_drop = ['pdb_code', 'filename', 'allele', 'peptide', 'pep_len', 'CA_rmsd', 'BB_rmsd', 'HA_rmsd']

# LEAVE K PDBs OUT!
df_test_list = []
for outer in tqdm(range(1, 7), position=0):
    with open('./pickle_splits/LKPO_Outer_' + str(outer) +'.pkl', 'rb') as f:
        outer_split_pdb_list = pkl.load(f)
    df_train = df[df["pdb_code"].isin(outer_split_pdb_list)]
    y_train = df_train[RMSD_type].values
    df_train = df_train.drop(columns=column_drop)
    df_test = Pandora_df[~Pandora_df["pdb_code"].isin(outer_split_pdb_list)]
    y_test = df_test[RMSD_type].values
    df_test_ranker = df_test.drop(columns=Pandora_column_drop)
    if objective == 'regr':
        ranker = xgb.XGBRegressor(tree_method="hist", eta=learning_rate_LKPO, gamma=gamma_LKPO, max_depth=max_depth_LKPO,
                                  min_child_weight=min_child_weight_LKPO, max_delta_step=max_delta_step_LKPO,
                                  subsample=subsample_LKPO, reg_lambda=lambda_l2_LKPO, alpha=alpha_l1_LKPO, nthread=num_cores)
        ranker.fit(df_train, y_train)
    elif objective == 'rank:pairwise':
        ranker = xgb.XGBRanker(tree_method="hist", eta=learning_rate_LKPO, gamma=gamma_LKPO, max_depth=max_depth_LKPO,
                               min_child_weight=min_child_weight_LKPO, max_delta_step=max_delta_step_LKPO,
                               subsample=subsample_LKPO, reg_lambda=lambda_l2_LKPO, alpha=alpha_l1_LKPO,
                               lambdarank_num_pair_per_sample=lambda_num_pair_per_sample_LKPO, objective=objective,
                               lambdarank_pair_method="mean", nthread=num_cores)
        ranker.fit(df_train, -y_train)
    else:
        scaler = MinMaxScaler()
        y_train_01 = scaler.fit_transform((-y_train).reshape(-1, 1))
        ranker = xgb.XGBRanker(tree_method="hist", eta=learning_rate_LKPO, gamma=gamma_LKPO, max_depth=max_depth_LKPO,
                               min_child_weight=min_child_weight_LKPO, max_delta_step=max_delta_step_LKPO,
                               subsample=subsample_LKPO, reg_lambda=lambda_l2_LKPO, alpha=alpha_l1_LKPO,
                               lambdarank_num_pair_per_sample=lambda_num_pair_per_sample_LKPO, objective=objective,
                               lambdarank_pair_method="mean", ndcg_exp_gain=False, nthread=num_cores)
        ranker.fit(df_train, y_train_01)
    df_test_test_ranker = df_test.drop(columns=Pandora_column_drop)
    scores = ranker.predict(df_test_test_ranker)
    df_test = df_test.assign(Score = scores)
    df_test_list.append(df_test)
final_df = pd.concat(df_test_list)[['qid', 'pdb_code', 'filename', 'peptide', 'allele', 'Score']]
final_df.to_csv("./Pandora_features/res_v5/LKPO_Pandora_openmm.csv")

# LEAVE K ALLELES OUT!
df_test_list = []
for outer in tqdm(range(1, 7), position=0):
    with open('./pickle_splits/LKAO_Outer_' + str(outer) +'.pkl', 'rb') as f:
        outer_split_pdb_list = pkl.load(f)
    df_train = df[df["allele"].isin(outer_split_pdb_list)]
    y_train = df_train[RMSD_type].values
    df_train = df_train.drop(columns=column_drop)
    df_test = Pandora_df[~Pandora_df["allele"].isin(outer_split_pdb_list)]
    y_test = df_test[RMSD_type].values
    df_test_ranker = df_test.drop(columns=Pandora_column_drop)
    if objective == 'regr':
        ranker = xgb.XGBRegressor(tree_method="hist", eta=learning_rate_LKAO, gamma=gamma_LKAO, max_depth=max_depth_LKAO,
                                  min_child_weight=min_child_weight_LKAO, max_delta_step=max_delta_step_LKAO,
                                  subsample=subsample_LKAO, reg_lambda=lambda_l2_LKAO, alpha=alpha_l1_LKAO, nthread=num_cores)
        ranker.fit(df_train, y_train)
    elif objective == 'rank:pairwise':
        ranker = xgb.XGBRanker(tree_method="hist", eta=learning_rate_LKAO, gamma=gamma_LKAO, max_depth=max_depth_LKAO,
                               min_child_weight=min_child_weight_LKAO, max_delta_step=max_delta_step_LKAO,
                               subsample=subsample_LKAO, reg_lambda=lambda_l2_LKAO, alpha=alpha_l1_LKAO,
                               lambdarank_num_pair_per_sample=lambda_num_pair_per_sample_LKAO, objective=objective,
                               lambdarank_pair_method="mean", nthread=num_cores)
        ranker.fit(df_train, -y_train)
    else:
        scaler = MinMaxScaler()
        y_train_01 = scaler.fit_transform((-y_train).reshape(-1, 1))
        ranker = xgb.XGBRanker(tree_method="hist", eta=learning_rate_LKAO, gamma=gamma_LKAO, max_depth=max_depth_LKAO,
                               min_child_weight=min_child_weight_LKAO, max_delta_step=max_delta_step_LKAO,
                               subsample=subsample_LKAO, reg_lambda=lambda_l2_LKAO, alpha=alpha_l1_LKAO,
                               lambdarank_num_pair_per_sample=lambda_num_pair_per_sample_LKAO, objective=objective,
                               lambdarank_pair_method="mean", ndcg_exp_gain=False, nthread=num_cores)
        ranker.fit(df_train, y_train_01)
    df_test_test_ranker = df_test.drop(columns=Pandora_column_drop)
    scores = ranker.predict(df_test_test_ranker)
    df_test = df_test.assign(Score = scores)
    df_test_list.append(df_test)
final_df = pd.concat(df_test_list)[['qid', 'pdb_code', 'filename', 'peptide', 'allele', 'Score']]
final_df.to_csv("./Pandora_features/res_v5/LKAO_Pandora_openmm.csv")

# LEAVE ONE LENGTH OUT!
df_test_list = []
for outer in tqdm(range(1, 7), position=0):
    with open('./pickle_splits/LOLO_Outer_' + str(outer) +'.pkl', 'rb') as f:
        outer_split_pdb_list = pkl.load(f)
    df_train = df[df["pep_len"].isin(outer_split_pdb_list)]
    y_train = df_train[RMSD_type].values
    df_train = df_train.drop(columns=column_drop)
    df_test = Pandora_df[~Pandora_df["pep_len"].isin(outer_split_pdb_list)]
    y_test = df_test[RMSD_type].values
    df_test_ranker = df_test.drop(columns=Pandora_column_drop)
    if objective == 'regr':
        ranker = xgb.XGBRegressor(tree_method="hist", eta=learning_rate_LOLO, gamma=gamma_LOLO, max_depth=max_depth_LOLO,
                                  min_child_weight=min_child_weight_LOLO, max_delta_step=max_delta_step_LOLO,
                                  subsample=subsample_LOLO, reg_lambda=lambda_l2_LOLO, alpha=alpha_l1_LOLO, nthread=num_cores)
        ranker.fit(df_train, y_train)
    elif objective == 'rank:pairwise':
        ranker = xgb.XGBRanker(tree_method="hist", eta=learning_rate_LOLO, gamma=gamma_LOLO, max_depth=max_depth_LOLO,
                               min_child_weight=min_child_weight_LOLO, max_delta_step=max_delta_step_LOLO,
                               subsample=subsample_LOLO, reg_lambda=lambda_l2_LOLO, alpha=alpha_l1_LOLO,
                               lambdarank_num_pair_per_sample=lambda_num_pair_per_sample_LOLO, objective=objective,
                               lambdarank_pair_method="mean", nthread=num_cores)
        ranker.fit(df_train, -y_train)
    else:
        scaler = MinMaxScaler()
        y_train_01 = scaler.fit_transform((-y_train).reshape(-1, 1))
        ranker = xgb.XGBRanker(tree_method="hist", eta=learning_rate_LOLO, gamma=gamma_LOLO, max_depth=max_depth_LOLO,
                               min_child_weight=min_child_weight_LOLO, max_delta_step=max_delta_step_LOLO,
                               subsample=subsample_LOLO, reg_lambda=lambda_l2_LOLO, alpha=alpha_l1_LOLO,
                               lambdarank_num_pair_per_sample=lambda_num_pair_per_sample_LOLO, objective=objective,
                               lambdarank_pair_method="mean", ndcg_exp_gain=False, nthread=num_cores)
        ranker.fit(df_train, y_train_01)
    df_test_test_ranker = df_test.drop(columns=Pandora_column_drop)
    scores = ranker.predict(df_test_test_ranker)
    df_test = df_test.assign(Score = scores)
    df_test_list.append(df_test)
final_df = pd.concat(df_test_list)[['qid', 'pdb_code', 'filename', 'peptide', 'allele', 'Score']]
final_df.to_csv("./Pandora_features/res_v5/LOLO_Pandora_openmm.csv")

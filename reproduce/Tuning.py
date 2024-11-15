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
	parser.add_argument("--learning_rate", type=float, default=0.2, help='Learning rate')
	parser.add_argument("--gamma", type=float, default=0.0, help='Minimum loss reduction required to make a further partition on a leaf node of a tree')
	parser.add_argument("--max_depth", type=int, default=8, help='Maximum depth of a tree. Increasing this value will make the model more complex and more likely to overfit.')
	parser.add_argument("--min_child_weight", type=float, default=1, help='Higher values prevent a model from learning relations which might be highly specific to the particular sample selected for a tree.')
	parser.add_argument("--max_delta_step", type=float, default=0, help='If this is set to a positive value, it can help making the update step more conservative')
	parser.add_argument("--subsample", type=float, default=1, help='Setting it to 0.5 means that XGBoost would randomly sample half of the training data prior to growing trees')
	parser.add_argument("--lambda_l2", type=float, default=1.0, help='L2 norm')
	parser.add_argument("--alpha_l1", type=float, default=0, help='L1 norm')
	parser.add_argument("--lambda_num_pair_per_sample", type=int, default=5, help='Ranking-only parameter')
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
learning_rate = args.learning_rate
gamma = args.gamma
max_depth = args.max_depth
min_child_weight = args.min_child_weight
max_delta_step = args.max_delta_step
subsample = args.subsample
lambda_l2 = args.lambda_l2
alpha_l1 = args.alpha_l1
lambda_num_pair_per_sample = args.lambda_num_pair_per_sample
effort = features + '_' + RMSD_type + '_' + reduced_features + '_' + padding + '_' + redundancy + '_' + objective + '_lr' + \
         str(learning_rate) + '_gamma' + str(gamma) + '_max_depth' + str(max_depth) + '_min_child_weight' + \
         str(min_child_weight) + '_max_delta_step' + str(max_delta_step) + '_subsample' + str(subsample) + \
         '_l2norm' + str(lambda_l2) + '_l1norm' + str(alpha_l1) + \
         '_lambda_num_pair_per_sample' + str(lambda_num_pair_per_sample)

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

# Load the feature file
df = pd.read_csv(feature_file)
df = df.merge(rmsd_scores, how='inner', on=['pdb_code', 'filename'])

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

# Add peptide length
df.insert(5, "pep_len", df['peptide'].str.len())

column_drop = ['pdb_code', 'filename', 'allele', 'peptide', 'pep_len', 'CA_rmsd', 'BB_rmsd', 'HA_rmsd', 'dscore']

# LEAVE K PDBs OUT!
outer_spearman_LKPO_list = []
outer_min_diff_LKPO_list = [] 
for outer in tqdm(range(1, 7), position=0):
	df_test_list = []
	peptide_list = []
	MHC_list = []
	min_diff_list = []
	spearman_list = []
	with open('./pickle_splits/LKPO_Outer_' + str(outer) +'.pkl', 'rb') as f:
		outer_split_pdb_list = pkl.load(f)
	df_outer = df[df["pdb_code"].isin(outer_split_pdb_list)]
	for inner in tqdm(range(1, 6), position=1, leave=True):
		with open('./pickle_splits/LKPO_Inner_' + str(outer) + '_' + str(inner) +'.pkl', 'rb') as f:
			inner_split_pdb_list = pkl.load(f)
		df_train = df_outer[df_outer["pdb_code"].isin(inner_split_pdb_list)]
		y_train = df_train[RMSD_type].values
		df_train = df_train.drop(columns=column_drop)
		df_test = df_outer[~df_outer["pdb_code"].isin(inner_split_pdb_list)]
		df_test = pd.merge(df_test, non_redundand)
		y_test = df_test[RMSD_type].values
		df_test_ranker = df_test.drop(columns=column_drop)
		if objective == 'regr':
			ranker = xgb.XGBRegressor(tree_method="hist", eta=learning_rate, gamma=gamma, max_depth=max_depth,
                                      min_child_weight=min_child_weight, max_delta_step=max_delta_step,
                                      subsample=subsample, reg_lambda=lambda_l2, alpha=alpha_l1, nthread=num_cores)
			ranker.fit(df_train, y_train)
		elif objective == 'rank:pairwise':
			ranker = xgb.XGBRanker(tree_method="hist", eta=learning_rate, gamma=gamma, max_depth=max_depth,
                                   min_child_weight=min_child_weight, max_delta_step=max_delta_step,
                                   subsample=subsample, reg_lambda=lambda_l2, alpha=alpha_l1,
                                   lambdarank_num_pair_per_sample=lambda_num_pair_per_sample, objective=objective,
                                   lambdarank_pair_method="mean", nthread=num_cores)
			ranker.fit(df_train, -y_train)
		else:
			scaler = MinMaxScaler()
			y_train_01 = scaler.fit_transform((-y_train).reshape(-1, 1))
			ranker = xgb.XGBRanker(tree_method="hist", eta=learning_rate, gamma=gamma, max_depth=max_depth,
                                   min_child_weight=min_child_weight, max_delta_step=max_delta_step,
                                   subsample=subsample, reg_lambda=lambda_l2, alpha=alpha_l1,
                                   lambdarank_num_pair_per_sample=lambda_num_pair_per_sample, objective=objective,
                                   lambdarank_pair_method="mean", ndcg_exp_gain=False, nthread=num_cores)
			ranker.fit(df_train, y_train_01)
		df_test_test_ranker = df_test.drop(columns=column_drop)
		scores = ranker.predict(df_test_test_ranker)
		df_test = df_test.assign(Score = scores)
		for idx in pd.unique(df_test["qid"]):
			df_test_test = df_test[df_test["qid"] == idx]
			peptide_list.append(pd.unique(df_test_test['peptide'])[0])
			MHC_list.append(pd.unique(df_test_test['allele'])[0])
			if objective == 'regr':
				min_diff_list.append(df_test_test[RMSD_type].values[np.argmin(df_test_test['Score'])] - df_test_test[RMSD_type].min())
				res = stats.spearmanr(df_test_test[RMSD_type], df_test_test['Score']).correlation
				if (np.isnan(res)):
					spearman_list.append(1/df_test_test[RMSD_type].shape[0])
				else:
					spearman_list.append(res)
			else:
				min_diff_list.append(df_test_test[RMSD_type].values[np.argmax(df_test_test['Score'])] - df_test_test[RMSD_type].min())
				res = stats.spearmanr(df_test_test[RMSD_type], -df_test_test['Score']).correlation
				if (np.isnan(res)):
					spearman_list.append(1/df_test_test[RMSD_type].shape[0])
				else:
					spearman_list.append(res)
	outer_spearman_LKPO_list.append(np.mean(spearman_list))
	outer_min_diff_LKPO_list.append(np.mean(min_diff_list))

print(outer_spearman_LKPO_list)
print(outer_min_diff_LKPO_list)

# LEAVE K ALLELES OUT!
outer_spearman_LKAO_list = []
outer_min_diff_LKAO_list = [] 
for outer in tqdm(range(1, 7), position=0):
	df_test_list = []
	peptide_list = []
	MHC_list = []
	min_diff_list = []
	spearman_list = []
	with open('./pickle_splits/LKAO_Outer_' + str(outer) +'.pkl', 'rb') as f:
		outer_split_pdb_list = pkl.load(f)
	df_outer = df[df["allele"].isin(outer_split_pdb_list)]
	for inner in tqdm(range(1, 6), position=1, leave=True):
		with open('./pickle_splits/LKAO_Inner_' + str(outer) + '_' + str(inner) +'.pkl', 'rb') as f:
			inner_split_pdb_list = pkl.load(f)
		df_train = df_outer[df_outer["allele"].isin(inner_split_pdb_list)]
		y_train = df_train[RMSD_type].values
		df_train = df_train.drop(columns=column_drop)
		df_test = df_outer[~df_outer["allele"].isin(inner_split_pdb_list)]
		df_test = pd.merge(df_test, non_redundand)
		y_test = df_test[RMSD_type].values
		df_test_ranker = df_test.drop(columns=column_drop)
		if objective == 'regr':
			ranker = xgb.XGBRegressor(tree_method="hist", eta=learning_rate, gamma=gamma, max_depth=max_depth,
                                      min_child_weight=min_child_weight, max_delta_step=max_delta_step,
                                      subsample=subsample, reg_lambda=lambda_l2, alpha=alpha_l1, nthread=num_cores)
			ranker.fit(df_train, y_train)
		elif objective == 'rank:pairwise':
			ranker = xgb.XGBRanker(tree_method="hist", eta=learning_rate, gamma=gamma, max_depth=max_depth,
                                   min_child_weight=min_child_weight, max_delta_step=max_delta_step,
                                   subsample=subsample, reg_lambda=lambda_l2, alpha=alpha_l1,
                                   lambdarank_num_pair_per_sample=lambda_num_pair_per_sample, objective=objective,
                                   lambdarank_pair_method="mean", nthread=num_cores)
			ranker.fit(df_train, -y_train)
		else:
			scaler = MinMaxScaler()
			y_train_01 = scaler.fit_transform((-y_train).reshape(-1, 1))
			ranker = xgb.XGBRanker(tree_method="hist", eta=learning_rate, gamma=gamma, max_depth=max_depth,
                                   min_child_weight=min_child_weight, max_delta_step=max_delta_step,
                                   subsample=subsample, reg_lambda=lambda_l2, alpha=alpha_l1,
                                   lambdarank_num_pair_per_sample=lambda_num_pair_per_sample, objective=objective,
                                   lambdarank_pair_method="mean", ndcg_exp_gain=False, nthread=num_cores)
			ranker.fit(df_train, y_train_01)
		df_test_test_ranker = df_test.drop(columns=column_drop)
		scores = ranker.predict(df_test_test_ranker)
		df_test = df_test.assign(Score = scores)
		for idx in pd.unique(df_test["qid"]):
			df_test_test = df_test[df_test["qid"] == idx]
			peptide_list.append(pd.unique(df_test_test['peptide'])[0])
			MHC_list.append(pd.unique(df_test_test['allele'])[0])
			if objective == 'regr':
				min_diff_list.append(df_test_test[RMSD_type].values[np.argmin(df_test_test['Score'])] - df_test_test[RMSD_type].min())
				res = stats.spearmanr(df_test_test[RMSD_type], df_test_test['Score']).correlation
				if (np.isnan(res)):
					spearman_list.append(1/df_test_test[RMSD_type].shape[0])
				else:
					spearman_list.append(res)
			else:
				min_diff_list.append(df_test_test[RMSD_type].values[np.argmax(df_test_test['Score'])] - df_test_test[RMSD_type].min())
				res = stats.spearmanr(df_test_test[RMSD_type], -df_test_test['Score']).correlation
				if (np.isnan(res)):
					spearman_list.append(1/df_test_test[RMSD_type].shape[0])
				else:
					spearman_list.append(res)
	outer_spearman_LKAO_list.append(np.mean(spearman_list))
	outer_min_diff_LKAO_list.append(np.mean(min_diff_list))

print(outer_spearman_LKAO_list)
print(outer_min_diff_LKAO_list)

# LEAVE ONE LENGTH OUT!
outer_spearman_LOLO_list = []
outer_min_diff_LOLO_list = [] 
for outer in tqdm(range(1, 7), position=0):
	df_test_list = []
	peptide_list = []
	MHC_list = []
	min_diff_list = []
	spearman_list = []
	with open('./pickle_splits/LOLO_Outer_' + str(outer) +'.pkl', 'rb') as f:
		outer_split_pdb_list = pkl.load(f)
	df_outer = df[df["pep_len"].isin(outer_split_pdb_list)]
	for inner in tqdm(range(1, 6), position=1, leave=True):
		with open('./pickle_splits/LOLO_Inner_' + str(outer) + '_' + str(inner) +'.pkl', 'rb') as f:
			inner_split_pdb_list = pkl.load(f)
		df_train = df_outer[df_outer["pep_len"].isin(inner_split_pdb_list)]
		y_train = df_train[RMSD_type].values
		df_train = df_train.drop(columns=column_drop)
		df_test = df_outer[~df_outer["pep_len"].isin(inner_split_pdb_list)]
		df_test = pd.merge(df_test, non_redundand)
		y_test = df_test[RMSD_type].values
		df_test_ranker = df_test.drop(columns=column_drop)
		if objective == 'regr':
			ranker = xgb.XGBRegressor(tree_method="hist", eta=learning_rate, gamma=gamma, max_depth=max_depth,
                                      min_child_weight=min_child_weight, max_delta_step=max_delta_step,
                                      subsample=subsample, reg_lambda=lambda_l2, alpha=alpha_l1, nthread=num_cores)
			ranker.fit(df_train, y_train)
		elif objective == 'rank:pairwise':
			ranker = xgb.XGBRanker(tree_method="hist", eta=learning_rate, gamma=gamma, max_depth=max_depth,
                                   min_child_weight=min_child_weight, max_delta_step=max_delta_step,
                                   subsample=subsample, reg_lambda=lambda_l2, alpha=alpha_l1,
                                   lambdarank_num_pair_per_sample=lambda_num_pair_per_sample, objective=objective,
                                   lambdarank_pair_method="mean", nthread=num_cores)
			ranker.fit(df_train, -y_train)
		else:
			scaler = MinMaxScaler()
			y_train_01 = scaler.fit_transform((-y_train).reshape(-1, 1))
			ranker = xgb.XGBRanker(tree_method="hist", eta=learning_rate, gamma=gamma, max_depth=max_depth,
                                   min_child_weight=min_child_weight, max_delta_step=max_delta_step,
                                   subsample=subsample, reg_lambda=lambda_l2, alpha=alpha_l1,
                                   lambdarank_num_pair_per_sample=lambda_num_pair_per_sample, objective=objective,
                                   lambdarank_pair_method="mean", ndcg_exp_gain=False, nthread=num_cores)
			ranker.fit(df_train, y_train_01)
		df_test_test_ranker = df_test.drop(columns=column_drop)
		scores = ranker.predict(df_test_test_ranker)
		df_test = df_test.assign(Score = scores)
		for idx in pd.unique(df_test["qid"]):
			df_test_test = df_test[df_test["qid"] == idx]
			peptide_list.append(pd.unique(df_test_test['peptide'])[0])
			MHC_list.append(pd.unique(df_test_test['allele'])[0])
			if objective == 'regr':
				min_diff_list.append(df_test_test[RMSD_type].values[np.argmin(df_test_test['Score'])] - df_test_test[RMSD_type].min())
				res = stats.spearmanr(df_test_test[RMSD_type], df_test_test['Score']).correlation
				if (np.isnan(res)):
					spearman_list.append(1/df_test_test[RMSD_type].shape[0])
				else:
					spearman_list.append(res)
			else:
				min_diff_list.append(df_test_test[RMSD_type].values[np.argmax(df_test_test['Score'])] - df_test_test[RMSD_type].min())
				res = stats.spearmanr(df_test_test[RMSD_type], -df_test_test['Score']).correlation
				if (np.isnan(res)):
					spearman_list.append(1/df_test_test[RMSD_type].shape[0])
				else:
					spearman_list.append(res)
	outer_spearman_LOLO_list.append(np.mean(spearman_list))
	outer_min_diff_LOLO_list.append(np.mean(min_diff_list))

print(outer_spearman_LOLO_list)
print(outer_min_diff_LOLO_list)

res = ",".join((features, RMSD_type, reduced_features, padding, redundancy, objective, str(learning_rate), str(gamma),
			    str(max_depth), str(min_child_weight), str(max_delta_step), str(subsample), str(lambda_l2), str(alpha_l1),
			    str(lambda_num_pair_per_sample), str(outer_spearman_LKPO_list[0]), str(outer_spearman_LKPO_list[1]), 
			    str(outer_spearman_LKPO_list[2]), str(outer_spearman_LKPO_list[3]), str(outer_spearman_LKPO_list[4]), 
			    str(outer_spearman_LKPO_list[5]), str(outer_min_diff_LKPO_list[0]), str(outer_min_diff_LKPO_list[1]), 
			    str(outer_min_diff_LKPO_list[2]), str(outer_min_diff_LKPO_list[3]), str(outer_min_diff_LKPO_list[4]), 
			    str(outer_min_diff_LKPO_list[5]), str(outer_spearman_LKAO_list[0]), str(outer_spearman_LKAO_list[1]), 
			    str(outer_spearman_LKAO_list[2]), str(outer_spearman_LKAO_list[3]), str(outer_spearman_LKAO_list[4]), 
			    str(outer_spearman_LKAO_list[5]), str(outer_min_diff_LKAO_list[0]), str(outer_min_diff_LKAO_list[1]), 
			    str(outer_min_diff_LKAO_list[2]), str(outer_min_diff_LKAO_list[3]), str(outer_min_diff_LKAO_list[4]), 
			    str(outer_min_diff_LKAO_list[5]), str(outer_spearman_LOLO_list[0]), str(outer_spearman_LOLO_list[1]), 
			    str(outer_spearman_LOLO_list[2]), str(outer_spearman_LOLO_list[3]), str(outer_spearman_LOLO_list[4]), 
			    str(outer_spearman_LOLO_list[5]), str(outer_min_diff_LOLO_list[0]), str(outer_min_diff_LOLO_list[1]), 
			    str(outer_min_diff_LOLO_list[2]), str(outer_min_diff_LOLO_list[3]), str(outer_min_diff_LOLO_list[4]), 
			    str(outer_min_diff_LOLO_list[5])))

with open('./resht/' + effort + ".txt","w+") as f:
	f.write(res)

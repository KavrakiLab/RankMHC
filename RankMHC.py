from tqdm import tqdm
import pickle as pkl
import os
import sys
from mpire import WorkerPool

import numpy as np
import pandas as pd

import xgboost as xgb

from feature_extraction.template_identification import find_template
from feature_extraction.feature_extraction import extract_features, MHCFlurry_padding, calculate_sasa
from helper_scripts import argparser
from helper_scripts.RankMHC_helper_funcs import peptide_processing, MHC_processing, initialize_dir, copy_batch_of_files, pretty_print_analytics, sanitize_PDBs, determine_sequence_and_allotype

# 0. ARGUMENTS
parser = argparser.RankMHC_parser()
args = parser.parse_args()
structure_input = args.input[0] + '/'
features = args.features
feat_transform = args.feat_transform
feat_terms = args.feat_terms
redundancy = args.redundancy
objective = args.objective
sasa_method = args.sasa_method
num_cores = args.num_cores
peptide = args.peptide
MHC = args.MHC

# 0.1 Directory to store any intermediate files
filestore = args.dir
initialize_dir(filestore)

# 1. SANITIZE AND CHECK INPUT
initialize_dir([filestore + '/input/', filestore + '/sanitized_input/', filestore + '/fasta_files/'])

# 1.1 Check input

allowed_models_dir = {('peptide_only', 'pooling', 'full', 'redundant', 'pairwise'): 'RankMHC', ('peptide_only', 'pooling', 'full', 'redundant', 'ndcg'): 'other_models/RankMHC_ndcg', 
					  ('peptide_only', 'pooling', 'full', 'redundant', 'regr'): 'other_models/RankMHC_regr', ('peptide_only', 'pooling', 'reduced', 'redundant', 'pairwise'): 'other_models/RankMHC_reduced', 
					  ('peptide_only', 'pooling', 'reduced', 'nonredundant', 'pairwise'): 'other_models/RankMHC_nonredundant', ('peptide_only', 'padding', 'full', 'redundant', 'pairwise'): 'other_models/RankMHC_padding', 
				  	  ('peptide_HLA', 'pooling', 'full', 'redundant', 'pairwise'): 'other_models/RankMHC_HLA', ('pairwise', 'pooling', 'full', 'redundant', 'pairwise'): 'other_models/RankMHC_pairwise'}
RankMHC_model = None
for key, value in allowed_models_dir.items():
	if (features, feat_transform, feat_terms, redundancy, objective) == key:
		RankMHC_model = value
if not RankMHC_model:
	print("RankMHC model request does not exist based on the input passed. Check your arguments! Aborting...")
	sys.exit(0) 				  

# 1.2 Determine peptide-MHC info based on the input
copy_batch_of_files(structure_input,
					filestore + '/input/',
					query="")
sanitize_PDBs(filestore + '/input/', filestore)
if peptide == '' or MHC == '':
	print("Either peptide or MHC fields were empty! Now trying to determine peptide sequence and MHC allotype...")
	(peptide, MHC) = determine_sequence_and_allotype(filestore + '/sanitized_input/', filestore)
	print("Peptide and MHC allele found!!")
else:
	peptide_processing(peptide)
	MHC_processing(MHC)
print(peptide)
print(MHC)

# 2. FEATURIZE INPUT
initialize_dir([filestore + '/aligned_files/', filestore + '/sasa-rsa_files/'])

# 2.1 Find template
template = find_template(peptide, MHC)
print("Template found:")
print(template)

# 2.2 Calculate padding and load the contact list if user asks for it:
if feat_transform == "padding":
	padding_list = MHCFlurry_padding(peptide, len(peptide))
	print("Calculated padding:")
	print(padding_list)
else:
	padding_list = None
if features != "peptide_only":
	with open('./helper_files/total_contact_list.npy', 'rb') as f:
		contact_list = np.load(f)
		contact_list = list(zip(contact_list[:, 0], contact_list[:, 1]))
else:
	contact_list = None

structure_files_to_be_processed = os.listdir(filestore + '/sanitized_input/')
arg_list = []
for structure_file in structure_files_to_be_processed:
	structure_path = filestore + '/sanitized_input/' + structure_file
	arg_list.append((peptide, MHC, structure_path, structure_file, template, filestore, features, feat_transform, feat_terms, padding_list, contact_list))

# 2.3 Extract Rosetta features
with WorkerPool(n_jobs=num_cores) as pool:
	feature_list = pool.map(extract_features, arg_list, progress_bar=True)
features = pd.concat(feature_list)

# Sequential processing for debugging (parallel above)
#for argument in arg_list:
#	print(argument)
#	features = extract_features(argument[0], argument[1], argument[2], argument[3], argument[4], argument[5], argument[6], argument[7], argument[8], argument[9], argument[10])
#	print(features.head(1))
#	print(features.columns)
#	input()

# 2.4 Extract SASA/RSA features (this probably cannot happen in parallel due to race conditions with naccess)
structure_files_to_be_processed = os.listdir(filestore + '/sanitized_input/')
rsa_features_list = []
print("Calculating SASA/RSA...")
for structure_file in tqdm(structure_files_to_be_processed):
	structure_path = filestore + '/sanitized_input/' + structure_file
	index_list = features[features['pdb_code'] == structure_file]['index_list'][0]
	rsa_features_list.append(calculate_sasa(peptide, filestore, structure_path, filestore + '/sasa-rsa_files/', structure_file, index_list, padding_list, sasa_method, feat_transform))
rsa_features = pd.concat(rsa_features_list, axis = 0)
print("SASA/RSA calculation done!")

# 3 RANK INPUT

# 3.1 Prepare input
features = features.merge(rsa_features, on='pdb_code', how='inner').sort_values(by=['pdb_code'])
pdb_codes = features['pdb_code'].to_list()
features = features.drop(['pdb_code', 'filename', 'peptide', 'allele', 'index_list'], axis=1)

# 3.2 Load and rank
pred_list = []
for outer in tqdm(range(1, 7), position=0):
	with open('./models/' + RankMHC_model + '/ranker_' + str(outer) + '.pkl', 'rb') as f:
		model = pkl.load(f)
	pred_list.append(model.predict(features))
preds = np.mean(np.array([pred_list]), axis=1).tolist()[0]
results_csv = pretty_print_analytics(pd.DataFrame(data={'File name':pdb_codes, 'RankMHC Prediction':preds}), filestore)
print("\n\nEnd of RankMHC")
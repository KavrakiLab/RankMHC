import pandas as pd
import pickle as pkl
import numpy as np
from sklearn.model_selection import KFold

cv_file = pd.read_csv("./cv_file.csv")
cv_file['pep_len'] = cv_file['peptide'].str.len()
pdbs = pd.unique(cv_file['pdb_code'])
alleles = pd.unique(cv_file['MHC'])
lengths = pd.unique(cv_file['pep_len'])
kf = KFold(n_splits=6, shuffle=True, random_state=42)
for i, (train_index, test_index) in enumerate(kf.split(lengths)):
	print(f"Fold {i}:")
	print(f"  Train: {lengths[train_index]}")
	print(f"  Test: {lengths[test_index]}")
	inner_lengths = lengths[train_index]
	with open('./pickle_splits/LOLO_Outer_' + str(i + 1) +'.pkl', 'wb') as f:
		pkl.dump(inner_lengths, f)
	kf_inner = KFold(n_splits=5, shuffle=True, random_state=42)
	for j, (train_index_inner, val_index_inner) in enumerate(kf_inner.split(inner_lengths)):
		print(f"Inner Fold {j}:")
		print(f"Inner Train: {inner_lengths[train_index_inner]}")
		print(f"Inner Test: {inner_lengths[val_index_inner]}")
		with open('./pickle_splits/LOLO_Inner_' + str(i + 1) + '_' + str(j + 1) +'.pkl', 'wb') as f:
			pkl.dump(inner_lengths[train_index_inner], f)

for i, (train_index, test_index) in enumerate(kf.split(alleles)):
	print(f"Fold {i}:")
	print(f"  Train: {alleles[train_index]}")
	print(f"  Test: {alleles[test_index]}")
	inner_alleles = alleles[train_index]
	with open('./pickle_splits/LKAO_Outer_' + str(i + 1) +'.pkl', 'wb') as f:
		pkl.dump(inner_alleles, f)
	kf_inner = KFold(n_splits=5, shuffle=True, random_state=42)
	for j, (train_index_inner, val_index_inner) in enumerate(kf_inner.split(inner_alleles)):
		print(f"Inner Fold {j}:")
		print(f"Inner Train: {inner_alleles[train_index_inner]}")
		print(f"Inner Test: {inner_alleles[val_index_inner]}")
		with open('./pickle_splits/LKAO_Inner_' + str(i + 1) + '_' + str(j + 1) +'.pkl', 'wb') as f:
			pkl.dump(inner_alleles[train_index_inner], f)

for i, (train_index, test_index) in enumerate(kf.split(pdbs)):
	print(f"Fold {i}:")
	print(f"  Train: {pdbs[train_index]}")
	print(f"  Test: {pdbs[test_index]}")
	inner_pdbs = pdbs[train_index]
	with open('./pickle_splits/LKPO_Outer_' + str(i + 1) +'.pkl', 'wb') as f:
		pkl.dump(inner_pdbs, f)
	kf_inner = KFold(n_splits=5, shuffle=True, random_state=42)
	for j, (train_index_inner, val_index_inner) in enumerate(kf_inner.split(inner_pdbs)):
		print(f"Inner Fold {j}:")
		print(f"Inner Train: {inner_pdbs[train_index_inner]}")
		print(f"Inner Test: {inner_pdbs[val_index_inner]}")
		with open('./pickle_splits/LKPO_Inner_' + str(i + 1) + '_' + str(j + 1) +'.pkl', 'wb') as f:
			pkl.dump(inner_pdbs[train_index_inner], f)	
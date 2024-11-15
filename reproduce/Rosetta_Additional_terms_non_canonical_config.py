import pymol2
from Bio.PDB import PDBParser
from Bio import Align
from Bio import pairwise2
import pyrosetta
from pyrosetta import *

import os
from tqdm import tqdm
from copy import deepcopy
from itertools import zip_longest
import argparse

import pandas as pd
import numpy as np 

import seaborn as sns
import matplotlib.pyplot as plt

# util functions for extrating the data
def init_rosetta():
	pyrosetta.init("-mute all")

def create_scoring_function():

	cen_sfxn = ScoreFunction()

	# Standard ref2015 terms
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.fa_atr, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.fa_rep, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.fa_sol, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.fa_intra_sol_xover4, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.lk_ball_wtd, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.fa_intra_rep, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.fa_elec, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.pro_close, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.hbond_sr_bb, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.hbond_lr_bb, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.hbond_bb_sc, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.dslf_fa13, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.rama_prepro, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.omega, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.p_aa_pp, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.fa_dun, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.yhh_planarity, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.fa_intra_atr, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.fa_intra_sol, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.lk_ball, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.lk_ball_iso, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.hbond, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.mm_lj_intra_rep, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.mm_lj_intra_atr, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.mm_lj_inter_rep, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.mm_lj_inter_atr, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.mm_twist, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.mm_bend, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.mm_stretch, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.lk_costheta, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.lk_polar, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.lk_nonpolar, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.fa_elec_bb_bb, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.fa_elec_bb_sc, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.fa_elec_sc_sc, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.peptide_bond, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.pcs2, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.cen_pair_motifs, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.fa_elec_aro_aro, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.fa_elec_aro_all, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.ch_bond, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.ch_bond_bb_bb, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.ch_bond_sc_sc, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.ch_bond_bb_sc, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.rama2b, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.cenpack, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.cenpack_smooth, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.gauss, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.fa_pair_aro_aro, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.fa_pair_aro_pol, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.fa_pair_pol_pol, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.occ_sol_exact, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.pair, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.cen_pair_smooth, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.interchain_pair, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.interchain_vdw, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.gb_elec, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.envsmooth, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.PB_elec, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.cen_env_smooth, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.cbeta_smooth, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.env, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.cbeta, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.DFIRE, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.rg, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.sa, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.cart_bonded, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.cart_bonded_angle, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.cart_bonded_length, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.neigh_vect, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.neigh_count, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.neigh_vect_raw, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.interface_dd_pair, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.goap, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.goap_dist, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.goap_angle, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.pack_stat, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.surface, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.hpatch, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.fa_dun_dev, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.fa_dun_rot, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.fa_dun_semi, 1.0)
	cen_sfxn.set_weight(pyrosetta.rosetta.core.scoring.co, 1.0)

	return cen_sfxn

def extract_fullc_energies(pose, score_func):
	#load rosett
	score_func(pose)
	#get energies
	res_ene = pose.energies().total_energies_array()
	return res_ene

def extract_ppp_energies(pose, pep_len, score_func, amino_acids):
	#load rosetta
	score_func(pose)
	#get energies
	res_ene = pose.energies().residue_total_energies_array()
	peptide_ene = res_ene[amino_acids]
	return peptide_ene

def extract_rsa_feats(pdb_code, filename, amino_acids):
	rsa = pd.read_csv("./RSA_csvs_ensemble/" + str(pdb_code) + "/" + filename, header = None)[[3, 4, 5, 6]]
	return np.array(rsa.iloc[amino_acids, :].mean())

def align(template, result, result_file_name):

	p1 = pymol2.PyMOL()
	p1.start()
	p1.cmd.load(result, "mobile")
	p1.cmd.load(template, "ref")

	p1.cmd.align("mobile & chain A", "ref & chain A")

	p1.cmd.save(result_file_name, "mobile")
	#p1.cmd.save(template_file_name, "ref")
	p1.stop()

def anchor_alignment(seq1, seq2, anchor_diff1, anchor_diff2):
	
	# Actual alignment part
	extra_in_beginning = ''.join('-'*abs(anchor_diff1))
	extra_in_end = ''.join('-'*abs(anchor_diff2))
	if anchor_diff1 > 0:
		temp_seq1 = extra_in_beginning + seq1
		temp_seq2 = seq2
	else:
		temp_seq1 = seq1
		temp_seq2 = extra_in_beginning + seq2
	if anchor_diff2 > 0:
		temp_seq1 = temp_seq1 + extra_in_end
	else:
		temp_seq2 = temp_seq2 + extra_in_end

	# Add to anything:
	diff = len(temp_seq1) - len(temp_seq2)
	if diff >= 0:
		temp_seq2 = temp_seq2 + ''.join('-'*diff)
	else:
		temp_seq1 = temp_seq1 + ''.join('-'*abs(diff))

	# Remove redundancies
	[temp_seq1, temp_seq2] = [''.join(k) for k in zip(*[i for i in zip_longest(temp_seq1, temp_seq2, fillvalue = "") if (i[0] != '-' or i[1] != '-')])]

	return temp_seq1, temp_seq2

def score_sequences(seq1, seq2, anchor_1, anchor_2, matrix, gap_penalty, norm):

	# Initial checks
	anchor_1_idx = anchor_1 + (len(seq2) - len(seq2.lstrip('-'))) - 1
	anchor_2_idx = anchor_2 + (len(seq2) - len(seq2.lstrip('-'))) - 1
	if (seq2[anchor_1_idx] != '-') and (seq1[anchor_1_idx] == '-'): # When N-anchor is not covered
		return -1000
	if (seq2[anchor_2_idx] != '-') and (seq1[anchor_2_idx] == '-'): # When C-anchor is not covered
		return -1000
	if (anchor_1 >= 3) and ((len(seq2) - len(seq2.lstrip('-'))) > 0): # When N-anchor is found at pos 3, then we cannot add aas before
		return -1000
	if (anchor_1 == 2) and ((len(seq2) - len(seq2.lstrip('-'))) > 1): # When N-anchor is found at pos 2, then we cannot add more than 1 aa before it.
		return -1000
	if (anchor_1 == 1) and ((len(seq2) - len(seq2.lstrip('-'))) > 2): # When N-anchor is found at pos 1, then we cannot add more than 2 aa before it.
		return -1000

	# Actual scoring
	score = 0
	num_gaps = 0
	for A,B in zip(seq1,seq2):
		diag = ('-'==A) or ('-'==B)
		score += 0 if diag else matrix[A,B]
		num_gaps += 1 if diag else 0
	score = score/norm
	return score + num_gaps*gap_penalty

def MHCscore_feature_extraction_parser():
	parser = argparse.ArgumentParser(description="Feature Extraction Parser", 
									 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('--pdb_code', type=str, nargs=1, help='The PDB code that will be used')
	return parser

# Parse args
parser = MHCscore_feature_extraction_parser()
args = parser.parse_args()
if not (args.pdb_code):
    parser.error('No action requested, add --pdb_code')
pdb_code = args.pdb_code[0]

# Initialize Rosetta
init_rosetta()

# Bring all needed arguments
template_db = pd.read_csv("./Template_DB.csv")[['pdb_code', 'MHC', 'peptide', 'peptide_length', 'Major_anchors']]
template_db = template_db[template_db['pdb_code'] == pdb_code]
MHC_of_interest = template_db['MHC'].values[0]
peptide = template_db['peptide'].values[0]
peptide_anchors = template_db['Major_anchors'].values[0]
	
print(pdb_code)
print(MHC_of_interest)
print(peptide)
print(peptide_anchors)
input()

#Pre-process anchors first
peptide_anchors = [int(pos) for pos in peptide_anchors.split(",")]
blosum_62 = Align.substitution_matrices.load("BLOSUM62")

template_db = pd.read_csv("./Template_DB.csv")[['pdb_code', 'MHC', 'peptide', 'peptide_length', 'Major_anchors']]
template_db = template_db[template_db['peptide_length'] == 9]
template_db = template_db[template_db['Major_anchors'] == "2,9"]
template_db['Major_anchor_1'] = template_db['Major_anchors'].apply(lambda x: x.split(",")).str[0].astype(int)
template_db['Major_anchor_2'] = template_db['Major_anchors'].apply(lambda x: x.split(",")).str[1].astype(int)
template_db['Anchor_diff_1'] = template_db['Major_anchor_1'] - peptide_anchors[0]
template_db['Anchor_diff_2'] = template_db['Major_anchor_2'] - template_db['peptide_length'] + peptide_anchors[1] - len(peptide)
	
# Peptide similarity
template_sequences = template_db['peptide'].tolist()
Anchor_diff_1 = template_db['Anchor_diff_1'].tolist()
Anchor_diff_2 = template_db['Anchor_diff_2'].tolist()
anchor_1_list = template_db['Major_anchor_1'].tolist()
anchor_2_list = template_db['Major_anchor_2'].tolist()
self_score = score_sequences(peptide, peptide, 0, 0, matrix=blosum_62, gap_penalty=0, norm=1)
score_list = []
for i, template_sequence in enumerate(template_sequences):
	temp_sequence_in_question, temp_template_sequence = anchor_alignment(peptide, template_sequence, 
																		 Anchor_diff_1[i], Anchor_diff_2[i])
	score_list.append(score_sequences(temp_sequence_in_question, temp_template_sequence, anchor_1_list[i], 
									  anchor_2_list[i], matrix=blosum_62, gap_penalty=0, norm=self_score))
template_db['Peptide_similarity'] = score_list

# MHC similarity
similarity_matrix = pd.read_csv("./9mer_similarity.csv")
similarity_matrix = similarity_matrix[similarity_matrix['Allele'] == MHC_of_interest].T.reset_index().iloc[1:]
similarity_matrix.columns = ['MHC', 'MHC_similarity']
template_db = template_db.merge(similarity_matrix, on=['MHC'])
template_db = template_db[template_db['pdb_code'] != pdb_code]
	
# Final Choice
template_db['Similarity_score'] = 0.5*template_db['MHC_similarity'] + 0.5*template_db['Peptide_similarity']
print(template_db)
input()
template = template_db[template_db['Similarity_score'] == template_db['Similarity_score'].max()].sample(n=1)['pdb_code'].values[0]
print("Template chosen:")
print(template)

# Structural Alignment
template_file_name = './pMHCI_AC/' + template
result_file_name = './leftover_files/' + pdb_code

file_list = [x for x in os.listdir('./Structures/' + pdb_code + '/results/6_final_conformations/')]
file_list.sort()
peptide_only_df_list = []
df_list = []
for filename in tqdm(file_list):

	align(template_file_name, 
		  './Structures/' + pdb_code + '/results/6_final_conformations/' + filename, 
		  result_file_name)

	parser = PDBParser()
	template_structure = parser.get_structure('template', template_file_name)
	template_model = template_structure[0]
	template_chain = template_model['C']
	result_structure = parser.get_structure('result', result_file_name)
	result_model = result_structure[0]
	result_chain = result_model['C']
	a = ["C" + residue.get_resname() + str(residue.id[1]) for residue in result_model['A'].get_residues()]
	c = ["C" + residue.get_resname() + str(residue.id[1]) for residue in result_chain.get_residues()]
	residue_names = np.array(a + c)

	index_list = []
	for i, residue1 in enumerate(result_chain):
		min_dist = float("inf")
		min_j = 0
		for j, residue2 in enumerate(template_chain):
			try:
				if (residue1.get_resname() == "GLY") or (residue2.get_resname() == "GLY"): 
					dist = abs(residue1['CA'] - residue2['CA'])
				else:
					dist = abs(residue1['CB'] - residue2['CB'])
			except KeyError:
				## no CA atom, e.g. for H_NAG
				dist = float("inf")
			if dist < min_dist:
				min_dist = dist
				min_j = j + 1
		index_list.append(min_j)
	index_list = np.array(index_list)
	print("Calculated indexes:")
	print(index_list)

	# Global features
	pose = pyrosetta.pose_from_pdb(result_file_name)
	score_func = create_scoring_function()
	res_ene = extract_fullc_energies(pose, score_func)
	res_ene_names = res_ene.dtype.names
	res_dict = {'pdb_code': [pdb_code],
		        'filename': [filename],
		        'allele': [MHC_of_interest], 
	            'peptide': [peptide]
		   	   }
	for i, x in enumerate(res_ene[0]):
		if res_ene_names[i] not in res_dict: res_dict[res_ene_names[i]] = []
		res_dict[res_ene_names[i]].append(res_ene[0][i])
		
	# Per-peptide (ppp) features
	feat_dict = {}
	for i in range(1, 10):
		amino_acids = np.argwhere(index_list == i).tolist()	
		amino_acids = [amino_acid_2 for amino_acid in amino_acids for amino_acid_2 in amino_acid]
		amino_acids = np.array([pose.total_residue() - len(peptide) + amino_acid for amino_acid in amino_acids])
		if not amino_acids:
			peptide_ene = np.zeros(85)
		else:		
			peptide_ene = extract_ppp_energies(pose, len(peptide), score_func, amino_acids)
			peptide_ene = peptide_ene.view((float, len(peptide_ene.dtype.names)))
			peptide_ene = np.mean(peptide_ene, axis=0)
		feat_dict[i] = peptide_ene
	peptide_feat = np.array([feat_dict[i] for i in range(1, 10)])
	peptide_feature_names = ['{}_{}'.format(feat, i) for i in range(1, 10) for feat in list(res_dict.keys())[4:]]
	peptide_feat = pd.DataFrame(peptide_feat.flatten()).T
	peptide_feat.columns = peptide_feature_names

	# RSA feats
	feat_dict = {}
	for i in range(1, 10):
		amino_acids = np.argwhere(index_list == i).tolist()	
		amino_acids = [amino_acid_2 for amino_acid in amino_acids for amino_acid_2 in amino_acid]
		if not amino_acids:
			peptide_rsa = np.zeros(4)
		else:
			peptide_rsa = extract_rsa_feats(pdb_code, filename, amino_acids)
		feat_dict[i] = peptide_rsa
	rsa_feat = np.array([feat_dict[i] for i in range(1, 10)])
	rsa_feature_names = ['{}_{}'.format(feat, i) for i in range(1, 10) for feat in ['SASA_all', 'RSA_all',
																					'SASA_side', 'RSA_side']]
	rsa_feat = pd.DataFrame(rsa_feat.flatten()).T
	rsa_feat.columns = rsa_feature_names

	# HLA feats
	HLA_ene = extract_ppp_energies(pose, len(peptide), score_func, list(range(0, 180)))
	HLA_ene = HLA_ene.view((float, len(HLA_ene.dtype.names)))
	HLA_feature_names = ['{}_{}'.format(feat, i) for i in range(1, 181) for feat in list(res_dict.keys())[4:]]
	HLA_feat = pd.DataFrame(HLA_ene.flatten()).T
	HLA_feat.columns = HLA_feature_names

	# Construct final DF
	peptide_only_df_list.append(pd.concat([pd.DataFrame(res_dict), peptide_feat, rsa_feat], axis = 1))
	df_list.append(pd.concat([pd.DataFrame(res_dict), peptide_feat, rsa_feat, HLA_feat], axis = 1))

# Construct final DF
peptide_only_final_df = pd.concat(peptide_only_df_list, axis = 0)
final_df = pd.concat(df_list, axis = 0)
final_df.to_csv('./Rosetta_features_big_non_canonical_peptide/' + pdb_code + ".csv")
final_df.to_csv('./Rosetta_features_big_non_canonical_HLA/' + pdb_code + ".csv")
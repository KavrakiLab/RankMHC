import pymol2
from Bio.PDB import PDBParser
import pyrosetta
from pyrosetta import *

import pandas as pd
import numpy as np 
import math

from subprocess import call
import sys

from ortools.linear_solver import pywraplp

from feature_extraction.scoring_terms_definition import create_scoring_function, extract_pairwise_scores

def init_rosetta():
	pyrosetta.init("-mute all")

def MHCFlurry_padding(peptide, pep_len, left_edge = 4, right_edge = 4, min_length = 8, max_length = 15):

	# Routine based on MHCFlurry Neutral Amino Acid 'X' padding
	middle_length = max_length - left_edge - right_edge
	num_null = max_length - pep_len
	num_null_left = int(math.ceil(num_null/2))
	num_null_right = num_null - num_null_left
	num_middle_filled = middle_length - num_null
	padding_lists = [[1,2,3,4], [-1]*num_null_left, list(range(5, 5 + num_middle_filled)), 
					 [-1]*num_null_right, list(range(5 + num_middle_filled, pep_len + 1))]

	return [elem for pad_list in padding_lists for elem in pad_list]

def extract_energies(peptide, res_dict, pose, score_func, index_list, eightmer_assignment, feature_selection, feat_transform, feat_terms, padding_list, contact_list, feat_dim, routine_pass):

	# Score using pyrosetta
	score_func(pose)

	# Global features 
	res_ene = pose.energies().total_energies_array()
	res_ene_names = res_ene.dtype.names
	for i, x in enumerate(res_ene[0]):
		if res_ene_names[i] not in res_dict: res_dict[res_ene_names[i]] = []
		res_dict[res_ene_names[i]].append(res_ene[0][i])
	global_feat = pd.DataFrame(res_dict)
	global_feat.columns = ['global_' + col if col not in ['pdb_code', 'filename', 'allele', 'peptide'] else col for col in global_feat.columns]

	
	# Per-amino acid features based on indexes
	aa_ene = pose.energies().residue_total_energies_array() # Per-amino acid features
	aa_ene_names = aa_ene.dtype.names

	if (feature_selection != 'pairwise') and (feat_transform == 'pooling'):
		opt_dict = {}
		for i in range(1, 10):
			amino_acids = np.argwhere(index_list == i).tolist()	
			amino_acids = [amino_acid_2 for amino_acid in amino_acids for amino_acid_2 in amino_acid]
			amino_acids = np.array([pose.total_residue() - len(peptide) + amino_acid for amino_acid in amino_acids])
			if not amino_acids.size:
				peptide_ene = aa_ene[np.array([pose.total_residue() - len(peptide) + eightmer_assignment - 1])]
				peptide_ene = peptide_ene.view((float, len(peptide_ene.dtype.names)))
				peptide_ene = np.mean(peptide_ene, axis=0)
				peptide_ene = np.zeros(feat_dim)
				opt_dict[i] = peptide_ene
			else:
				peptide_ene = aa_ene[amino_acids]
				peptide_ene = peptide_ene.view((float, len(peptide_ene.dtype.names)))
				peptide_ene = np.mean(peptide_ene, axis=0)
				opt_dict[i] = peptide_ene
		opt_feat = np.array([opt_dict[i] for i in range(1, 10)])
		opt_feature_names = ['pep_{}_{}'.format(feat, i) for i in range(1, 10) for feat in list(res_dict.keys())[4:]]
		opt_feat = pd.DataFrame(opt_feat.flatten()).T
		opt_feat.columns = opt_feature_names

	# Padding option
	if (feature_selection != 'pairwise') and (feat_transform == 'padding'):
		pad_dict = {}
		for i in range(1, 16):
			if padding_list[i - 1] == -1:
				peptide_ene = np.zeros(feat_dim)
			else:
				amino_acids = np.array([pose.total_residue() - len(peptide) + padding_list[i - 1] - 1])
				peptide_ene = aa_ene[amino_acids]
				peptide_ene = peptide_ene.view((float, len(peptide_ene.dtype.names)))
				peptide_ene = np.mean(peptide_ene, axis=0)
			pad_dict[i] = peptide_ene
		pad_feat = np.array([pad_dict[i] for i in range(1, 16)])
		pad_feature_names = ['pad_{}_{}'.format(feat, i) for i in range(1, 16) for feat in list(res_dict.keys())[4:]]
		pad_feat = pd.DataFrame(pad_feat.flatten()).T
		pad_feat.columns = pad_feature_names

	# HLA feats option
	if (feature_selection == 'peptide_HLA'): 
		hla_res_list = list(set([contact[0] for contact in contact_list]))
		hla_res_list = [hla_res - 1 for hla_res in hla_res_list]
		HLA_ene = aa_ene[hla_res_list]
		HLA_ene = HLA_ene.view((float, len(HLA_ene.dtype.names)))
		HLA_feature_names = ['HLA_{}_{}'.format(feat, i + 1) for i in hla_res_list for feat in list(res_dict.keys())[4:]]
		HLA_feat = pd.DataFrame(HLA_ene.flatten()).T
		HLA_feat.columns = HLA_feature_names

	# Energy graph contact features (pairwise option)
	if feature_selection == 'pairwise':

		pairwise_names = ['fa_atr', 'fa_rep', 'fa_sol', 'fa_intra_sol_xover4', 'lk_ball_wtd', 'fa_intra_rep', 'fa_elec',
				  'pro_close', 'hbond_sr_bb', 'hbond_lr_bb', 'hbond_bb_sc', 'hbond_sc', 'dslf_fa13', 'rama_prepro', 'omega',
				  'p_aa_pp', 'fa_dun', 'yhh_planarity', 'ref', 'fa_intra_atr', 'fa_intra_sol', 'lk_ball', 'lk_ball_iso',
				  'hbond', 'mm_lj_intra_rep', 'mm_lj_intra_atr', 'mm_lj_inter_rep', 'mm_lj_inter_atr', 'mm_twist',
				  'mm_bend', 'mm_stretch', 'lk_costheta', 'lk_polar', 'lk_nonpolar', 'fa_elec_bb_bb',
				  'fa_elec_bb_sc', 'fa_elec_sc_sc', 'peptide_bond', 'pcs2', 'cen_pair_motifs', 'fa_elec_aro_aro',
				  'fa_elec_aro_all', 'ch_bond', 'ch_bond_bb_bb', 'ch_bond_sc_sc', 'ch_bond_bb_sc', 'rama2b',
				  'cenpack', 'cenpack_smooth', 'gauss', 'fa_pair_aro_aro', 'fa_pair_aro_pol', 'fa_pair_pol_pol',
				  'occ_sol_exact', 'pair', 'cen_pair_smooth', 'interchain_pair', 'interchain_vdw', 'gb_elec', 
				  'envsmooth', 'cen_env_smooth', 'cbeta_smooth', 'env', 'cbeta', 'DFIRE', 'rg', 'sa', 'cart_bonded',
				  'cart_bonded_angle', 'cart_bonded_length', 'neigh_vect', 'neigh_count', 'neigh_vect_raw', 
				  'interface_dd_pair', 'goap', 'goap_dist', 'goap_angle', 'pack_stat', 'surface', 'hpatch', 
				  'fa_dun_dev', 'fa_dun_rot', 'fa_dun_semi', 'co']

		if feat_dim > 20:
			energy_dict = {}
			energy_graph = pose.energies().energy_graph()
			for contact in contact_list:
				amino_acids = np.argwhere(index_list == contact[1]).tolist()	
				amino_acids = [amino_acid_2 + 1 for amino_acid in amino_acids for amino_acid_2 in amino_acid]
				amino_acids = np.array([pose.total_residue() - len(peptide) + amino_acid for amino_acid in amino_acids])
				if not amino_acids.size:
					try:
						peptide_ene = extract_pairwise_scores(energy_graph.find_edge(contact[0], 
													  	  	  pose.total_residue() - len(peptide) + eightmer_assignment).fill_energy_map())
					except AttributeError:
						peptide_ene = np.zeros(feat_dim - 2)
					peptide_ene = peptide_ene.reshape(1, -1)
					peptide_ene = np.mean(peptide_ene, axis=0)
					peptide_ene = np.zeros(feat_dim - 2)
					energy_dict[(contact[0], contact[1])] = peptide_ene
				else:
					peptide_ene = np.empty((1,feat_dim - 2))
					for amino_acid in amino_acids:
						try:
							e_feats = extract_pairwise_scores(energy_graph.find_edge(contact[0], amino_acid).fill_energy_map())
							e_feats = e_feats.reshape(1, -1)
						except AttributeError:
							e_feats = np.zeros(feat_dim - 2)
							e_feats = e_feats.reshape(1, -1)
						peptide_ene = np.append(peptide_ene, e_feats, axis = 0)
					peptide_ene = np.mean(peptide_ene[1:], axis=0)
					energy_dict[(contact[0], contact[1])] = peptide_ene
			energy_feat = np.array([energy_dict[(contact[0], contact[1])] for contact in contact_list])
			energy_feature_names = ['pair_{}_{}'.format(feat, (str(contact[0]) + '-' + str(contact[1]))) for contact in contact_list for feat in pairwise_names]
			energy_feat = pd.DataFrame(energy_feat.flatten()).T
			energy_feat.columns = energy_feature_names
		else:
			energy_feat = None

	# Concat and clean features
	feat_name_dir = './feature_extraction/feature_names/' + feature_selection + '_' + feat_transform + '_' + feat_terms + '_names.txt'
	if (feature_selection == 'peptide_only') and (feat_transform == 'padding'):
		feat = pd.concat([global_feat, pad_feat], axis = 1)
	if (feature_selection == 'peptide_only') and (feat_transform == 'pooling'):
		feat = pd.concat([global_feat, opt_feat], axis = 1)
	if (feature_selection == 'peptide_HLA') and (feat_transform == 'pooling'):
		feat = pd.concat([global_feat, opt_feat, HLA_feat], axis = 1)
	if (feature_selection == 'pairwise'):
		feat = pd.concat([global_feat, energy_feat], axis = 1)
	with open(feat_name_dir, 'r') as f:
		feat_names = f.read().replace('\n', '').split(",")
	feat = feat.loc[:, feat.columns.isin(feat_names)]

	# Final cleaning depending on the routine called
	if (routine_pass == "first"):
		feat = pd.concat([pd.DataFrame(res_dict)[['pdb_code', 'filename', 'allele', 'peptide']], feat.filter(regex='total_score')], axis = 1)
	else:
		feat = feat[[col for col in feat.columns if "global_total_score" not in col]]
		feat = feat[[col for col in feat.columns if "pep_total_score_" not in col]]
		feat = feat[[col for col in feat.columns if "pad_total_score_" not in col]]
		feat = feat[[col for col in feat.columns if "HLA_total_score_" not in col]]
		feat = feat[[col for col in feat.columns if "pair_total_score_" not in col]]
		feat = feat.drop(['pdb_code', 'filename', 'allele', 'peptide'], axis=1)
	return feat

def extract_rsa_feats(pdb_file, amino_acids, rsa_path):
	rsa = pd.read_csv(rsa_path + pdb_file.split(".")[0] + '_rsa_sasa.csv', header = None)[[3, 4, 5, 6]]
	try:
		rsa = np.array(rsa.iloc[amino_acids, :].mean())
		return rsa
	except IndexError:
		print(pdb_file)
		print(amino_acids)
		print(rsa_path)
		print(rsa)
		return np.array(rsa.iloc[amino_acids, :].mean())

def calculate_sasa(peptide, filestore, pdb_file_absolute_path, pdb_file_wanted_path, pdb_file, index_list, padding_list, sasa_method, feat_transform):

	# Call naccess/freesasa
	if sasa_method == "naccess":
		call(["bash ./feature_extraction/run_naccess.sh "  + pdb_file_absolute_path + " " + pdb_file_wanted_path + " " + pdb_file.split(".")[0] + " " + peptide + " > " + filestore + "/sasa-rsa_files/logs.log 2>&1"], shell=True)
	else:
		call(["bash ./feature_extraction/run_freesasa.sh " + pdb_file_absolute_path + " " + pdb_file_wanted_path + " " + pdb_file.split(".")[0] + " " + peptide + " > " + filestore + "/sasa-rsa_files/logs.log 2>&1"], shell=True)

	# Per-amino acid feature sasa/rsa
	rsa_dict = {}
	if (feat_transform == 'pooling'):
		for i in range(1, 10):
			amino_acids = np.argwhere(index_list == i).tolist()	
			amino_acids = [amino_acid_2 for amino_acid in amino_acids for amino_acid_2 in amino_acid]
			if not amino_acids:
				peptide_rsa = np.zeros(4)
				rsa_dict[i] = peptide_rsa
			else:
				peptide_rsa = extract_rsa_feats(pdb_file, amino_acids, filestore + '/sasa-rsa_files/')
				rsa_dict[i] = peptide_rsa
		rsa_feat = np.array([rsa_dict[i] for i in range(1, 10)])
		rsa_feature_names = ['rsa_{}_{}'.format(feat, i) for i in range(1, 10) for feat in ['SASA_all', 'RSA_all', 'SASA_side', 'RSA_side']]
		rsa_feat = pd.DataFrame(rsa_feat.flatten()).T
		rsa_feat.columns = rsa_feature_names
	else:
		for i in range(1, 16):
			if padding_list[i - 1] == -1:
				peptide_rsa = np.zeros(4)
			else:
				amino_acids = np.array([padding_list[i - 1] - 1])
				peptide_rsa = extract_rsa_feats(pdb_file, amino_acids, filestore + '/sasa-rsa_files/')
			rsa_dict[i] = peptide_rsa
		rsa_feat = np.array([rsa_dict[i] for i in range(1, 16)])
		rsa_feature_names = ['radsa_{}_{}'.format(feat, i) for i in range(1, 16) for feat in ['SASA_all', 'RSA_all', 'SASA_side', 'RSA_side']]
		rsa_feat = pd.DataFrame(rsa_feat.flatten()).T
		rsa_feat.columns = rsa_feature_names
	rsa_feat['pdb_code'] = pdb_file
	return rsa_feat

def align(template, result, result_file_name):

	p1 = pymol2.PyMOL()
	p1.start()
	p1.cmd.load(result, "mobile")
	p1.cmd.load(template, "ref")

	p1.cmd.align("mobile & chain A", "ref & chain A")

	p1.cmd.save(result_file_name, "mobile")
	#p1.cmd.save(template_file_name, "ref")
	p1.stop()

def index_assignment(result_chain, template_chain, pep_len):

	# Calculate all-vs-all residue distances (costs)
	costs = []
	for i, residue1 in enumerate(result_chain):
		per_residue_costs = []
		for j, residue2 in enumerate(template_chain):
			try:
				if (residue1.get_resname() == "GLY") or (residue2.get_resname() == "GLY"): 
					dist = np.sqrt(residue1['CA'] - residue2['CA'])
				else:
					dist = np.sqrt(residue1['CB'] - residue2['CB'])
			except KeyError:
				## no CA atom, e.g. for H_NAG
				dist = float("inf")
			per_residue_costs.append(dist)
		costs.append(per_residue_costs)

	num_workers = len(costs)
	num_tasks = len(costs[0])

	solver = pywraplp.Solver.CreateSolver("SCIP")
	if not solver:
		print("Error in solver")
		sys.exit(1)

	x = {}
	for i in range(num_workers):
		for j in range(num_tasks):
			x[i, j] = solver.IntVar(0, 1, "")

	# Each amino acid is assigned to exactly one bin.
	if pep_len > 8:
		for i in range(num_workers):
			solver.Add(solver.Sum([x[i, j] for j in range(num_tasks)]) == 1)
	else:
		for i in range(num_workers):
			solver.Add(solver.Sum([x[i, j] for j in range(num_tasks)]) >= 1)

	# Each bin is assigned to at least one amino acid (to avoid zeroes and padding).
	for j in range(num_tasks):
		solver.Add(solver.Sum([x[i, j] for i in range(num_workers)]) >= 1)

	objective_terms = []
	for i in range(num_workers):
		for j in range(num_tasks):
			objective_terms.append(costs[i][j] * x[i, j])
	solver.Minimize(solver.Sum(objective_terms))

	status = solver.Solve()

	index_list = []
	if status == pywraplp.Solver.OPTIMAL or status == pywraplp.Solver.FEASIBLE:
		if pep_len > 8:
			for i in range(num_workers):
				for j in range(num_tasks):
					# Test if x[i,j] is 1 (with tolerance for floating point arithmetic).
					if x[i, j].solution_value() > 0.5:
						index_list.append(j + 1)
			eightmer_assignment = None
		else:
			previous_worker = -1
			for i in range(num_workers):
				for j in range(num_tasks):
					# Test if x[i,j] is 1 (with tolerance for floating point arithmetic).
					if x[i, j].solution_value() > 0.5:
						if previous_worker == i:
							eightmer_assignment = i + 1
						else:
							index_list.append(j + 1)
							previous_worker = i
	else:
		print("No COP solution found.")
		sys.exit(1)
	return np.array(index_list), eightmer_assignment

def extract_features(peptide, MHC_of_interest, pdb_file_path, pdb_file, template, filestore, feature_selection, feat_transform, feat_terms, padding_list, contact_list):

	template_file_name = './templates/' + template
	result_file_name = filestore + '/aligned_files/' + pdb_file

	# 1. Align the structure with the nonamer template
	align(template_file_name, 
		  pdb_file_path, 
		  result_file_name)

	# 2. Solve the assignment problem and find the optimal matching
	parser = PDBParser(QUIET=True)
	template_structure = parser.get_structure('template', template_file_name)
	template_model = template_structure[0]
	template_chain = template_model['C']
	result_structure = parser.get_structure('result', result_file_name)
	result_model = result_structure[0]
	result_chain = result_model['C']
	a = ["C" + residue.get_resname() + str(residue.id[1]) for residue in result_model['A'].get_residues()]
	c = ["C" + residue.get_resname() + str(residue.id[1]) for residue in result_chain.get_residues()]
	residue_names = np.array(a + c)

	index_list, eightmer_assignment = index_assignment(result_chain, template_chain, len(peptide))

	# Re-specify the file
	result_file_name = filestore + '/sanitized_input/' + pdb_file

	# Initialize Rosetta
	init_rosetta()

	# Extract ref2015 total features
	res_dict = {'pdb_code': [pdb_file],
				'filename': [pdb_file_path],
				'allele': [MHC_of_interest], 
				'peptide': [peptide]
		   	   }
	pose = pyrosetta.pose_from_pdb(result_file_name)
	scorefxn = pyrosetta.get_fa_scorefxn()
	ref2015_feat = extract_energies(peptide, res_dict, pose, scorefxn, index_list, eightmer_assignment, feature_selection, feat_transform, feat_terms, padding_list, contact_list, feat_dim = 20, routine_pass = "first")

	# Extract custom additional features
	res_dict = {'pdb_code': [pdb_file],
				'filename': [pdb_file_path],
				'allele': [MHC_of_interest], 
				'peptide': [peptide]
		   	   }
	score_func = create_scoring_function()
	if feat_terms == "full":
		custom_feat = extract_energies(peptide, res_dict, pose, score_func, index_list, eightmer_assignment, feature_selection, feat_transform, feat_terms, padding_list, contact_list, feat_dim = 86, routine_pass = "second")
	else:
		custom_feat = extract_energies(peptide, res_dict, pose, scorefxn, index_list, eightmer_assignment, feature_selection, feat_transform, feat_terms, padding_list, contact_list, feat_dim = 20, routine_pass = "second")

	# Extract rsa features
	#rsa_feat = extract_rsa(peptide, filestore, result_file_name, filestore + '/sasa-rsa_files/', pdb_file, index_list, padding_list, sasa_method, feat_transform)

	# Construct final DFs (all versions)
	#final_feat = pd.concat([ref2015_feat, custom_feat, rsa_feat], axis = 1)
	final_feat = pd.concat([ref2015_feat, custom_feat], axis = 1)
	final_feat['index_list'] = [index_list]
	return final_feat
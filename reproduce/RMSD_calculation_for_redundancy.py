import sys
import os
from mpire import WorkerPool

import pymol2
from biopandas.pdb import PandasPdb


import pandas as pd
import numpy as np
import random
import hdbscan
from sklearn.metrics.pairwise import pairwise_distances

def align(template, result, result_file_name):

	p1 = pymol2.PyMOL()
	p1.start()
	p1.cmd.load(result, "mobile")
	p1.cmd.load(template, "ref")

	p1.cmd.align("mobile & chain A", "ref & chain A")

	p1.cmd.save(result_file_name, "mobile")
	#p1.cmd.save(template_file_name, "ref")
	p1.stop()

def calculate_rmsd(template_file_name, result_file_name, pep_len):

	# Load template file
	template_peptide = PandasPdb()
	template_peptide.read_pdb(template_file_name)
	template_peptide.df['ATOM'] = template_peptide.df['ATOM'][template_peptide.df['ATOM']['chain_id'] == 'C']

	peptide_result = PandasPdb()
	peptide_result.read_pdb(result_file_name)
	peptide_result.df['ATOM'] = peptide_result.df['ATOM'][peptide_result.df['ATOM']['chain_id'] == 'C']

	peptide_result.df['ATOM'] = peptide_result.df['ATOM'].merge(template_peptide.df['ATOM'][['atom_name','residue_number']], on=['atom_name','residue_number'])
	template_peptide.df['ATOM'] = template_peptide.df['ATOM'].merge(peptide_result.df['ATOM'][['atom_name','residue_number']], on=['atom_name','residue_number'])

	max_rmsd = 0
	for pos in range(1, pep_len + 1):

		temp_template = template_peptide.df['ATOM'][template_peptide.df['ATOM']['residue_number'] == pos]
		temp_template = temp_template.set_index('atom_name')

		res = peptide_result.df['ATOM'][peptide_result.df['ATOM']['residue_number'] == pos]
		temp_template = temp_template.reindex(index=res['atom_name'])
		temp_template = temp_template.reset_index()

		rmsd = PandasPdb.rmsd(temp_template, res, s='hydrogen', invert = True) # all atoms, excluding hydrogens
		if rmsd > max_rmsd:
			max_rmsd = rmsd

	return rmsd

def rmsd_calc(pdb_code, pep_len):

	template_file_name = './pMHCI_AC/' + pdb_code 
	result_file_name = './temp_files/' + pdb_code

	file_list = [x for x in os.listdir('./Structures/' + pdb_code + '/results/6_final_conformations/')]
	file_list.sort() 
	
	rmsd_dict = {}
	for filename in file_list:
		rmsd_dict[filename] = {}
		for filename2 in file_list:

			align('./Structures/' + pdb_code + '/results/6_final_conformations/' + filename, 
			  	  './Structures/' + pdb_code + '/results/6_final_conformations/' + filename2, 
			  	  result_file_name)

			rmsd = calculate_rmsd('./Structures/' + pdb_code + '/results/6_final_conformations/' + filename, 
								  result_file_name, pep_len)
			rmsd_dict[filename][filename2] = rmsd

	rmsd_df = pd.DataFrame.from_dict(rmsd_dict)
	distance_matrix = pairwise_distances(rmsd_df)
	hdb = hdbscan.HDBSCAN(min_cluster_size=2, metric='precomputed')
	hdb.fit(distance_matrix)
	index_list = []
	for cluster in range(0, hdb.labels_.max() + 1):
		index_list.append([random.choice(np.where(hdb.labels_ == cluster)[0].tolist())])
	index_list.append(np.where(hdb.labels_ == -1)[0].tolist())
	index_list = [item for sublist in index_list for item in sublist]
	pdb_codes = rmsd_df.iloc[index_list].index.tolist()
	final_df = pd.DataFrame({'pdb_code':[pdb_code]*len(pdb_codes), 'filename':pdb_codes})

	return final_df

if __name__ == "__main__":

	cv_file = pd.read_csv("./cv_file.csv")
	pdb_list = cv_file['pdb_code'].tolist()
	pep_len_list = cv_file['peptide'].str.len().tolist()

	arg_list = []
	for i, pdb_code in enumerate(pdb_list):
		arg_list.append((pdb_code, pep_len_list[i]))

	#for i, pdb_file in enumerate(pdb_list):
	#	df_list.append(rmsd_calc(pdb_file, pep_len_list[i]))
	
	with WorkerPool(n_jobs=6) as pool:
		results = pool.map(rmsd_calc, arg_list, progress_bar=True)
	final_df = pd.concat(results, ignore_index=True)
	final_df.to_csv("Non_rendundant_structs.csv", index=False)
import os
import shutil
from collections import Counter
from tqdm import tqdm

from Bio import Align
from Bio import SeqIO

from pdbtools import pdb_rplchain, pdb_selchain, pdb_tofasta

import numpy as np
import pandas as pd
from biopandas.pdb import PandasPdb

standard_three_to_one_letter_code = {'ARG':'R', 'HIS':'H', 'LYS':'K', 'ASP':'D', 'GLU':'E', \
									 'SER':'S', 'THR':'T', 'ASN':'N', 'GLN':'Q', 'CYS':'C', \
									 'GLY':'G', 'PRO':'P', 'ALA':'A', 'VAL':'V', 'ILE':'I', \
									 'LEU':'L', 'MET':'M', 'PHE':'F', 'TYR':'Y', 'TRP':'W'}

def peptide_processing(peptide):
	sequence_list = list(peptide)
	for amino_acid in sequence_list:
		if (amino_acid not in standard_three_to_one_letter_code.values()):
			print("The provided amino acid " + str(amino_acid) + " in the peptide sequence is wrong! Aborting...")
			sys.exit(0)

def MHC_processing(allotype):
	for seq_record in SeqIO.parse("./helper_files/MHC_data.fasta", "fasta"):
		if seq_record.id == allotype:
			return
	print("The provided MHC allotype does not exist in our sequence DB and is not supported. Check the MHC naming, or omit the argument so that MHC allotype can be determined from the pdb file itself.")
	sys.exit(0)
	 
def initialize_dir(dir_name):
	if type(dir_name) == str:
		if os.path.exists(dir_name):
			for root, dirs, files in os.walk(dir_name):
				for f in files:
					os.unlink(os.path.join(root, f))
				for d in dirs:
					shutil.rmtree(os.path.join(root, d))
		else:
			os.umask(0)
			os.makedirs(dir_name,mode=0o777)
	elif type(dir_name) == list:
		for dir in dir_name:
			initialize_dir(dir)

def apply_function_to_file(func, input_filename, output_filename="", **kwargs):
	if output_filename == "": output_filename = input_filename

	overwritten = func(input_filename, **{key: value for key, value in kwargs.items() if key in func.__code__.co_varnames})
	overwritten = ''.join(overwritten)

	with open(output_filename, 'w') as output:
		output.write(overwritten)
	return output_filename

def copy_batch_of_files(src, dst, query):
	files = os.listdir(src)
	for f in files:
		if (query in f): shutil.copy(src + f, dst)

def filter_chains(pdb_file, chains, dst):
	filtered = pdb_selchain.run(open(pdb_file, 'r'), chains)
	with open(dst, 'w') as file:
		file.write(''.join(filtered))
	file.close()

def pdb_to_fasta(pdb_file, dst):
	fastaded = pdb_tofasta.run(open(pdb_file, 'r'), multi=True)
	with open(dst, 'w') as file:
		file.write(''.join(fastaded))
	file.close()

def remove_remarks_and_others_from_pdb(pdb_file, records=('ATOM', 'TER', 'END ')): 
	fhandle = open(pdb_file, 'r')
	for line in fhandle:
		if line.startswith(records):
			yield line
	fhandle.close()

def replace_HETATM(pdb_file):
	fhandle = open(pdb_file, 'r')
	for line in fhandle:
		if line.startswith("HETATM"):
			yield "ATOM  " + line[6:]
			continue
		yield line
	fhandle.close()

def check_chains(pdb_file):

	ppdb_peptide = PandasPdb()
	ppdb_peptide.read_pdb(pdb_file)
	pdb_df_peptide = ppdb_peptide.df['ATOM']
	chains = pd.unique(pdb_df_peptide['chain_id']).tolist()
	return chains

def retrieve_sequences_from_fasta(fasta_file):
	fasta_fields = []
	for seq_record in SeqIO.parse(fasta_file, "fasta"):
		fasta_fields.append(seq_record.seq)
	return(str(fasta_fields[1]), str(fasta_fields[0]))

def pretty_print_analytics(results_csv, filestore, verbose=True):

	results_csv.sort_values(by=['RankMHC Prediction'], inplace=True, ascending=False)
	results_csv.to_csv(filestore + '/results.csv', index=False)
	no_of_final_conformations = results_csv.shape[0]
	print("Total number of overall conformations: " + str(no_of_final_conformations))
	print("\nAnalytics:")
	print(results_csv.to_markdown(index=False))
	return results_csv

def sanitize_PDBs(structure_dir, filestore):
	structures = [f for f in os.listdir(structure_dir) if os.path.isfile(os.path.join(structure_dir, f))]
	for structure in structures:
		structure_file = structure_dir + structure
		apply_function_to_file(replace_HETATM, structure_file)		
		apply_function_to_file(remove_remarks_and_others_from_pdb, structure_file)
		chains = check_chains(structure_file)
		if (Counter(chains) == Counter(['A', 'B', 'C'])) or (Counter(chains) == Counter(['A', 'C'])):
			filter_chains(structure_file, ("A", "C"), filestore + '/sanitized_input/' + structure)
		else:
			print("Chains in your PDB files should be named as A, B and C! Aborting...")
			sys.exit(0)

def determine_sequence_and_allotype(structure_dir, filestore):
	structures = [f for f in os.listdir(structure_dir) if os.path.isfile(os.path.join(structure_dir, f))]
	peptide_list = []
	MHC_list = []
	for structure in structures:
		structure_file = structure_dir + structure
		pdb_to_fasta(structure_file, filestore + '/fasta_files/' + structure)
		(peptide, MHC) = retrieve_sequences_from_fasta(filestore + '/fasta_files/' + structure)
		peptide_list.append(peptide)
		MHC_list.append(MHC)
	if(len(list(set(peptide_list))) != 1 or len(list(set(MHC_list))) != 1):
		print("Different peptide or MHC in different files, these should be the same for proper ranking! Aborting...")
		sys.exit(0)
	best_record_allotype = None
	best_score = -float("inf")
	aligner = Align.PairwiseAligner()
	aligner.open_gap_score = -0.5
	aligner.extend_gap_score = 0.0
	for seq_record in tqdm(SeqIO.parse("./helper_files/MHC_data.fasta", "fasta")):		
		alignments = aligner.align(MHC, str(seq_record.seq))
		if(alignments[0].score > best_score):
			best_score = alignments[0].score
			best_record_allotype = seq_record.id
	return (peptide, best_record_allotype)
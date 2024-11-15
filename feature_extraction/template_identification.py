import pandas as pd
import numpy as np

import os
from itertools import zip_longest

from Bio.Align import substitution_matrices

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

def extract_anchors_PMBEC(peptide, MHC, frequencies):
	
	# Load and set the SMM matrix
	smm_matrix = pd.read_csv('./helper_files/PMBEC/' + MHC + '.txt', sep = '\t', skiprows=1, header = None, nrows=20).transpose()
	new_header = smm_matrix.iloc[0] #grab the first row for the header
	smm_matrix = smm_matrix[1:] #take the data less the header row
	smm_matrix.columns = new_header #set the header row as the df header
	smm_matrix = smm_matrix.reset_index().rename(columns={'index': 'position'})

	matrix = frequencies[(frequencies['allele'] == MHC) & (frequencies['length'] == 9)]
	
	# N-termini_candidate_1
	first_part_of_peptide = peptide[:3]
	pep_sequence = list(first_part_of_peptide)

	potential_pos_11 = smm_matrix[smm_matrix['position'] == 2][pep_sequence[0]].values[0]
	inertia_pos_11 = smm_matrix[smm_matrix['position'] == 1][pep_sequence[0]].values[0]
	will_11 = potential_pos_11 - inertia_pos_11
	potential_12 = smm_matrix[smm_matrix['position'] == 3][pep_sequence[1]].values[0]
	inertia_12 = smm_matrix[smm_matrix['position'] == 2][pep_sequence[1]].values[0]
	will_12 = potential_12 - inertia_12
	potential_pos_13 = smm_matrix[smm_matrix['position'] == 4][pep_sequence[2]].values[0]
	inertia_pos_13 = smm_matrix[smm_matrix['position'] == 3][pep_sequence[2]].values[0]
	will_13 = potential_pos_13 - inertia_pos_13
	
	# N-termini_candidate_3
	first_part = matrix[(matrix['position'].isin([1,2,3,4,5]))]
	anchor_1 = first_part.iloc[:, 5:].max(1).argmax() + 1
	if anchor_1 in [1, 4, 5]: # This is because it cannot happen
		anchor_1 = 2
	pep_sequence = list(peptide)

	potential_32 = smm_matrix[smm_matrix['position'] == anchor_1 - 1][pep_sequence[anchor_1 - 1]].values[0]
	inertia_32 = smm_matrix[smm_matrix['position'] == anchor_1][pep_sequence[anchor_1 - 1]].values[0]
	will_32 = potential_32 - inertia_32
	potential_pos_33 = smm_matrix[smm_matrix['position'] == anchor_1][pep_sequence[anchor_1]].values[0]
	inertia_pos_33 = smm_matrix[smm_matrix['position'] == anchor_1 + 1][pep_sequence[anchor_1]].values[0]
	will_33 = potential_pos_33 - inertia_pos_33
	potential_pos_35 = smm_matrix[smm_matrix['position'] == 4][pep_sequence[4]].values[0]
	inertia_pos_35 = smm_matrix[smm_matrix['position'] == 5][pep_sequence[4]].values[0]
	filling_pos_35 = smm_matrix[smm_matrix['position'] == 5][pep_sequence[5]].values[0]
	will_35 = potential_pos_35 - inertia_pos_35 + filling_pos_35
	
	# C-termini_candidate
	second_part_of_peptide = peptide[7:]
	second_part = matrix[(matrix['position'] >= 8) & (matrix['position'] <= len(peptide))]
	pep_sequence = list(second_part_of_peptide)
	stability_c = smm_matrix[smm_matrix['position'] == 9][pep_sequence[len(peptide) - 8]].values[0]
	min_will_c = float("inf")
	arg_will_c = 0
	for pos in range(8, len(peptide)):
		potential_pos_c = smm_matrix[smm_matrix['position'] == 9][pep_sequence[pos - 8]].values[0]
		intertia_pos_c = smm_matrix[smm_matrix['position'] == 8][pep_sequence[pos - 8]].values[0]
		will_pos_c = potential_pos_c - intertia_pos_c
		if min_will_c > will_pos_c:
			min_will_c = will_pos_c
			arg_will_c = pos

	anchor_1 = "2"
	if (will_11 < -0.25) and (will_12 < 0.25) and (will_13 < 0):
		anchor_1 = "1"
	if (will_32 < 0.0) and (will_33 < -0.5) and (will_35 < 0.5):
		anchor_1 = "3"
	anchor_2 = str(len(peptide))
	if (stability_c > 0.0) and (min_will_c < -0.25):
		anchor_2 = str(arg_will_c)
	return ",".join([anchor_1, anchor_2])

def predict_anchors_PMBEC(peptide, MHC):

	frequencies = pd.read_csv("./helper_files/mhcflurry.ba.frequency_matrices.csv")
	frequencies = frequencies[(frequencies['cutoff_fraction'] == 0.01)]
	frequencies['X'] = np.zeros(frequencies.shape[0])

	frequencies_alleles = os.listdir('./helper_files/PMBEC/') 
	frequencies_alleles = [x.split('.')[0] for x in frequencies_alleles]
	
	if (len(peptide) <= 7) or MHC not in frequencies_alleles:
		anchors = "2," + str(len(peptide))
		anchor_status = "Not Known"
	else:
		anchors = extract_anchors_PMBEC(peptide, MHC, frequencies)
		anchor_status = "Known"
	return anchors, anchor_status

def find_template(peptide, MHC_of_interest):

	# Determine anchors of peptide-MHC pair
	peptide_anchors, anchor_status = predict_anchors_PMBEC(peptide, MHC_of_interest)

	# Print all arguments
	print("Predicted Anchors:")
	print(peptide_anchors)

	# Pre-process anchors first
	peptide_anchors = [int(pos) for pos in peptide_anchors.split(",")]
	blosum_62 = substitution_matrices.load("BLOSUM62")

	# Continue with the alignment and template choice
	template_db = pd.read_csv("./helper_files/Template_DB.csv")[['pdb_code', 'MHC', 'peptide', 'peptide_length', 'Major_anchors']]
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
	similarity_matrix = pd.read_csv("./helper_files/9mer_similarity.csv")
	similarity_matrix = similarity_matrix[similarity_matrix['Allele'] == MHC_of_interest].T.reset_index().iloc[1:]
	similarity_matrix.columns = ['MHC', 'MHC_similarity']
	template_db = template_db.merge(similarity_matrix, on=['MHC'])

	# Final Choice
	template_db['Similarity_score'] = 0.5*template_db['MHC_similarity'] + 0.5*template_db['Peptide_similarity']
	template = template_db[template_db['Similarity_score'] == template_db['Similarity_score'].max()].sample(n=1)['pdb_code'].values[0]

	return template
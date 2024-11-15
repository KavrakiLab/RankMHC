import argparse

def RankMHC_parser():

	parser = argparse.ArgumentParser(description="RankMHC parser", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('input', type=str, nargs=1, help='Relative (or absolute)? Directory that all structures are kept')
	parser.add_argument('--peptide', type=str, default='', help='Peptide input sequence')
	parser.add_argument('--MHC', type=str, default='', help='MHC input sequence')
	parser.add_argument("--dir", type=str, default='intermediate_files', help='Location for all the intermediate files')
	parser.add_argument('--features', type=str, default='peptide_only', choices=['peptide_only', 'peptide_HLA', 'pairwise'], help='Different feature sets to use')
	parser.add_argument('--feat_transform', type=str, default='pooling', choices=['pooling', 'padding'], help='Transformation that brings all the features to the same length')
	parser.add_argument("--feat_terms", type=str, default='full', choices=['full', 'reduced'], help='Either use the full set of features or the ref2015 set (reduced)')
	parser.add_argument("--redundancy", type=str, default='redundant', choices=['non_redundant', 'redundant'], help='Use the models that were trained on data that included structurally redundant features')
	parser.add_argument("--objective", type=str, default='pairwise', choices=['ndcg', 'pairwise', 'regr'], help='Use the models based on a specific training objective (ranking, regression)')
	parser.add_argument("--sasa_method", type=str, default='naccess', choices=['naccess', 'freesasa'], help='Choose method to calculate SASA and RSA values per residue')
	parser.add_argument("--num_cores", type=int, default=8, help='Number of cores to use')
	
	return parser
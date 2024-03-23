import sys
sys.path.append("../Utils/")

import argparse as ap
from Utils_EditPaths import *

if __name__ == "__main__":
	parser = ap.ArgumentParser(description = "This is a script to sample sequences from the edit paths between different taxa pairs."+ 
											 "These sequences are used to perform classifier sensitivity testing. This script assumes, "+
											 "Get_EditPaths_Batches.py is run and has directories containing fasta files of edit paths and  "+
											 "a ccorresponding summary file of the edit paths. This script retruns a fasta file of worst case "+
											 "adversarial sequences and further samples multiple sequences between taxa pairs at various points "+
											 "from the decision boundary upon passing the correct parameters. Get_EditPaths_Batches.py produces "+
											 "<directory>/Edit_Paths/<taxa_family>/ containing fasta files of edit paths and "+
											 "<directory>/Edit_Paths_Summary/<taxa_family>/ containing summary files of edit paths for every taxon "+
											 "at the family in the database.")

	parser.add_argument("-e", "--edit_path_directory", help="A directory containing fasta files of edit paths.", required=True)
	parser.add_argument("-s", "--edit_path_summary", help="A directory containing summary file of edit paths.", required=True)
	parser.add_argument("-o", "--output_directory", help="Location of the output directory", required=True)
	parser.add_argument("-u", "--summary_path", help = "Location to write worst-case adversarial sequences", required = True)
	parser.add_argument("--sample_multiple_seqs", help = "Flag to sample multiple sequences at various distances from the decision boundary.", 
						action='store_true')
	parser.add_argument("-c", "--num_seqs", help = "Number of sequences to sample if --sample_multiple_seqs is set to True.", default = "5", required = False)
	parser.add_argument("--all", help = "Get all adversarial sequences.", action='store_true')
	
	args = parser.parse_args()

	edit_path_sequences_dir = args.edit_path_directory
	edit_path_summary_dir = args.edit_path_summary
	outdir = args.output_directory
	sample_multiple_seq = args.sample_multiple_seqs
	seq_count = int(args.num_seqs)
	summary_file = args.summary_path
	all = args.all

	assert isdir(edit_path_sequences_dir),  "Check input directories. "+edit_path_sequences_dir+" not found"
	assert isdir(edit_path_summary_dir),  "Check input directories. "+edit_path_summary_dir+" not found"
	assert (sample_multiple_seq == True and seq_count > 0) or (sample_multiple_seq == False), "Invalid choice of num_seqs. num_seqs has to be an integer > 0."

	if not isdir(outdir):
		try: mkdir(outdir)
		except OSError:
			print("Cannot create ", outdir)
			sys.exit(1)

	files = listdir(edit_path_sequences_dir)
	adversarial_sequences = {}
	Worstcase_Sequences_Summary = []
	for f in files:
		if (f.endswith(".fasta")) and not(f.startswith(".")):
			edit_path_seq = edit_path_sequences_dir+f
			edit_summary_file = edit_path_summary_dir + f.replace(".fasta", ".csv")
			if isfile(edit_path_seq) and isfile(edit_summary_file):
				W = Get_Worst_Case_Adversarial_Sequences(edit_summary_file)
				W = W+(f.replace(".fasta",""),)
				sa, sb = W[0].replace("/","-"), W[1].replace("/","-")
				sampled_seqs = Sample_Adversarial_Examples(edit_path_seq, W, sample_multiple_seq, seq_count, all)
				try: adversarial_sequences[(sa, sb)].update(sampled_seqs)
				except KeyError: adversarial_sequences[(sa, sb)] = sampled_seqs
				Worstcase_Sequences_Summary.append(W)
			else:
				if isfile(edit_path_seq) == "False": print(edit_path_seq, " not found")
				if isfile(edit_summary_file) == "False": print(edit_summary_file, " not found")

	for k in adversarial_sequences:
		out_file = outdir+k[0]+'_'+k[1]+'.fasta'
		Write_DB(adversarial_sequences[k], out_file)
	
	df_adversarial_seqs_summary = pd.DataFrame(Worstcase_Sequences_Summary, columns = ['GA','GB','Path','DA','DB','Edits','Edit_Path'])
	df_adversarial_seqs_summary.to_csv(summary_file, sep = "\t")
import sys
sys.path.append("../Utils/")

import argparse as ap
from Utils_EditPaths import *

if __name__ == "__main__":
	parser = ap.ArgumentParser(description = "This is a script to generate parameters for sampling edit paths "+ 
											 "between pairs of taxa that have the \"family\" level classification. This "+
											 "script takes a fasta file of sequences as input and runs the RDP classifier and "+
											 "returns a set of text files that contain the parameters for running Get_EditPaths_Batches.py "+
											 "concurrently on the SLURM cluster as an array job.")

	parser.add_argument("-i", "--input", help="A fasta file of sequences.", required=True)
	parser.add_argument("-o", "--output", help="Directory to write outputs to.", required=True)
	parser.add_argument("-rdp", "--rdp_classifier", help="Directory containing the RDP classifier program.", required=True)
	parser.add_argument("-db", "--rdp_DB", help = "Directory containing the trained RDP classifier.", required=True)
	parser.add_argument("-pre","--prefix", help = "Prefix for the files containing the parameters. (Default = Batch)", required = False, default = "Batch")
	parser.add_argument("-b","--Batch_Size", help = "Size of the batch. (Default = 1000)", required = False, default = "1000")

	parser.add_argument("-c","--confidence", help = "Minimum score of the classifier to consider for generating edit paths. (Default = 0.50)", required = False, default = "0.50")
	parser.add_argument("-NTAX","--number_of_taxa", help = "Maximum number of taxa to consider, to enable tractability. (Default = 100)", required = False, default = "100")
	
	args = parser.parse_args()

	fasta_file = args.input
	out_dir = args.output
	Prog_RDP = args.rdp_classifier #'/fs/cbcb-scratch/hsmurali/RDP_Outlier_Analysis/RDPTools/'
	RDP_DB = args.rdp_DB #'/fs/cbcb-lab/mpop/projects/RDP_Outlier_Analysis/RDP_Database/rRNAClassifier.properties'
	pre = args.prefix

	confidence = float(args.confidence)
	NTAX = int(args.number_of_taxa)
	batch_size = int(args.Batch_Size)

	if not isdir(out_dir):
		try:
			mkdir(out_dir)
		except OSError:
			print('Unable to create '+out_dir)
			sys.exit(1)

	head, tail = split(fasta_file)
	RDP_out = out_dir+'/'+tail.replace(".fasta",".tax")
	Run_Classifier((Prog_RDP, fasta_file, RDP_DB, RDP_out))

	fasta_file_seqs = out_dir + 'Sequences/'
	if not isdir(fasta_file_seqs):
		try: mkdir(fasta_file_seqs)
		except OSError:
			print('Unable to create '+fasta_file_seqs)
			sys.exit(1)

	blast_out_dir = out_dir + 'BLAST_Hits/'
	if not isdir(blast_out_dir):
		try: mkdir(blast_out_dir)
		except OSError:
			print('Unable to create '+blast_out_dir)
			sys.exit(1)

	batch_out_dir = out_dir + 'Batch_Parameters/'
	if not isdir(batch_out_dir):
		try: mkdir(batch_out_dir)
		except OSError:
			print('Unable to create '+batch_out_dir)
			sys.exit(1)

	BLAST_Jobs = Select_Candidates_Pairwise_BLAST(RDP_out, fasta_file, fasta_file_seqs, blast_out_dir, confidence, 'family', NTAX)
	ctr = 0
	batch_id = 0
	parameters = []

	for b in BLAST_Jobs:
		parameters.append(str(b)+'\n')
		ctr += 1
		if ctr >= batch_size:
			buf = open(batch_out_dir+pre+"."+str(batch_id)+".txt","w")
			buf.writelines(parameters)
			buf.close()
			batch_id += 1
			ctr = 0
			parameters = []
	buf = open(batch_out_dir+pre+"."+str(batch_id)+".txt","w")
	buf.writelines(parameters)
	buf.close()
import sys
sys.path.append("../Utils/")

import argparse as ap
from Utils_Transitive_Annotation import *

if __name__ == "__main__":
	parser = ap.ArgumentParser(description = "This is a script to split a fasta file into smaller batch of sequences. This "+
											 "is particularly helpful to run the RDP classifier in multiple threads.")
	parser.add_argument("-q", "--query", help="A fasta file of query sequences.", required=True)
	parser.add_argument("-d", "--directory", help="directory to write the split up fasta file.", required=True)
	parser.add_argument("-b", "--batch_size", help="size of the smaller fasta files", required=False, default = "10")

	args = parser.parse_args()
	assert isfile(args.query), "Incorrect input file."

	output_dir = args.directory
	if not isdir(output_dir):
		try: mkdir(output_dir)
		except OSError:
			print("Unable to create ",output_dir)
			sys.exit(1)

	batches = int(args.batch_size)
	seqs = Load_Database(args.query, keep_header = False)
	print(len(seqs))
	ctr = 0
	batch = 0
	buf = {}
	for s in seqs:
		buf[s] = seqs[s][:500]
		ctr += 1
		if ctr == batches:
			Write_DB(buf, output_dir+"Batch."+str(batch)+".fasta")
			ctr = 0
			batch+=1
			buf = {}
	Write_DB(buf, output_dir+"Batch."+str(batch)+".fasta")

	
	
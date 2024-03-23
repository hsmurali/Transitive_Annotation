import sys
sys.path.append("../Utils/")

import argparse as ap
from Utils_Transitive_Annotation import *

if __name__ == "__main__":
	parser = ap.ArgumentParser(description = "This is a script to sample query sequences to be used as controls.")
	parser.add_argument("-q", "--query", help="A fasta file of sequences or a directory of fasta files containing split up query sequences.", required=True)
	parser.add_argument("-db", "--rdp_DB", help = "Directory containing the trained RDP database.", required=True)
	parser.add_argument("-rdp", "--rdp_classifier", help="Directory containing the RDP classifier program.", required=True)
	parser.add_argument("-a", "--adversarial_sequences", help="Directory containing fasta files of adversarial sequences.", required=True)
	parser.add_argument("-taxa","--taxa", help = "Taxa file for training the RDP classifier", required = True)
	parser.add_argument("-o", "--output", help="Directory to write outputs to.", required=True)

	args = parser.parse_args()

	RDP_path = args.rdp_classifier #'/fs/cbcb-scratch/hsmurali/RDP_Outlier_Analysis/RDPTools/'
	database = args.rdp_DB #'/fs/cbcb-lab/mpop/projects/RDP_Outlier_Analysis/RDP_Database/trainset18_062020.fa'
	taxa_path = args.taxa #'/fs/cbcb-lab/mpop/projects/RDP_Outlier_Analysis/RDP_Database/trainset18_db_taxid.txt'
	fna_path = args.query
	adv_seqs = args.adversarial_sequences
	outdir = args.output
	RDP_out_path = outdir+'RDP.tax.out'

	if not isdir(outdir):
		mkdir(outdir)

	DB = Load_Database(database, keep_header = True)
	d_taxa = Extract_Taxonomic_Annotation(list(DB.keys()))
	#Taxa_Counts, Taxa_Dict = Load_Taxa_Information(list(DB.keys()))
	
	#Train_Classifier(RDP_path, database, taxa_path, outdir) 
	#Run_Classifier((RDP_path, fna_path, outdir+'rRNAClassifier.properties', RDP_out_path))

	Query_Seqs = Load_Database(fna_path)
	df = Parse_RDP_Output(RDP_out_path)

	seqs_of_interest = []
	for f in listdir(adv_seqs):
		if not isdir(adv_seqs+f+'/'): continue
		for s in listdir(adv_seqs+f+'/'):
			if s.startswith('.'): continue
			if s.endswith(".csv"):
				splits = s.replace(".csv","").split('_')
				seqs_of_interest += [splits[0]+'_'+splits[1], splits[2]+'_'+splits[3]]
	seqs_of_interest = list(set(seqs_of_interest))
	df_filter = df.loc[seqs_of_interest]

	df_filter = df_filter.sample(1000)
	sequences = df_filter.index.tolist()
	genus = df_filter['genus'].tolist()
	p_genus = df_filter['p-genus'].tolist()

	Adversarial_Dataset = {}

	for j in range(len(sequences)):
		g = genus[j].replace(" ","")
		Adversarial_Dataset['Control_Seq_'+sequences[j]+'_'+g+'_'+str(p_genus[j])+'\t'+d_taxa[g.replace("-","/")]] =  Reverse_Compliment(Query_Seqs[sequences[j]])
	Write_DB(Adversarial_Dataset, outdir+'/Control_Seqs.fasta')
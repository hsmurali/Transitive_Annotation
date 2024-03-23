import sys
sys.path.append("../Utils/")

import shutil
import argparse as ap
from copy import deepcopy
from Utils_Transitive_Annotation import *
		
if __name__ == "__main__":
	parser = ap.ArgumentParser(description = "This is a script to perform classfier senstivity training by adding "+
											 "in transitively annotated query sequences.")

	parser.add_argument("-q", "--query", help="A fasta file of sequences or a directory of fasta files containing split up query sequences.", required=True)
	parser.add_argument("-o", "--output", help="Directory to write outputs to.", required=True)
	parser.add_argument("-d", "--RDP_wo_adversarial", help = "Results of running the RDP classifer without any adversarial example. Required if --run_vanilla is disabled.",
						required = False, default = "")
	parser.add_argument("-s","--sample", help="A fasta file to add query sequences as control", required = True)
	parser.add_argument("-rdp", "--rdp_classifier", help="Directory containing the RDP classifier program.", required=True)
	parser.add_argument("-db", "--rdp_DB", help = "Directory containing the trained RDP database.", required=True)
	parser.add_argument("-taxa","--taxa", help = "Taxa file for training the RDP classifier", required = True)
	parser.add_argument("--split", help="Flag to indicate if a directory of split query sequences is passed.", action='store_true')	
	parser.add_argument("--reverse", help = "Add the reverse compliment of adversarial sequence to the training data.", action='store_true')
	parser.add_argument("--run_vanilla", help = "Flag to indicate if the RDP classifier is run without adversarial sequences",action='store_true')
	parser.add_argument("-nt", "--threads", help = "Number of threads to run the RDP classifier", default = "32", required = False)
	args = parser.parse_args()

	RDP_path = args.rdp_classifier #'/fs/cbcb-scratch/hsmurali/RDP_Outlier_Analysis/RDPTools/'
	database = args.rdp_DB #'/fs/cbcb-lab/mpop/projects/RDP_Outlier_Analysis/RDP_Database/trainset18_062020.fa'
	taxa_path = args.taxa #'/fs/cbcb-lab/mpop/projects/RDP_Outlier_Analysis/RDP_Database/trainset18_db_taxid.txt'
	split = args.split

	vanilla = args.run_vanilla
	threads = int(args.threads)
	reverse = args.reverse
	sample_ = args.sample
	query_path = args.query #'/fs/cbcb-lab/mpop/projects/RDP_Outlier_Analysis/Query_Split/'
	Classifier_WO_Adversarial = args.RDP_wo_adversarial
	Ad_DB_Dir = args.output

	assert split == True and isdir(query_path), "Split enabled but query directory not found"
	assert isdir(RDP_path) == True, "RDP classifier path is incorrect"
	assert isfile(database) == True, "RDP Database path is incorrect"
	assert isfile(taxa_path) == True, "Taxa file for RDP classifier is incorrect"
	assert (vanilla == False and isfile(Classifier_WO_Adversarial)) or (vanilla == True), "Incorrect RDP_wo_adversarial file"

	if not isdir(Ad_DB_Dir):
		try: mkdir(Ad_DB_Dir)
		except OSError:
			print("Error creating ",Ad_DB_Dir)
			sys.exit(1)

	####Loading the RDP Database
	DB = Load_Database(database, keep_header = True)
	Taxa = Extract_Taxonomic_Annotation(list(DB.keys()))
	Taxa_Counts, Taxa_Dict = Load_Taxa_Information(list(DB.keys()))

	####Loading the output from running the RDP classifier without adversarial examples. 
	if vanilla == True:
		Classifier_WO_Adversarial = Ad_DB_Dir+"WO_Adversarial_Training.tax.out"
		if split == False:
			Train_Classifier(RDP_path, database, taxa_path, Ad_DB_Dir)
			args = (RDP_path, query_path, Ad_DB_Dir+'rRNAClassifier.properties', Classifier_WO_Adversarial)
			Run_Classifier(args)
		elif split == True:
			Train_Classifier(RDP_path, database, taxa_path, Ad_DB_Dir)
			print('Done')
			Run_Classifier_Multiple_Threads(RDP_path, Ad_DB_Dir+'/TMP/', query_path, Ad_DB_Dir+'rRNAClassifier.properties', num_threads = threads)
			Merge_Files(Ad_DB_Dir+'/TMP/', Classifier_WO_Adversarial)
	try: df_wo_adv = Parse_RDP_Output(Classifier_WO_Adversarial)
	except:
		print("Error loading the output of RDP...")
		sys.exit(1)

	Ad_DB = Load_Database(sample_, keep_header = True)
	seq_keys = list(Ad_DB.keys())

	w = open(Ad_DB_Dir+'Adverserial_Trained.summary','w')
	out = pd.DataFrame()

	for k in seq_keys:
		DB_Adverserial = deepcopy(DB)
		rand_key = k.split('\t')[0]
		GA = k.split('\t')[-1].split(";")[-1]

		if not isfile(Ad_DB_Dir+'Adverserial_Trained_'+rand_key+".tax.out"):
			if isdir(Ad_DB_Dir+'/TMP/'): shutil.rmtree(Ad_DB_Dir+'/TMP/') 
			
			try: rand_seq = Ad_DB[k]
			except KeyError:
				print(k, "Key not found!!!")
				sys.exit(1)

			rand_key = rand_key.replace(" ","").replace("/","-")
			if reverse: S = {k:Reverse_Compliment(rand_seq)}
			else: S = {k:rand_seq}
			DB_Adverserial = Update_Database(DB_Adverserial, S) 
			Write_DB(DB_Adverserial, Ad_DB_Dir+'Adverserial_DB.fasta')
			Train_Classifier(RDP_path, Ad_DB_Dir+'Adverserial_DB.fasta', taxa_path, Ad_DB_Dir)

			if split == True:
				Run_Classifier_Multiple_Threads(RDP_path, Ad_DB_Dir+'/TMP/', query_path, Ad_DB_Dir+'rRNAClassifier.properties', num_threads = threads)
				Merge_Files(Ad_DB_Dir+'/TMP/', Ad_DB_Dir+'Adverserial_Trained_'+rand_key+".tax.out")
			elif split == False:
				args = (RDP_path, query_path, Ad_DB_Dir+'rRNAClassifier.properties', Ad_DB_Dir+'Adverserial_Trained_'+rand_key+".tax.out")
				Run_Classifier(args)

			remove(Ad_DB_Dir+'Adverserial_DB.fasta')

		rand_key = rand_key.replace(" ","").replace("/","-")
		df_adv = Parse_RDP_Output(Ad_DB_Dir+'Adverserial_Trained_'+rand_key+".tax.out")
		summary = Compare_Classifiers(df_wo_adv, df_adv, GA, "", Taxa_Counts, Taxa_Dict)
		summary['SeqID'] = rand_key
		df_label_flips = Summarize_Label_Flips(df_wo_adv, df_adv, rand_key, Taxa_Dict)
		out = out.append(df_label_flips, ignore_index = True)
		w.write(str(summary)+'\n')

	w.close()
	out.to_csv(Ad_DB_Dir+'Label_Flips.txt', sep = "\t")
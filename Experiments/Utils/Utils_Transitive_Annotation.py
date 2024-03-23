import sys
import subprocess
import numpy as np
import pandas as pd
import multiprocessing
from os import popen, listdir, remove, mkdir 
from os.path import isfile, split, realpath, isdir,dirname

def Calculate_LCA(taxa_list, rdp_taxa):
	taxa_filtered = []
	for t in taxa_list:
		taxa_filtered.append(rdp_taxa[t.replace("/","-").replace(" ","")])
	LCA_found = False
	while not (LCA_found):
		for l in ['G','F','O','C','P','K','D']:
			diff = False
			for i in range(0, len(taxa_filtered)-1):
				try:
					if taxa_filtered[i][l] != taxa_filtered[i-1][l]:
						diff = True
						break
				except KeyError:
					diff = True
					break
			if diff == False:
				LCA = l
				LCA_found = True
				break
		if l == 'D' and LCA_found == False:
			LCA = 'D'
			LCA_found = True
	return(LCA)

def Parse_RDP_Output(filedir):
	op = []
	try:
		lines = open(filedir).readlines()
	except FileNotFoundError:
		return pd.DataFrame()
	for l in lines:
		l = l.replace('\n','')
		l = l.split("\t")
		d = {'Seq-ID':l[0]}
		for i in range(2, len(l), 3):
			try:
				d[l[i+1]] = l[i].replace(" ","").replace("/","-")
				d['p-'+l[i+1]] = float(l[i+2])
			except IndexError:
				d['p-genus'] = 0.0
				print(l)
		op.append(d)
	df = pd.DataFrame(op)
	df = df.set_index('Seq-ID')
	return df

def Diff(row):
	if row['genus(wo-adv)'] == row['genus(adv)']:
		return 0
	return 1

def Load_Taxa_Information(D):
	Taxa_Counts = {}
	Taxa_Dict = {}
	levels = ['D','K','P','C','O','F','G']

	for r in D:
		r = r.split('\t')
		taxa = r[1].replace("\n","").replace("/","-").replace(" ","").split(";")
		d = dict(zip(levels[:len(taxa)], taxa))
		genus = taxa[-1]
		Taxa_Dict[genus] = d
		try:
			Taxa_Counts[genus] += 1
		except KeyError: 
			Taxa_Counts[genus] = 1
	return Taxa_Counts, Taxa_Dict
	
def Compare_Classifiers(df_wo_adv, df_adv, GA, GB, Taxa_Counts, Taxa_Dict):
	df_wo_adv = df_wo_adv[['genus','p-genus']].rename(columns = {'genus':'genus(wo-adv)', 'p-genus':'p-genus(wo-adv)'})
	df_adv = df_adv[['genus', 'p-genus']].rename(columns = {'genus':'genus(adv)', 'p-genus':'p-genus(adv)'})
	df = df_wo_adv.join(df_adv)
	df['Diff'] = df.apply(Diff, axis = 1)
	df_filt = df[df['Diff'] == 1]
	n_taxa_flipped = len(df_filt['genus(wo-adv)'].unique())
	flipped_taxa_list = list(df_filt['genus(wo-adv)'].unique())
	Largest_Taxa = df_filt.groupby(['genus(wo-adv)']).count().reset_index()
	max_val = Largest_Taxa['Diff'].max()
	try:
		Largest_Taxa = Largest_Taxa[Largest_Taxa['Diff'] == max_val].iloc[0]['genus(wo-adv)']
	except IndexError:
		Largest_Taxa = ""
		max_val = -1
	
	d = {'GA':GA, 'GB':GB,'Flips':df['Diff'].sum(), 'N_Flipped_Taxa':n_taxa_flipped, 
		 'Flipped_Taxa_List':flipped_taxa_list, 'Largest-Flipped-Taxa':Largest_Taxa, 
		 'n-Largest-Taxa-Counts':max_val, 'LCA_Flipped_Taxa':Calculate_LCA(flipped_taxa_list, Taxa_Dict)}

	try:
		d['GA_Taxa_Counts(DB)'] = Taxa_Counts[GA.replace("/","-")]
	except KeyError:
		print('GA', GA, GB)

	try:
		d['GB_Taxa_Counts(DB)'] = Taxa_Counts[GB.replace("/","-")]
	except KeyError:
		print('GB', GA, GB)
	
	df_adv_GA = df_adv[df_adv['genus(adv)'] == GA]
	df_adv_GB = df_adv[df_adv['genus(adv)'] == GB]
	
	df_wo_adv_GA = df_wo_adv[df_wo_adv['genus(wo-adv)'] == GA]
	df_wo_adv_GB = df_wo_adv[df_wo_adv['genus(wo-adv)'] == GB]
	
	if GA != "" and GB != "":
		d .update({'GA':GA, 'GB':GB, '#Seqs-GA(wo-adv)':len(df_wo_adv_GA), 
				   '#Seqs-GB(wo-adv)':len(df_wo_adv_GB),
				   'Avg-P-GA(wo-adv)':df_wo_adv_GA['p-genus(wo-adv)'].mean(),
				   'Avg-P-GB(wo-adv)':df_wo_adv_GB['p-genus(wo-adv)'].mean(),
				   '#Seqs-GA(adv)':len(df_adv_GA), '#Seqs-GB(adv)':len(df_adv_GB),
				   'Avg-P-GA(adv)':df_adv_GA['p-genus(adv)'].mean(),
				   'Avg-P-GB(adv)':df_adv_GB['p-genus(adv)'].mean()})
	return d

def Load_Database(seq_path, keep_header = False):
	d = {}
	header, seq = "", ""
	fasta_sequences = open(seq_path).readlines()
	for s in fasta_sequences:
		if s.startswith(">"):
			d[header] = seq.replace("-","")
			header = s.replace("\n","").replace(">","")
			if keep_header == False:
				header = header.split(' ')[0]
			seq = ""
		else:
			seq += s.replace("\n","")
	d[header] = seq.replace("-","")
	header = s.replace("\n","").replace(">","")
	del d[""]
	return d

def Update_Database(DB, seqs):
	for s in seqs:
		DB[s] = seqs[s]
	return DB

def Write_DB(DB, out_path):
	o = open(out_path,'w')
	for d in DB:
		o.write('>'+d+'\n')
		o.write(DB[d]+'\n')
	o.close()

def Sample(query_path, df_rdp, n):
	D = {}
	if isdir(query_path):
		files = listdir(query_path)
		for f in files:
			d = Load_Database(query_path+f, keep_header = True)
			D = Update_Database(D, d)
	elif isfile(query_path):
		D = Load_Database(query_path, keep_header = True)
	df_frac = df_rdp.sample(n = n)
	sequence_ids = df_frac.index.tolist()
	taxa = df_frac['genus'].tolist()
	
	adversarial_controls = {}
	for i in range(len(sequence_ids)):
		s = sequence_ids[i]
		g = taxa[i]
		adversarial_controls[(s,g)] = D[s]
	return adversarial_controls

def Extract_Taxonomic_Annotation(seq_headers):
	d = {}
	for s in seq_headers:
		try:
			t = s.split("\t")[1]
			genus = t.split(";")[-1]
			d[genus.replace(" ","")] = t
		except IndexError:
			print(s)
	return d

def Run_Classifier(args):
	RDP_Path = args[0] 
	query_path = args[1] 
	classifier_path = args[2] 
	out_file = args[3]

	is_RDP = isdir(RDP_Path)
	is_query = isfile(query_path)
	is_classifier = isfile(classifier_path)

	assert is_RDP == True, "Classfier missing from, " + str(RDP_Path)
	assert is_query == True, "Query file missing from, " + str(query_path)
	assert is_classifier == True,"Classifier file missing from, " + str(classifier_path)
		
	cmd = 'java  -Xmx10g -jar '+RDP_Path+'classifier.jar classify -o '+out_file+' -q '+query_path+' -t '+classifier_path
	subprocess.Popen(cmd,shell=True).wait()
	
	assert isfile(out_file) == True, "Error running the RDP Classifier"

def Merge_Files(data_dir, out_file):
	files = listdir(data_dir)
	out_buf = []
	for f in files:
		lines = open(data_dir+f).readlines()
		for l in lines:
			if l.endswith("\n"): 
				out_buf.append(l)
			else: 
				out_buf.append(l+"\n")
	b = open(out_file,'w')
	b.writelines(out_buf)
	b.close()

	
	command = "rm -rf "+data_dir
	subprocess.Popen(command,shell=True).wait()
	assert isdir(out_file) == False, "Error removing "+data_dir
	
def Run_Classifier_Multiple_Threads(RDP_Path, out_dir, data_dir, classifier_path, num_threads = 32):
	command_args = []
	seqs = listdir(data_dir)

	if not isdir(out_dir):
		mkdir(out_dir)
	for s in seqs:
		if s.endswith(".fasta") or s.endswith(".fa") or s.endswith(".fna"):
			arg = (RDP_Path, data_dir+s, classifier_path, out_dir+s.replace(".fasta",".out").replace(".fa",".out").replace(".fna",".out"))
			command_args.append(arg)
	pool = multiprocessing.Pool(int(num_threads))
	result = pool.map(func=Run_Classifier, iterable=command_args)
	pool.close()
	pool.join()

def Train_Classifier(RDP_Path, db_path, taxa_path, out_path):
	is_RDP = isdir(RDP_Path)
	is_db = isfile(db_path)
	is_taxa = isfile(taxa_path)

	if not isdir(out_path):
		mkdir(out_path)
		
	assert is_RDP == True, "Classfier missing from, " + str(RDP_Path)
	assert is_db == True, "Query file missing from, " + str(db_path)
	assert is_taxa == True,"Classifier file missing from, " + str(taxa_path)

	cmd = 'java -Xmx10g -jar '+RDP_Path+'classifier.jar train -o '+out_path+' -s '+db_path+' -t '+taxa_path
	subprocess.Popen(cmd,shell=True).wait()
	
	assert isdir(out_path) == True, "Error training the RDP classifier"

	f = open(out_path+"/rRNAClassifier.properties",'w')
	
	assert isfile(out_path+'bergeyTrainingTree.xml'), "bergeyTrainingTree.xml missing"
	f.write('bergeyTree=bergeyTrainingTree.xml\n\n')

	assert isfile(out_path+'genus_wordConditionalProbList.txt'), "genus_wordConditionalProbList.txt missing"
	f.write("probabilityList=genus_wordConditionalProbList.txt\n\n")
	
	assert isfile(out_path+'wordConditionalProbIndexArr.txt'), "wordConditionalProbIndexArr.txt missing"
	f.write('probabilityIndex=wordConditionalProbIndexArr.txt\n\n')
	
	assert isfile(out_path+'logWordPrior.txt'), "logWordPrior.txt missing"
	f.write('wordPrior=logWordPrior.txt\n\n')
	f.write('classifierVersion=RDP Naive Bayesian rRNA Classifier Version 2.5, May 2012\n\n') 
	
	f.close()

def Reverse_Compliment(S):
	Rev = {'A':'T','T':'A','G':'C','C':'G','-':'N','N':'N'}
	r = ""
	for s in S:
		r += Rev[s.upper()]
	return r[::-1]

def Get_LCA_Flips(row, rdp_taxa):
	taxa_1 = row['genus(wo-adv)']
	taxa_2 = row['genus(adv)']
	return Calculate_LCA([taxa_1, taxa_2], rdp_taxa)

def Summarize_Label_Flips(df_wo_adv, df_adv, seq_id, rdp_taxa):
	df_wo_adv = df_wo_adv[['genus']].rename(columns={'genus':'genus(wo-adv)'})
	df_adv = df_adv[['genus']].rename(columns={'genus':'genus(adv)'})
	df_joined = df_wo_adv.join(df_adv)
	df_joined['Count'] = 0
	df_joined = df_joined.groupby(['genus(wo-adv)','genus(adv)']).count()
	df_joined['Seq_ID'] = seq_id
	df_joined = df_joined.reset_index()
	df_joined['LCA'] = df_joined.apply(Get_LCA_Flips, rdp_taxa = rdp_taxa, axis = 1)
	return df_joined

def Extract_Parameter(S, filt):
	splits = S.split('_')
	for s in splits:
		if s.startswith(filt+"="):
			s = s.split('-')[0]
			return int(s.replace(filt+"=",""))

def Extract_SA_SB(Seq_Id):
	S = Seq_Id.split('-')
	splits = S[0].split('_')
	SA = splits[0]+'_'+splits[1]
	SB = splits[2]+'_'+splits[3]
	return pd.Series({'SA':SA,'SB':SB})
from Utils_Transitive_Annotation import *
import re
import random
import json

def Validate_BLAST_Hits(jsonOut, align_thresh):
	jf = None
	flag = True
	try: jf = open(jsonOut, "r")
	except OSError as err: 
		print("Cannot open json: " + err.strerror)
		flag = False

	aln = None
	try: aln = json.load(jf)
	except json.JSONDecodeError as err:
		print("Could not decode json input at line {l}: {e}".format(l=err.lineno, e=err.msg))
		flag = False

	alen, qlen = 0, int(aln["BlastOutput2"][0]["report"]["results"]["bl2seq"][0]["query_len"])
	
	try:
		if aln["BlastOutput2"][0]["report"]["results"]["bl2seq"][0]["message"] == "No hits found":
			flag = False
	except KeyError:
		 alen = int(aln["BlastOutput2"][0]["report"]["results"]["bl2seq"][0]["hits"][0]["hsps"][0]["align_len"])

	if alen == 0 or qlen / alen < align_thresh: flag = False
		
	firstHit = aln["BlastOutput2"][0]["report"]["results"]["bl2seq"][0]["hits"][0]["hsps"][0]
	query, hit = firstHit["qseq"], firstHit["hseq"]

	if len(query) != len(hit):
		print("Query and hit must have same length\n", query,"\n", hit)
		flag = False
	return flag, firstHit

def Select_Candidates_Pairwise_BLAST(rdp_path, fasta_file, fasta_output_dir, blast_output_dir,confidence, level, NTAX):
	db = Load_Database(fasta_file, False)
	df = Parse_RDP_Output(rdp_path)
	
	if level not in df.columns.tolist():
		print(level, 'not valid...')
		sys.exit(1)
	df_genus = df[[level, 'genus', 'p-genus']]
	df_genus = df_genus[df_genus['p-genus'] >= confidence].reset_index()
	taxa = df_genus.groupby(['family'])['Seq-ID'].apply(list).to_dict()
	seqAnn = dict(zip(df_genus['Seq-ID'].tolist(), df_genus['genus'].tolist()))
	BLAST_Jobs = []

	for rank in taxa.keys():
		try: 
			mkdir(fasta_output_dir+rank+'/')
			mkdir(blast_output_dir+rank+'/')
		except OSError:
			print("Could not create output dir for taxon " + rank + ": ")
			pass

		if len(taxa[rank]) > NTAX:  # make sure things are tractable
			taxa[rank] = random.sample(taxa[rank], NTAX)

		for i in range(len(taxa[rank])-1):
			for j in range(i+1, len(taxa[rank])):
				if not seqAnn[taxa[rank][i]] == seqAnn[taxa[rank][j]]:
					d = {taxa[rank][i]: db[taxa[rank][i]]}
					fasta_1_path = fasta_output_dir+rank+'/'+taxa[rank][i]+'.fasta'
					Write_DB(d, fasta_1_path)

					d = {taxa[rank][j]: db[taxa[rank][j]]}
					fasta_2_path = fasta_output_dir+rank+'/'+taxa[rank][j]+'.fasta'
					Write_DB(d, fasta_2_path)
					
					jsonOut = blast_output_dir + "/" + rank + "/" + taxa[rank][i] + "_" + taxa[rank][j] + ".json"
					BLAST_Jobs.append((fasta_1_path, fasta_2_path, rank, jsonOut))
	return BLAST_Jobs

def Run_BLAST(fasta_1, fasta_2, jsonOut):
	assert isfile(fasta_1) == True, fasta_1+" missing!!!"
	assert isfile(fasta_2) == True, fasta_2+" missing!!!"

	try: subprocess.run(["blastn", "-query", fasta_1, "-subject", fasta_2, "-outfmt", "15", "-out", jsonOut], check=True)
	except subprocess.CalledProcessError as err:
		print("Could not run blastn: " + err.output)
		sys.exit(1)

def editPath(firstHit, output_path, npath):
	query, hit = firstHit["qseq"], firstHit["hseq"]
	
	mutations = []
	for i in range(len(query)):
		if query[i] != hit[i]: # found mutation
			mutations.append([i, hit[i]])

	mutationOrder = list(range(len(mutations)))
	savedMutations = mutations.copy()

	n = 1
	out_d = {}
	for path in range(1, npath + 1):
		mutations = savedMutations.copy()
		random.shuffle(mutationOrder)
		currSeq = query

		for i in range(len(mutationOrder)):
			n += 1
			if mutations[mutationOrder[i]][1] != '-': 
				currSeq = currSeq[0:mutations[mutationOrder[i]][0]] + mutations[mutationOrder[i]][1] + currSeq[mutations[mutationOrder[i]][0] + 1:]  
			else: 
				currSeq = currSeq[0:mutations[mutationOrder[i]][0]] + currSeq[mutations[mutationOrder[i]][0] + 1:]
				for j in range(mutationOrder[i] + 1, len(mutations)):
					mutations[j] = [mutations[j][0] - 1, mutations[j][1]]
			seqid = "seq{n} p={p} d={d} pos={loc}\n".format(n =n, p= path, d = i, loc=mutations[mutationOrder[i]][0])
			out_d[seqid] = currSeq
		if currSeq != hit.replace('-', ''):
			print("Something went wrong\n",query,'\n',hit,'\n',currSeq,hit.replace('-', ''))
			return
	Write_DB(out_d, output_path)

def Summarize_Paths(edit_fasta_path, rdp_out_edit_path, output_file):
	infasta = open(edit_fasta_path, "r")
	df = Parse_RDP_Output(rdp_out_edit_path)
	df = df[['genus', 'p-genus']].reset_index()
	taxa = dict(zip(df['Seq-ID'].tolist(), zip(df['genus'].tolist(), df['p-genus'].tolist())))
	paths = []
	for line in infasta:
		ma = re.match(r"^>(\S+) p=(\d+) d=\d+ pos=(\d+)", line)
		if ma:
			try: paths[int(ma.group(2)) - 1].append([ma.group(1), ma.group(3)])
			except IndexError: paths.append([[ma.group(1), ma.group(3)]])
	out_f = open(output_file, 'w')
	for i in range(len(paths)): # for each path
		line1 = str(i + 1)
		line2 = " "
		line3 = " "
		for j in range(len(paths[i])): # for each "edit"
			line1 = line1 + "\t" + paths[i][j][1]  # position of edit
			line2 = line2 + "\t" + taxa[paths[i][j][0]][0] # tax label
			line3 = line3 + "\t" + str(taxa[paths[i][j][0]][1]) # confidence
		out_f.write(line1 + "\n")
		out_f.write(line2 + "\n")
		out_f.write(line3 + "\n")
	out_f.close()

def Get_Worst_Case_Adversarial_Sequences(edit_summary_file):
	lines = open(edit_summary_file).readlines()
	j_max = 5000
	max_header = ()
	d = {}
	
	for i in range(0, len(lines), 3):
		positions = lines[i].replace("\n","").split("\t")
		classes = lines[i+1].replace("\n","").split("\t")[1:]
		s1, s2 = classes[0], classes[-1]
		classes = [classes[0]]+classes+[classes[-1]]

		for j in range(0, len(classes)-1):
			if classes[j] != classes[j+1]: break
		if i == 0: seq_id = j+1
		else: seq_id = int(i/3)*(len(positions)-1)+j+1

		if j < j_max:
			j_max = j
			da = j-1
			db = len(classes)-da
			max_header=((int(i/3)+1), da, db, len(classes)-2)
	
	return (s1, s2) + max_header

def Sample_Adversarial_Examples(edit_fasta_file, worst_case_sequence, sample_multiple_seq, seq_count, all):
	fasta_seqs = Load_Database(edit_fasta_file, keep_header = True)
	head, tail = split(edit_fasta_file)
	k = tail.replace(".fasta","")

	min_d = worst_case_sequence[3]
	p = worst_case_sequence[2]
	sa, sb = worst_case_sequence[0], worst_case_sequence[1]

	fasta_ids = []
	sampled_seqs = {}
	for m in fasta_seqs.keys():
		ed = int(m.split(' ')[2].replace("d=",""))
		if 'p='+str(p) in m and ed < min_d: 
			fasta_ids.append(m)
		elif ed == min_d and 'p='+str(p) in m: 
			sampled_seqs[k+'-'+m.replace(" ","_")+"_delta=0"+'_'+sa+'-'+sb] = fasta_seqs[m]
	
	if (sample_multiple_seq == True) or (all == True):
		if all == True: c = len(fasta_ids)
		elif sample_multiple_seq == True: c = min(seq_count, len(fasta_ids))
		rand_ids = np.random.randint(0, len(fasta_ids), c)
		for i in rand_ids:
			seq_id = fasta_ids[i]
			try:
				delta = min_d-int(seq_id.split(" ")[2].replace("d=",""))
				sampled_seqs[k+'-'+seq_id.replace(" ","_")+"_delta="+str(delta)+'_'+sa+'-'+sb] = fasta_seqs[seq_id]
			except KeyError:
				print(seq_id, " not found...")
	return sampled_seqs

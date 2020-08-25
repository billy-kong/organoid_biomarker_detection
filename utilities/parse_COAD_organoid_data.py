##### parse organoid data ###
import numpy as np
import scipy.stats as stat
from collections import defaultdict
import time, os
execfile('pathway_utilities.py', globals())
gene2uniprot, uniprot2gene = geneID2uniprot(), uniprot2geneID()


#
## drug response
def return_COAD_organoid_drug_response_IC50():
	'''
	output = { sample : { drug : IC50 } }
	output_drugList = [list of drugs tested ]
	---
	median IC50 values from di/triplicates are returned
	'''
	output = {} # { sample : { drug : IC50 } }
	drugList = set() # [drugs]

	tmp = {} # { sample : { drug : [list of IC50]}}
	fi_directory = '../data/organoid_COAD/drug_response'
	f = open('%s/supp_table_S2b.txt' %fi_directory, 'r')
	for line in f.xreadlines():
		line = line.strip().split('\t')
		if not 'Organoid' in line[0]:
			if len(line)>5:
				sample, drug, IC50 = line[0].upper(), line[2].upper(), float(line[4])
				if not sample in tmp:
					tmp[sample] = defaultdict(list)
				tmp[sample][drug].append(IC50)
				drugList.add(drug)
	f.close()

	for sample in tmp:
		output[sample] = {}
		for drug in tmp[sample]:
			output[sample][drug] = np.median(tmp[sample][drug])			
	return output, list(drugList)


#
## expression
def parse_GPL16686():
	'''
	output = { probe ID : gene ID }
	'''
	output = {} # { probe ID : gene ID }
	
	print 'importing probe ID - refseq ID, ', time.ctime()
	tmp = {} # { probe ID : RefSeq ID }
	fi_directory = '../data/organoid_COAD/expression/GPL16686'
	f = open('%s/GPL16686_family.soft' %fi_directory, 'r')
	for line in f.xreadlines():
		line = line.strip().split('\t')
		if (not '!' in line[0]) and (len(line)>5):
			probeID, RefSeqID = line[0], line[5]
			if len(RefSeqID) > 0:
				tmp[probeID] = RefSeqID
	f.close()

	# RefSeq ID to gene ID
	print 'importing refseq ID - gene ID, ', time.ctime()
	refseq_gene = {}
	f = open('%s/geneID_RefSeqID_biomart_20190916.txt' %fi_directory, 'r')
	for line in f.xreadlines():
		line = line.strip().split('\t')
		if ('ENSG' in line[0]) and (len(line)==3):
			geneID, refseqID = line[1], line[2]
			refseq_gene[refseqID] = geneID
	f.close()

	print 'output dictionary, ', time.ctime()
	for pID in tmp:
		rID = tmp[pID]
		if rID in refseq_gene:
			gene = refseq_gene[rID]
			if gene in gene2uniprot:
				gID = uniprot2gene[gene2uniprot[gene]]
				output[pID] = gID

	return output


def return_COAD_2015_cell_organoid_RMA_normalized_expression():
	'''
	Returns median gene expression
	-------------------------------------
	output = { sample : { gene : exp } }
	output2 = { sample : { uniprot : exp } }
	'''
	output = {} # { sample : { gene : exp } }
	output2 = {} # { sample : { uniprot : exp } }
	tmpExp = {} # { sample : { gene : [ list of expressions ] } }
	geneList = set()

	fi_directory = '../data/organoid_organoid/expression/GSE64392/GSE64392_series_matrix.txt'
	
	fiList = os.listdir(fi_directory)

	# make expression file
	if not 'geneID_expression_median.txt' in fiList:
		probe_gene = parse_GPL16686() # { probe ID : gene }
		f = open('%s/GSE64392_series_matrix.txt' %fi_directory, 'r')
		for line in f.xreadlines():
			line = line.strip().split('\t')
			if '!Sample_title' in line[0]:
				sampleList = line[1:]
			elif '!Sample_source_name_ch1' in line[0]:
				typeList = line[1:] # normal/tumor info
			elif (not '!' in line[0]) and (not 'ID_REF' in line[0]) and (not line == ''):
				probeID, expList = line[0], line[1:]
				for index, exp in enumerate(expList):
					exp = float(exp)
					sample, sampleType = sampleList[index], typeList[index]
					sample = sample.replace('"','').replace('t', '').upper()

					if (probeID in probe_gene) and ('colon carcinoma organoids' in sampleType):
						geneID = probe_gene[probeID]
						if not sample in tmpExp:
							tmpExp[sample] = defaultdict(list)
						tmpExp[sample][geneID].append(exp)
		f.close()

		for sample in tmpExp:
			output[sample] = {}
			output2[sample] = {}
			for gene in tmpExp[sample]:
				median = np.median(tmpExp[sample][gene])
				output[sample][gene] = median

				uniprot = gene2uniprot[gene]
				output2[sample][uniprot] = median
				geneList.add(gene)
		geneList = list(geneList)

		# make gene expression file
		fo = open('%s/geneID_expression_median.txt'%fi_directory,'w')
		sampleList = output.keys()
		sampleList.sort()
		print >> fo, '\t'.join(['geneID', 'uniprotID']) + '\t' + '\t'.join(map(str, sampleList))
		for gene in geneList:
			uniprot = gene2uniprot[gene]
			tmp_output = [gene, uniprot]
			for sample in sampleList:
				tmp_output.append(output[sample][gene])
			print >> fo, '\t'.join(map(str, tmp_output))
		fo.close()
	
	else:
		f = open('%s/geneID_expression_median.txt'%fi_directory, 'r')
		for line in f.xreadlines():
			line = line.strip().split('\t')
			if 'geneID' in line[0]:
				sampleList = line[2:]
			else:
				gene, uniprot, expList = line[0], line[1], line[2:]
				for index, sample in enumerate(sampleList):
					if not sample in output:
						output[sample] = {}
					if not sample in output2:
						output2[sample] = {}
					output[sample][gene] = float(expList[index])
					output2[sample][uniprot] = float(expList[index])
		f.close()
	return output, output2
				

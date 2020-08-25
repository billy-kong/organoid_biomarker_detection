#### parse patient expression data ###
import pandas as pd
from collections import defaultdict
import gzip
import os, time, json
import scipy.stats as stat
import numpy as np
execfile('pathway_utilities.py', globals())
ensg2gene, gene2uniprot, uniprot2gene = ensembl2geneID(), geneID2uniprot(), uniprot2geneID()

# expression
def parse_TCGA_log2_FPKM(cancer_type):
	'''
	output = { pat : { gene : log2(FPKM+1) } }
	output2 = { pat : { uniprot : log2(FPKM+1) } }
	'''
	output, output2 = {}, {} # { pat : { gene : log2(FPKM+1) } }, { pat : { uniprot : log2(FPKM+1) } }
	if 'coad' in cancer_type.lower():
		output, output2 = parse_TCGA_COAD_log2_FPKM_expression()
	return output, output2


## ===================================================================================
## COAD
def parse_TCGA_COAD_caseid_fileName_barcodeid():
	fi_directory = '../data/TCGA_COAD'#
	output = {} # { case ID : { 'barcode', 'file name' } }
	fileName_barcode = {} # { 'file name' : 'barcode' }

	# caseID and barcodeID
	f = open('%s/clinical.tsv'%fi_directory, 'r')
	for line in f.xreadlines():
		line = line.strip().split('\t')
		if not 'case_id' in line[0]:
			case_id, barcode_id = line[0], line[1]
			if not case_id in output:
				output[case_id] = {}
			output[case_id]['barcode'] = barcode_id
	f.close()

	# caseID and fileName
	with open('%s/files.20190910_COAD_FPKM_UQ.json'%fi_directory) as json_file:
		json_data = json.load(json_file)
		for i in range(len(json_data)):
			case_id = json_data[i]['cases'][0]['case_id']
			file_name = json_data[i]['file_name']
			if case_id in output:
				output[case_id]['file_name'] = file_name

	# fileName_barcode
	for case_id in output:
		if ('barcode' in output[case_id]) and ('file_name' in output[case_id]):
			barcode, fileName = output[case_id]['barcode'], output[case_id]['file_name']
			fileName_barcode[fileName] = barcode

	if not 'COAD_caseID_fileName_barcodeID.txt' in os.listdir(fi_directory):
		fo = open('%s/COAD_caseID_fileName_barcodeID.txt' %fi_directory, 'w')
		print >> fo, '\t'.join(['caseID', 'fileName', 'barcodeID'])
		for caseID in output:
			if ('barcode' in output[case_id]) and ('file_name' in output[case_id]):	
				tmp = [caseID, output[caseID]['file_name'], output[caseID]['barcode']]
				print >> fo, '\t'.join(map(str, tmp))
		fo.close()

	return output, fileName_barcode


def parse_TCGA_COAD_log2_FPKM_expression():
	"""
	returns,
	{ pat : { gene : exp } }, { pat : { uniprot : exp } }

	"""
	fi_directory = '../data/TCGA_COAD'#
	output, output_uniprot = {}, {} # { pat : { gene : exp } }, { pat : { uniprot : exp } }
	_, fileName_barcode = parse_TCGA_COAD_caseid_fileName_barcodeid()
	
	f = open('%s/expression_LOG2_FPKM_UQ.txt' %fi_directory, 'r')
	for line in f.xreadlines():
		line = line.strip().split('\t')
		if 'gene' in line[0]:
			patList = line[2:]
		else:
			gene, uniprot, expList = line[0], line[1], line[2:]
			for index, pat in enumerate(patList):
				exp = float(expList[index])
				if not pat in output:
					output[pat] = {}
				if not pat in output_uniprot:
					output_uniprot[pat] = {}
				output[pat][gene] = exp
				output_uniprot[pat][uniprot] = exp
	f.close()
	return output, output_uniprot


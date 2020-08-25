############# parse data of preclinical models ##############
import os, time
import scipy.stats as stat
import numpy as np
from collections import defaultdict
import pandas as pd
execfile('pathway_utilities.py', globals())
execfile('parse_Drugbank.py', globals())
gene2uniprot, uniprot2gene = geneID2uniprot(), uniprot2geneID()
_, _, _, common_syn, syn_common = parse_Drugbank_drugbankID_synonyms()


# ssGSEA NES data
def parse_organoid_ssGSEA_NES(cancer_type, pathway_source):
	'''
	output = { sample : { pathway : NES } } 
	----
	Input
	cancer_type = 'coad', 'blca'
	pathway_source = 'reactome'
	'''
	if cancer_type.lower() == 'coad':
		execfile('parse_COAD_ssGSEA.py', globals())
		output = return_COAD_ssGSEA_NES( 'organoid', pathway_source )
	return output

## transcriptome data
def parse_organoid_transcriptome(cancer_type):
	'''
	output = { sample ID : { gene : exp } }
	output2 = { sample ID : { uniprot : exp } } 
	---------------------
	cancer_type = 'coad'
	'''
	if cancer_type.lower() == 'coad':
		output, output2 = parse_coad_organoid_transcriptome()
	return output, output2
	
# ----------------------------------------------------------------------------------------------	
## organoid drug response
def parse_organoid_drug_response_commonDrugID( cancer_type, response_unit ):
	'''
	output = { sample : { drug : drug response } } 
	drugList = [ drugs ]
	---
	response_unit == 'IC50', 'AUC'
	'''
	output, drugList = {}, []
	if cancer_type.lower() == 'coad':
		tmp_output, tmp_drugList = parse_coad_organoid_drug_response(response_unit)
	
	for sample in tmp_output:
		output[sample.upper()] = {}
		for drug in tmp_output[sample]:
			cDrug = drug
			if drug in syn_common:
				cDrug = syn_common[drug]
			output[sample.upper()][cDrug] = tmp_output[sample][drug]
			drugList.append(cDrug)
					
	return output, list(set(drugList))


## colorectal cancer organoid
def parse_coad_organoid_transcriptome():
	'''
	returns median gene expression
	------------------------------------------
	output = { sample : { gene : exp } }
	output2 = { sample : { uniprot : exp } }
	'''
	current_dir = os.getcwd()
	os.chdir('/home/junghokong/PROJECT/colon_cancer/code/1_drugResponsePrediction')
	execfile('parse_COAD_organoid_data.py', globals())
	output, output2 = return_COAD_2015_cell_organoid_RMA_normalized_expression()
	os.chdir(current_dir)
	return output, output2


def parse_coad_organoid_drug_response( response_unit ):
	'''
	output = { sample : { drug : drug response } }
	drugList = [ drugs ]
	---
	response_unit = 'IC50'
	---
	median IC50 values from di/triplicates are returned
	'''
	current_dir = os.getcwd()
	os.chdir('/home/junghokong/PROJECT/colon_cancer/code/1_drugResponsePrediction')
	execfile('parse_COAD_organoid_data.py', globals())
	if response_unit == 'IC50':
		output, drugList = return_COAD_organoid_drug_response_IC50()
	os.chdir(current_dir)
	return output, drugList

################ Parse Drug-Target Relationship Provided From < Drugbank > ################
import pandas as pd
from collections import defaultdict
import time, os, csv
execfile('pathway_utilities.py', globals())


def parse_Drugbank_Drug_Target_relationship():
    """
    Drug-target relationship provided by Drugbank
	output = { drug : [ list of targets in gene ID ] }
    """
    print 'all drug-drug targets (including candidate drug targets) from Drugbank data are imported'

    output = defaultdict(list) # { drug : [ list of targets ] }
    gene2uniprot, uniprot2gene = geneID2uniprot(), uniprot2geneID()

    # DRUGBANK Drug ID & Drug Name matching
    annoDic = {} # { drugbank ID : drugName }
    #fi_directory = '/home/junghokong/PROJECT/bladder_cancer/code/15_network_based_analysis/data/Drugbank/drugbank_all_drugbank_vocabulary_downloadDate_20181116.csv/drugbank_vocabulary.csv'
    f = open('../data/drugbank_vocabulary.csv', 'r')
    rdr = csv.reader(f)
    for line in rdr:
        if not 'DrugBank' in line[0]:
            Drugbank_ID, drugName = line[0], line[2].upper()
            annoDic[Drugbank_ID] = drugName
    f.close()

    # DRUGBANK Drug name & Drug Target (uniprot ID)
    tmpDic = defaultdict(list) # { Drug Name : [ list of drug targets ( in gene ID ) ] }
    #fi_directory = '/home/junghokong/PROJECT/bladder_cancer/code/15_network_based_analysis/data/Drugbank/drugbank_all_target_polypeptide_ids_downloadDate_20181113.csv/all.csv'
    f = open('../data/all.csv', 'r')
    rdr = csv.reader(f)
    for line in rdr:
        if not 'ID' == line[0]:
            uniprot_id, drug_id_list = line[5], line[12].split('; ')
            if uniprot_id in uniprot2gene:
                for drug_id in drug_id_list:
                    if drug_id in annoDic:
                        geneID = uniprot2gene[uniprot_id]
                        drugName = annoDic[drug_id].upper()
                        if not geneID in tmpDic[drugName]:
                            tmpDic[drugName].append(geneID)
    f.close()
    return tmpDic


def parse_Drugbank_Drug_Target_relationship_uniprot_commonDrugID():
	'''
    Drug-target relationship provided by Drugbank
	output = { common drug ID : [ list of targets in uniprot ID ] }
	'''
	output = defaultdict(list) # { common drug ID : [ list of targets in uniprot ID ] }
	gene2uniprot, uniprot2gene = geneID2uniprot(), uniprot2geneID()
	commonID_dbID, dbID_drug, drug_dbID, commonID_drug, drug_commonID = parse_Drugbank_drugbankID_synonyms() # dbID : Drugbank ID
	dt = parse_Drugbank_Drug_Target_relationship() # { drug : [ list of targets in gene ID ] }
	
	for drug in dt:
		common_drugID = drug
		if drug in drug_commonID:
			common_drugID = drug_commonID[drug]
		if not common_drugID in output:
			output[common_drugID] = []
		for gene in dt[drug]:
			if gene in gene2uniprot:
				uniprot = gene2uniprot[gene]
				output[common_drugID].append(uniprot)
	return output




def parse_Drugbank_drugbankID_synonyms():
	"""
	Drugbank ID and drug synonyms
	output1 = { common drug name : Drugbank ID }
	output2 = { Drugbank ID : [ list of synonyms and common drug name ] }
	output3 = { common drug name/drug synonym : Drugbank ID }
	output4 = { common drug name : [ list of synonyms ] }
	output5 = { synonym : common drug name }
	"""

    # DRUGBANK Drug ID & Drug Name matching
	output1 = {} # { common drug name : drugbank ID }
	annoDic = defaultdict(list) # { drugbank ID : [drugs] }
	output3 = {} # { drug : drugbank ID }
	output4 = defaultdict(list) # { common drug name : [ list of synonyms ] }
	output5 = {} # { synonym ; common drug name }

	#fi_directory = '/home/junghokong/PROJECT/bladder_cancer/code/15_network_based_analysis/data/Drugbank/drugbank_all_drugbank_vocabulary_downloadDate_20181116.csv/drugbank_vocabulary.csv'
	f = open('../data/drugbank_vocabulary.csv', 'r')
	rdr = csv.reader(f)
	for line in rdr:
		if not 'DrugBank' in line[0]:
			Drugbank_ID, drugName, syn_list = line[0], line[2].upper(), line[5].split(' | ')
			output1[drugName] = Drugbank_ID
			output3[drugName] = Drugbank_ID
			if not drugName in annoDic[Drugbank_ID]:
				annoDic[Drugbank_ID].append(drugName)
			
			for syn in syn_list:
				syn_drug = syn.upper()
				if len(syn_drug)>0:
					if not syn_drug in annoDic[Drugbank_ID]:
						annoDic[Drugbank_ID].append(syn_drug)
					if not syn_drug in output4[drugName]:
						output4[drugName].append(syn_drug)
					if not syn_drug.replace('-','') in output4[drugName]:
						output4[drugName].append(syn_drug.replace('-',''))
					if not syn_drug.replace(' ','') in output4[drugName]:
						output4[drugName].append(syn_drug.replace(' ',''))
					output5[syn_drug] = drugName
					output5[syn_drug.replace('-','')] = drugName
					output5[syn_drug.replace(' ','')] = drugName
					
					output3[syn_drug] = Drugbank_ID
	f.close()
	return output1, annoDic, output3, output4, output5
	
	

import pandas as pd
import os, time
from collections import defaultdict

def parse_ssGSEA_NES(cancer_type, source, pathway_source):
	'''
	returns ssGSEA NES values
	{ sample : { pathway : NES } }
	'''
	if 'coad' in cancer_type.lower():
		output = return_COAD_ssGSEA_NES( source, pathway_source )
	return output




## COAD
def return_COAD_ssGSEA_NES( source, pathway_source ):
	'''
	returns ssGSEA NES values
	{ sample : { pathway : NES } }
	'''
	output = {} # { sample : { pathway : NES } }
	fi_dir = '../python/results/COAD'
	folderList = os.listdir(fi_dir)
	if source in folderList:
		if os.path.isdir('%s/%s/%s_ssgsea_result' %(fi_dir, source, pathway_source.lower())) == True:
			f = open('%s/%s/%s_ssgsea_result/gseapy.samples.normalized.es.txt' %(fi_dir, source, pathway_source.lower()),'r')
			for line in f.xreadlines():
				line = line.strip().split('\t')
				if 'Term' in line[0]:
					sampleList = line[1:]
				if ( not '#' in line[0] ) and ( not 'Term' in line[0] ):
					pathway, NES_list = line[0], map(float, line[1:])
					for index, sample in enumerate(sampleList):
						if not sample in output:
							output[sample] = {}
						output[sample][pathway] = NES_list[index]
			f.close()
	return output				

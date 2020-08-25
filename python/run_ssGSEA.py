############### run single-sample GSEA (ssGSEA) for cancer patient or organoid data ###############
import gseapy
import gseapy as gp
import scipy.stats as stat
import numpy as np
import time, os
import pandas as pd
from collections import defaultdict
cur_dir = os.getcwd()
os.chdir('../utilities')
execfile('pathway_utilities.py', globals())
execfile('parse_preclinical_data.py', globals())
#execfile('parse_COAD_organoid_data.py', globals())
execfile('parse_patient_expression.py', globals())
execfile('parse_GDSC_cell_line_data.py', globals())
gene2uniprot, uniprot2gene = geneID2uniprot(), uniprot2geneID()
os.chdir(cur_dir)

## INITIALIZE
#======================
# INITIALIZE PARAMETERS
source = 'TCGA' # 'organoid', 'TCGA', 'GDSC'
cancer_type = 'COAD'
testing_pathway_list = ['REACTOME']

#==================
# IMPORT EXPRESSION
print 'importing expression for %s, ' %source, time.ctime()

expDic = {} # { sample ID : { gene in uniprot : exp } }
expDic_geneID = {} # { sample ID : { gene : exp } }
geneList, sampleList = [], []

if source.lower() == 'organoid':
	expDic_geneID, expDic = parse_organoid_transcriptome(cancer_type)

if source.upper() == 'TCGA':
	expDic_geneID, expDic = parse_TCGA_log2_FPKM(cancer_type)

sampleList = expDic.keys()
sampleList.sort()

for sample in expDic_geneID:
	if len(geneList) == 0:
		geneList = expDic_geneID[sample].keys()
	else:
		geneList = list(set(geneList).intersection(expDic_geneID[sample].keys()))




#========================================
# IMPORT PATHWAYS FOR ENRICHMENT ANALYSIS
print 'importing pathways, ', time.ctime()
reactome = reactome_genes_uniprot() # { pathway : [ gene list ] }
pathwayDic = {'reactome':reactome} # 


## PRINT ssGSEA RESULTS

#===============
# MAKE DIRECTORY
fo_directory = './results' 
dir_list = [cancer_type.upper(), source]
for d in dir_list:
	if os.path.isdir('%s/%s' %(fo_directory, d)) == False:
		os.mkdir('%s/%s' %(fo_directory, d))
	fo_directory = '%s/%s'%(fo_directory, d)

#=====================
# MAKE GSEA INPUT FILE
fiList = os.listdir(fo_directory)

# gene expression
if not 'expression.txt' in fiList:
    fo = open('%s/expression.txt' %(fo_directory), 'w')
    print >> fo, '\t'.join(['NAME', 'DESCRIPTION']) + '\t' + '\t'.join(sampleList)
    for gene in geneList:
        tmp = [gene, 'na']
        for sample in sampleList:
            tmp.append(expDic_geneID[sample][gene])
        print >> fo, '\t'.join(map(str, tmp))
    fo.close()


#=======
# ssGSEA
for testing_pathway in testing_pathway_list:
    if testing_pathway.lower() in pathwayDic:
        print 'running ssGSEA for %s ... , ' %testing_pathway.lower(), time.ctime()
        # gene sets for ssGSEA
        gene_sets = {}
        pw_list = []
        for pw in pathwayDic[testing_pathway.lower()]:
            for uniprot in pathwayDic[testing_pathway.lower()][pw]:
                if uniprot in uniprot2gene:
                    gene = uniprot2gene[uniprot]
                    if gene in geneList:
                        if not pw in gene_sets:
                            gene_sets[pw] = []
                        gene_sets[pw].append(gene)
        pw_list = gene_sets.keys()

        fo = open('%s/%s.gmt' %(fo_directory, testing_pathway.lower()), 'w')
        for pw in pw_list:
            print >> fo, pw + '\t' + '\t'.join(gene_sets[pw])
        fo.close()


        # ssGSEA
        ss = gp.ssgsea(data='%s/expression.txt'%(fo_directory), outdir='%s/%s_ssgsea_result'%(fo_directory, testing_pathway.lower()),
                       gene_sets='%s/%s.gmt' %(fo_directory, testing_pathway.lower()),
                       sample_norm_method='rank', permutation_num=0, no_plot=True, scale=True, min_size=2)


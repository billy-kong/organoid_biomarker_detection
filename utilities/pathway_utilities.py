## Providing pathway (reactome) genes, gene-UniProtID mapping, gene-Ensembl ID mapping
import pandas as pd
import numpy as np
import scipy.stats as stat
from collections import defaultdict
import os, time


# REACTOME genes
def reactome_genes(): # provide in a dictionary
    output = defaultdict(list)
    output_list = []
    f = open('../data/msigdb.v6.1.symbols.gmt.txt','r')
    for line in f.xreadlines():
        line = line.strip().split('\t')
        if 'REACTOME' in line[0]:
            reactome = line[0]
            output_list.append(reactome)
            for i in range(2, len(line)):
                gene = line[i]
                output[reactome].append(gene)
    f.close()
    return output

def reactome_genes_uniprot():
    output, reactome = defaultdict(list), reactome_genes()
    gene2uniprot = geneID2uniprot()
    for pathway in reactome:
        for gene in reactome[pathway]:
            if gene in gene2uniprot:
                uniprot = gene2uniprot[gene]
                if not uniprot in output[pathway]:
                    output[pathway].append(uniprot)
    return output


    
    
## gene annotation conversion utilities
def convert_geneList_to_uniprotList( input_geneList ):
    output = []
    for gene in input_geneList:
        if gene in gene2uniprot:
            output.append(gene2uniprot[gene])
    return list(set(output))

def convert_uniprotList_to_geneList( input_uniprotList ):
    output = []
    for uniprot in input_uniprotList:
        if uniprot in uniprot2gene:
            output.append(uniprot2gene[uniprot])
    return list(set(output))
    

## gene annotation    
# ensembl gene annotation
def annotation():
    geneID2ensembl, ensembl2geneID = defaultdict(set), {}
    df = pd.read_csv('../data/2017_07_31_biomart_protein_coding_genes.txt', sep='\t')
    for i in range(len(df)):
        geneID, ensembl = df['Gene name'][i], df['Gene stable ID'][i]
        geneID2ensembl[ geneID ].add( ensembl )
        ensembl2geneID[ ensembl ] = geneID

    for geneID in geneID2ensembl:
        geneID2ensembl[geneID] = list(geneID2ensembl[geneID])
    return geneID2ensembl, ensembl2geneID

def ensembl2geneID():
    output = {} # { ensembl : geneID }
    df = pd.read_csv('../data/2017_07_31_biomart_protein_coding_genes.txt', sep='\t')
    for i in range(len(df)):
        ensembl, gene = df['Gene stable ID'][i], df['Gene name'][i]
        output[ensembl] = gene
    return output

def geneID2uniprot():
    output = {} # { gene ID : uniprot ID }
    df = pd.read_csv('../data/uniprot_homoSapiens_multipleGeneName_20180802.tab', sep='\t')
    for i in range(len(df)):
        uniprot, geneList = df['Entry'][i], df['Gene names'][i]
        if pd.isnull(geneList) == False:
            geneList = geneList.split()
            for gene in geneList:
                output[gene] = uniprot
    return output

def uniprot2geneID():
    output = {} # { uniprot ID : gene ID }
    df = pd.read_csv('../data/uniprot_homoSapiens_multipleGeneName_20180802.tab', sep='\t')
    for i in range(len(df)):
        uniprot, geneList = df['Entry'][i], df['Gene names'][i]
        if pd.isnull(geneList) == False:
            geneList = geneList.split()
            gene = geneList[0]
            output[uniprot] = gene
    return output
gene2uniprot, uniprot2gene = geneID2uniprot(), uniprot2geneID()

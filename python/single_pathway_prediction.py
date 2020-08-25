## Return predictive performance of a pathway
import time, os, math, random
from collections import defaultdict
import pandas as pd
import numpy as np
import scipy.stats as stat
from sklearn.preprocessing import StandardScaler # zscore standardization
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import RidgeCV, Ridge
from sklearn.svm import SVR
import matplotlib.pyplot as plt
import lifelines; from lifelines.statistics import logrank_test; from lifelines import KaplanMeierFitter

## Initialize
test_types = {'coad':'FLUOROURACIL'}
ML_algorithm_list = ['Ridge', 'SVR', 'LinearRegression']
pathway_source = 'reactome'
n_jobs = 10
network = 'STRING_700'
zscore_cutoff = -1.2816
testing_pathway_rank = 1
tmp_dir = './results/single_pathway_predictions'

start_time = time.ctime()

## Parse proximal pathways
def return_proximal_pathways(network, pathway_source, zscore_cutoff):
	'''
	output = { drug : [ proximal pathways ] }
	'''
	output = defaultdict(list)
	if pathway_source.lower() == 'reactome':
		df = pd.read_csv('../data/coad_blca_organoid_drugs_zscore_result_reactome.txt', sep='\t')
		tmp_drugs = df.columns
		for tmp_drug in tmp_drugs:
			for drug in tmp_drug.split('_'):
				pathway_list = list(df.loc[df[tmp_drug]<=zscore_cutoff, :]['Pathway'])
				output[drug] = pathway_list
	return output


## Output file
fo_all = open('%s/single_pathway_predictions.txt'%(tmp_dir), 'w')
print >> fo_all, '\t'.join(['cancer_type', 'drug', 'ML', 'pathway_rank', 'pathway', '5yr_responder', '5yr_nonresponder', 'pvalue'])

fo_all_coef = open('%s/single_pathway_coefficients.txt'%(tmp_dir), 'w')
print >> fo_all_coef, '\t'.join(['cancer', 'drug','num_samples','ML', 'pathway', 'reg_coef', 'abs_reg_coef'])


## Analyze
for ML in ML_algorithm_list:
	# set directory
	fo_dir = tmp_dir
	if os.path.isdir('%s/%s'%(fo_dir, ML)) == False:
		os.mkdir('%s/%s'%(fo_dir, ML))
	fo_dir = '%s/%s'%(fo_dir, ML)

	# cancer specific analysis
	for cancer_type in test_types:
		drug = test_types[cancer_type]
		print '\n\n-----------------\ntesting for %s / %s / %s,  '%(cancer_type, drug, ML), time.ctime()


		## Import data
		c_dir = os.getcwd()
		os.chdir('../utilities')
		execfile('parse_patient_clinical_data.py', globals())
		execfile('parse_preclinical_model_data.py', globals())
		execfile('parse_ssGSEA.py', globals())
		nes_dic = parse_ssGSEA_NES(cancer_type, 'organoid', pathway_source) 
		response_dic, drugList = parse_organoid_drug_response_commonDrugID(cancer_type, 'IC50')
		network_dic = return_proximal_pathways(network, pathway_source, zscore_cutoff) # { drug : [ pathways ] }

		drugPat, patDrug = parse_TCGA_drug_treatment_data(cancer_type)
		surDic = parse_TCGA_survival_data_boolean_format(cancer_type)
		patNES = parse_ssGSEA_NES(cancer_type, 'TCGA', pathway_source)
		
		if pathway_source.lower() == 'reactome':
			feature_list = reactome_genes_uniprot()
			feature_list = feature_list.keys()
		os.chdir(c_dir)


		## feature_list // proximal pathways only
		for sample in nes_dic:
			feature_list = list(set(feature_list) & set(nes_dic[sample].keys()))
		for pat in list(set(surDic.keys())&set(patNES.keys())):
			feature_list = list(set(feature_list)&set(patNES[pat].keys()))
		feature_list = list(set(feature_list) & set(network_dic[drug]))

		## Predict drug response
		# scale expressions (organoid)
		expList, samples, responses = [], [], []
		for sample in list(set(nes_dic.keys())&set(response_dic.keys())):
			if drug in response_dic[sample]:
				samples.append(sample)#; samples = sorted(samples, reverse=True)
				responses.append(response_dic[sample][drug])
		for sample in samples:
			tmp = []
			for feature in feature_list:
				if feature in network_dic[drug]: # proximal pathways only
					tmp.append(nes_dic[sample][feature])
			expList.append(tmp)
		scaler = StandardScaler()
		scaler.fit(expList)
		scaled_expList = scaler.transform(expList) # scaled expression
		scaled_expList = np.array(scaled_expList)

		# regression (organoid)
		if ML == 'Ridge':
			regr = RidgeCV(cv=3, alphas=np.arange(0.1,1,0.1)).fit(scaled_expList, responses)
		if ML == 'SVR':
			regr = SVR(kernel='linear').fit(scaled_expList, responses)
		if ML == 'LinearRegression':
			regr = LinearRegression().fit(scaled_expList, responses)
		feature_importance = list(regr.coef_)
		if ML == 'SVR':
			feature_importance = feature_importance[0]

		# feature ranks
		coef_dic = {} # { feature : coefficient }
		abs_coef_dic = {}
		for feature, coef in zip(feature_list, feature_importance):
			# print regression coefficients
			print >> fo_all_coef, '\t'.join(map(str, [cancer_type, drug, len(samples), ML, feature, coef, np.abs(coef)]))
			# coefficents
			coef_dic[feature] = coef
			abs_coef_dic[feature] = np.abs(coef)
		r = {key: rank for rank, key in enumerate(sorted(set(abs_coef_dic.values()), reverse=True), 1)}
		feature_rank_dic = {k: r[v] for k,v in abs_coef_dic.items()}
		

		# scale expressions (patient)
		pat_expDic = {} # { pat : { feature : scaled expression } }
		pat_expList, pat_samples = [], []
		for pat in list(set(surDic.keys())&set(patNES.keys())):
			pat_samples.append(pat)
			tmp = []
			for feature in feature_list:
				tmp.append(patNES[pat][feature])
			pat_expList.append(tmp)
		scaler = StandardScaler()
		scaler.fit(pat_expList)
		scaled_pat_expList = scaler.transform(pat_expList)

		for p_index, pat in enumerate(pat_samples):
			pat_expDic[pat] = {}
			for f_index, feature in enumerate(feature_list):
				pat_expDic[pat][feature] = scaled_pat_expList[p_index][f_index]
		
		# single-pathway prediction
		features_used = []; coef_used = []
		for fi_index, fi in enumerate(feature_importance):
			feature = feature_list[fi_index]
			if testing_pathway_rank == feature_rank_dic[feature]:
				features_used.append(feature); coef_used.append(fi)

			
		# predicted drug response (patient)
		pred_response = {} # { pat : predicted drug response }
		month_dic, status_dic = defaultdict(list), defaultdict(list)
		fiveYear_dic = {} # { predicted response : 5 year survival }

		for pat in list(set(pat_expDic.keys())&set(drugPat[drug])&set(surDic.keys())):
			pred_r = 0
			for feature, coef in zip(features_used, coef_used):
				pred_r += coef * pat_expDic[pat][feature]
			pred_response[pat] = pred_r

		# classify patients
		response_cutoff = np.median(pred_response.values())
		for pat in pred_response:
			if pred_response[pat] <= response_cutoff:
				cls = 'Responder'
			else:
				cls = 'Nonresponder'
			month_dic[cls].append(surDic[pat]['months'])
			status_dic[cls].append(surDic[pat]['status'])
		
		# logrank Test
		results = logrank_test(month_dic['Responder'], month_dic['Nonresponder'], event_observed_A=status_dic['Responder'], event_observed_B=status_dic['Nonresponder'])
		pvalue = results.p_value
		
		for cls in month_dic:
			kmf = KaplanMeierFitter()
			kmf.fit(month_dic[cls], status_dic[cls])
			fiveYear_dic[cls] = kmf.predict(60)

		# draw survival plot
		f = plt.figure(figsize=(4,4))
		ax = f.add_subplot(1,1,1)
		plt.title('%s / %s / %s / %s\npvalue=%.4f\n'%(cancer_type, drug, ML, testing_pathway_rank, pvalue), fontsize=8)
		
		c1 = KaplanMeierFitter()
		ax = c1.fit(month_dic['Responder'], status_dic['Responder'], label='Responder (n=%s)'%len(month_dic['Responder'])).plot(ax=ax, ci_show=True, c='r')
		

		c2 = KaplanMeierFitter()
		ax = c1.fit(month_dic['Nonresponder'], status_dic['Nonresponder'], label='Nonresponder (n=%s)'%len(month_dic['Nonresponder'])).plot(ax=ax, ci_show=True, c='b')
		
		plt.xlabel('Survival (months)')
		plt.ylabel('Percent survival')
		ymin, ymax = 0, 1.1
		plt.ylim(ymin, ymax)
		plt.plot([60, 60], [ymin, ymax], c='k', linestyle='--')
		plt.tight_layout()
		plt.savefig('%s/%s_%s_rank_%s.jpg'%(fo_dir, cancer_type, drug, testing_pathway_rank), format='jpg')
		plt.savefig('%s/%s_%s_rank_%s.eps'%(fo_dir, cancer_type, drug, testing_pathway_rank), format='eps', dpi=300)
		plt.close()
		
		
fo_all.close()		
fo_all_coef.close()
print 'process complete, start time: %s - end time: %s ' %(start_time, time.ctime())

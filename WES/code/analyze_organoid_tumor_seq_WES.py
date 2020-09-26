############### analyze organoid/primary tumor mutations ###########
import openpyxl
import pandas as pd
from collections import defaultdict
import os, time
c_dir = os.getcwd()
execfile('./annotation_utilities.py', globals())
uniprot2gene, gene2uniprot = uniprot2geneID(), geneID2uniprot()

## 
# import WES variant calling data
mut_dic = {} # { sampleName : { gene : [ variant_class ] } }

fi_directory = '../data'
xlsxFile = '%s/Mutation_list.xlsx'%fi_directory 

sheetList = []
wb = openpyxl.load_workbook(xlsxFile) # bring sheet names using openpyxl
for i in wb.get_sheet_names():
	sheetList.append(i)

xlsx = pd.ExcelFile(xlsxFile)

for sampleName in sheetList:
	df = pd.read_excel(xlsx, sampleName)
	if not sampleName in mut_dic:
		mut_dic[sampleName] = defaultdict(list)
	for i in range(len(df)):
		uniprotID, ensemblID, geneID, variant_class, consequence, biotype, SIFT, polyphen, exon  = df['SWISSPROT'][i], df['Gene'][i], df['SYMBOL'][i], df['VARIANT_CLASS'][i], df['Consequence'][i], df['BIOTYPE'][i], df['SIFT'][i], df['PolyPhen'][i], df['EXON'][i]
		if (biotype == 'protein_coding'):# and ('deleterious' in SIFT) and ('damaging' in polyphen) and ('-' != exon):
			if uniprotID in uniprot2gene:
				geneID = uniprot2gene[uniprotID]
				if not variant_class in mut_dic[sampleName][geneID]:
					mut_dic[sampleName][geneID].append(variant_class)

# sample list
sampleList = sorted(mut_dic.keys())

##
# compare mutations between organoid and primary tumor

# mutation count
mutCount_dic = {} # { gene : mut count }

for sample in sampleList:
	for gene in mut_dic[sample]:
		if not gene in mutCount_dic:
			mutCount_dic[gene] = 0
		mutCount_dic[gene] += 1

geneList = sorted(mutCount_dic, key=mutCount_dic.get, reverse=True)


# mutation summary file
fo_directory = '../result'
fo = open('%s/mutation_variant_summary_allSamples.txt'%fo_directory, 'w')

tmp_sample = sampleList
print >> fo, 'gene' + '\t' + '\t'.join(tmp_sample)

for gene in geneList:
	tmp = [gene]
	for sample in sampleList:
		if gene in mut_dic[sample]:
			vcList = '_'.join(mut_dic[sample][gene])
			tmp.append(vcList)
		else:
			tmp.append('-')
	print >> fo, '\t'.join(map(str, tmp))
fo.close()

# mutation summary file / numeric
fo = open('%s/mutation_variant_summary_numeric_allSamples.txt'%fo_directory, 'w')
print >> fo, '\t'.join(map(str, ['#', '1:SNV', '2:INSERTION', '-1:DELETION', '0:NULL']))
print >> fo, 'gene' + '\t' + '\t'.join(tmp_sample)

for gene in geneList:
	tmp = [gene]
	for s in sampleList:
		if gene in mut_dic[s]:
			tmpList = []
			for mut in mut_dic[s][gene]:
				if 'SNV' == mut:
					tmpList.append(1)
				elif 'insertion' == mut:
					tmpList.append(2)
				elif 'deletion' == mut:
					tmpList.append(-1)
			vcList = '_'.join(map(str, tmpList))
			tmp.append(vcList)
		else:
			tmp.append('0')
	print >> fo, '\t'.join(map(str, tmp))
fo.close()


# mutation concordance
fo = open('%s/mutation_concordance_variant_level.txt'%fo_directory,'w')
print >> fo, '\t'.join(['Sample', 'Parental', 'Concordant', 'Parental tumor only', 'Organoid only'])

for sample in sampleList:
	if not 'parental' in sample:
		sample_type = sample.split(' ')[0]
		organoid, primary = sample, '%s parental tumour'%sample_type
		output = [sample, primary]
		tmp_dic = defaultdict(list)
		
		for s in [organoid, primary]:
			for gene in mut_dic[s]:
				for variant_class in mut_dic[s][gene]:
					tmp_dic[s].append('%s_%s'%(gene,variant_class))
		
		common = len(set(tmp_dic[organoid]).intersection(tmp_dic[primary]))
		total = len(set(tmp_dic[organoid]).union(tmp_dic[primary]))
		organoid_specific = len(tmp_dic[organoid])-common
		primary_specific = len(tmp_dic[primary])-common
		
		output.append( float(common)/float(total) )
		output.append( float(primary_specific)/float(total) )
		output.append( float(organoid_specific)/float(total) )
		print >> fo, '\t'.join(map(str, output))

fo.close()


fo = open('%s/mutation_concordance_gene_level.txt'%fo_directory,'w')
print >> fo, '\t'.join(['Sample', 'Parental', 'Concordant', 'Parental tumor only', 'Organoid only'])

for sample in sampleList:
        
	if not 'parental' in sample:
		sample_type = sample.split(' ')[0]
		organoid, primary = sample, '%s parental tumour'%sample_type
		output = [sample, primary]
		tmp_dic = defaultdict(list)
		
		for s in [organoid, primary]:
			for gene in mut_dic[s]:
				tmp_dic[s].append(gene)
		
		common = len(set(tmp_dic[organoid]).intersection(tmp_dic[primary]))
		total = len(set(tmp_dic[organoid]).union(tmp_dic[primary]))
		organoid_specific = len(tmp_dic[organoid])-common
		primary_specific = len(tmp_dic[primary])-common
		
		output.append( float(common)/float(total) )
		output.append( float(primary_specific)/float(total) )
		output.append( float(organoid_specific)/float(total) )
		print >> fo, '\t'.join(map(str, output))

fo.close()


# cross mutation concordance - variant level
all_sampleList = sampleList

fo = open('%s/cross_mutation_concordance_variant_level_allSamples.txt' %fo_directory, 'w')
print >> fo, '\t'.join(['Sample1', 'Sample2', 'Sample1_Sample2', 'Concordant', 'Sample1 only', 'Sample2 only'])

for s1_index, s1 in enumerate(all_sampleList):
	for s2_index, s2 in enumerate(all_sampleList):
		if s1_index<s2_index:
			output = [s1, s2, '%s_%s'%(s1,s2)]
			tmp_dic = defaultdict(list)

			for s in [s1, s2]:
				for gene in mut_dic[s]:
					for variant_class in mut_dic[s][gene]:
						tmp_dic[s].append('%s_%s'%(gene,variant_class))
			
			common = len(set(tmp_dic[s1]).intersection(tmp_dic[s2]))
			total = len(set(tmp_dic[s1]).union(tmp_dic[s2]))
			s1_specific = len(tmp_dic[s1])-common
			s2_specific = len(tmp_dic[s2])-common
					
			output.append( float(common)/float(total) )
			output.append( float(s1_specific)/float(total) )
			output.append( float(s2_specific)/float(total) )
			print >> fo, '\t'.join(map(str, output))

fo.close()


# cross mutation concordance - gene level
all_sampleList = sampleList

fo = open('%s/cross_mutation_concordance_gene_level_allSamples.txt' %fo_directory, 'w')
print >> fo, '\t'.join(['Sample1', 'Sample2', 'Sample1_Sample2', 'Concordant', 'Sample1 only', 'Sample2 only'])

for s1_index, s1 in enumerate(all_sampleList):
	for s2_index, s2 in enumerate(all_sampleList):
		if s1_index<s2_index:
			output = [s1, s2, '%s_%s'%(s1,s2)]
			tmp_dic = defaultdict(list)

			for s in [s1, s2]:
				for gene in mut_dic[s]:
					tmp_dic[s].append(gene)
			
			common = len(set(tmp_dic[s1]).intersection(tmp_dic[s2]))
			total = len(set(tmp_dic[s1]).union(tmp_dic[s2]))
			s1_specific = len(tmp_dic[s1])-common
			s2_specific = len(tmp_dic[s2])-common
					
			output.append( float(common)/float(total) )
			output.append( float(s1_specific)/float(total) )
			output.append( float(s2_specific)/float(total) )
			print >> fo, '\t'.join(map(str, output))

fo.close()

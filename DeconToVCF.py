import pandas as pd
import os
from copy import deepcopy
import csv
import datetime
import argparse


"""
Python script to convert DeCon output to VCF style template for further analysis.
Run script and give pathway to decon output folder (-d), ped file (-p), and output file (-o)
Author: Kalon Grimes

"""


"""
# set variables for GT

code copy:
if float(ratio) < low_gt:
	gt_list.append('0/0')
elif low_gt <= float(ratio) < mid_gt:
	gt_list.append('0/1')
elif mid_gt <= float(ratio) < high_gt:
	gt_list.append('1/1')
else:
	gt_list.append('0/1')
"""
low_gt = 0.4
mid_gt = 0.7
high_gt = 1.3


"""
Set variables for CN

code copy:
if float(rr) < low_cn:
	cn = '0'
elif low_cn <= float(rr) < mid_cn:
	cn = '1'
elif mid_cn <= float(rr) < high_cn:
	cn = '2'
else:
	cn = '3'
"""
low_cn = 0.4
mid_cn = 0.7
high_cn = 1.3



def get_CNV_table(raw_data_path):
	"""
	function to get CNVs reported seperately in the raw_data files, and put them into one
	dataframe which is filtered to take one GenomicsID per sampleID. If more than one CNV exists
	for a GenomicID region then the table selects the highest read quality (BF value)
	gt variables are defined at start of main script
	input: pathway to decon output raw_data folder
	output: filtered dataframe of GenomicsID-SampleID
	"""

	# create blank dataframe ready to append all the raw_data files into
	full_df_headers = ['CNV.ID', 'Sample', 'Correlation', 'N.comp', 'Start.b', 'End.b', 'CNV.type', 'N.exons', 'Start', 'End', 'Chromosome', 'Genomic.ID', 'BF', 'Reads.expected', 'Reads.observed', 'Reads.ratio', 'Gene', 'Custom.first', 'Custom.last']
	full_df = pd.DataFrame(columns = full_df_headers)

	# iterate through all raw_data files and append them to full_df
	file_list = os.listdir(raw_data_path)

	for file in file_list:
		if file[-7:] == 'all.txt':
			df1 = pd.read_csv(os.path.join(raw_data_path,file), sep = '\t')
			full_df = full_df.append(df1, ignore_index = True)

	# dedup across Sample, genomicID, CNV.type, and then take highest BF.
	full_df.sort_values(by = ['Genomic.ID', 'Sample', 'CNV.type', 'BF'], ascending = [True, True, True, True])
	full_df.drop_duplicates(subset = ['Genomic.ID', 'Sample', 'CNV.type'], inplace = True, ignore_index = True, keep = 'last')

	# seperate Sample column into two columns sampleID and runID
	sampleid_list = []
	runID_list = []
	for sample in full_df['Sample']:
		sampleid = sample.split('_')[4]
		runIDlist = sample.split('_')[0:4]
		runID = '_'.join(runIDlist)
		sampleid_list.append(sampleid)
		runID_list.append(runID)
	full_df['sampleid'] = sampleid_list
	full_df['runID'] = runID_list

	# add genotype column based off of reads ratio values (to be decided in discussion)
	gt_list = []
	for ratio in full_df['Reads.ratio']:
		if float(ratio) < low_gt:
			gt_list.append('0/0')
		elif low_gt <= float(ratio) < mid_gt:
			gt_list.append('0/1')
		elif mid_gt <= float(ratio) < high_gt:
			gt_list.append('1/1')
		else:
			gt_list.append('0/1')

	full_df['Genotype'] = gt_list

	return full_df



def get_sampleIDs(pedfile_path):
	"""
	function to get a list of sampleIDs from the PED file
	input: pathway to pedfile
	output: list of sampleIDs on the run
	"""

	# read ped file into dataframe and then convert column 2 into a list of sampleIDs
	sampleID_df = pd.read_csv(pedfile_path, sep = "\t", names = [1,2,3,4,5,6])
	sampleID_list = sampleID_df[2].tolist()

	# remove NTC value
	sampleID_list.remove('NTC')

	return sampleID_list



def get_vcf_dict(CNV_table, sampleID_list):

	"""
	function to create a nested dictionary of all the information for each genomic ID.
	format:
	{gen ID : {	sampleID-1 : GT:BF:CN:RE:RO:RR,
				sampleID-2 : GT:BF:CN:RE:RO:RR,
				etc...,
				meta : "{chrom}~{pos}~{ref}~{alt}~{qual}~{filt}~{info}~{form}" }
	}
	cn variables are defined at start of main script
	input: CNV_table and sampleID_list
	return: vcf_dict with all information to create csv file
	"""

	# create blank dicts
	vcf_dict = {}
	sample_dict = {}

	# initiate rownum counter
	rownum = 0

	# iterate through df and sample list to create dict template
	for sampleID in sampleID_list:
		sample_dict[sampleID] = ''
	# add meta option to sample_dict to add info later
	sample_dict['meta'] = ''

	for gen_id in CNV_table['Genomic.ID']:
		# print(gen_id)
		# print(rownum)

		# check if the gen id value is already a key in the dictionary
		if gen_id in vcf_dict.keys():
			# print("in dict already")
			for sample in vcf_dict[gen_id]:
				# print(sample)

				# see if the sample listed in the CNV table is the one linked to the genomicID
				if sample == (CNV_table.at[rownum, 'sampleid']):
					# print("sample matched")
					gt = CNV_table.at[rownum, 'Genotype']
					bf = CNV_table.at[rownum, 'BF']
					re = CNV_table.at[rownum, 'Reads.expected']
					ro = CNV_table.at[rownum, 'Reads.observed']
					rr = CNV_table.at[rownum, 'Reads.ratio']
					if float(rr) < low_cn:
						cn = '0'
					elif low_cn <= float(rr) < mid_cn:
						cn = '1'
					elif mid_cn <= float(rr) < high_cn:
						cn = '2'
					else:
						cn = '3'
					vcf_dict[gen_id][sample] = f'{gt}:{bf}:{cn}:{re}:{ro}:{rr}'

				else:
					# print("sample not matched")
					if vcf_dict[gen_id][sample] == '':
						vcf_dict[gen_id][sample] = './.:.:.:.:.'
		else:
			# new gen_id so need to add the sample_dict to it
			# print("sample not in dict, adding now")
			vcf_dict[gen_id] = deepcopy(sample_dict)
			for sample in vcf_dict[gen_id]:
				# print(sample)

				# see if the sample is meta, then populate with meta that easily splittable
				if sample == 'meta':
					# print("sample is meta")
					gen_id_nochr = gen_id[3:]
					chrom = CNV_table.at[rownum, 'Chromosome']
					pos = CNV_table.at[rownum, 'Start']
					end = CNV_table.at[rownum, 'End']
					ref = 'N'
					region_length = int(end) - int(pos)

					if CNV_table.at[rownum, 'CNV.type'] == 'deletion':
						alt = '<DEL>'
						ID = f'LOSS:{gen_id_nochr}'
						info = f'SVLEN=-{region_length};SVTYPE=CNV;END={end};REFLEN={region_length}'
					elif CNV_table.at[rownum, 'CNV.type'] == 'duplication':
						alt = '<DUP>'
						ID = f'GAIN:{gen_id_nochr}'
						info = f'SVLEN={region_length};SVTYPE=CNV;END={end};REFLEN={region_length}'
					else: 
						alt = '<UNK>'
						ID = f'UNK:{gen_id_nochr}'
						info = f'SVLEN={region_length};SVTYPE=CNV;END={end};REFLEN={region_length}'

					qual = '.'
					filt = 'PASS'
					form = 'GT:BF:CN:RE:RO:RR'
					vcf_dict[gen_id][sample] = f'{chrom},{pos},{ID},{ref},{alt},{qual},{filt},{info},{form}'
					# print(vcf_dict[gen_id][sample])

				# see if the sample listed in the CNV table is the one linked to the genomicID	
				elif sample == (CNV_table.at[rownum, 'sampleid']):
					# print("sample matched")
					gt = CNV_table.at[rownum, 'Genotype']
					bf = CNV_table.at[rownum, 'BF']
					re = CNV_table.at[rownum, 'Reads.expected']
					ro = CNV_table.at[rownum, 'Reads.observed']
					rr = CNV_table.at[rownum, 'Reads.ratio']
					if float(rr) < low_cn:
						cn = '2'
					elif low_cn <= float(rr) < mid_cn:
						cn = '2'
					elif mid_cn <= float(rr) < high_cn:
						cn = '2'
					else:
						cn = '3'
					vcf_dict[gen_id][sample] = f'{gt}:{bf}:{cn}:{re}:{ro}:{rr}'

				else:
					# print("sample not matched")
					if vcf_dict[gen_id][sample] == '':
						vcf_dict[gen_id][sample] = './.:.:.:.:.'
		rownum += 1

	return vcf_dict



def get_export_list(vcf_dict,sampleID_list):
	"""
	function to convert the vcf_dict nested format into a readable/exportable list of lines
	input: vcf_dict, list of samples
	output: list of lines in vcf format for each genomicID and all samples
	(samples are . if not got variant)
	"""

	# create empty export list
	export_list = []

	# create and export column header string
	sampleID_string = ','.join(sampleID_list)
	column_headers = f'#CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,{sampleID_string}'
	#print(column_headers)
	export_list.append(column_headers)


	for key in vcf_dict:

		# reset line variable
		line = ''

		# get meta information first
		line = vcf_dict[key]['meta']

		# get sample genotype field using specific sampleID callers not positionals (dict don't hold orders well)
		for sample in sampleID_list:
			genotype = vcf_dict[key][sample]
			line += f',{genotype}'

		# append line to export_list
		export_list.append(line)

	# convert csv to tab delim
	tab_export_list = []
	for line in export_list:
		line = line.replace(',','\t')
		tab_export_list.append(line)

	return tab_export_list



def get_vcf_header(run):
	"""
	hard-coded vcf header returned as a list, uses runID to populate run info and date fields
	"""

	datestr = run.split('_')[0]
	year = int(f'20{datestr[:2]}')
	month = int(datestr[2:4])
	day = int(datestr[4:])
	date = datetime.date(year, month, day)

	vcf_header = [
	'##fileformat=VCFv4.2',
	'##reference=GRCh37',
	'##contig=<ID=1,length=249250621>',
	'##contig=<ID=2,length=243199373>',
	'##contig=<ID=3,length=198022430>',
	'##contig=<ID=4,length=191154276>',
	'##contig=<ID=5,length=180915260>',
	'##contig=<ID=6,length=171115067>',
	'##contig=<ID=7,length=159138663>',
	'##contig=<ID=8,length=146364022>',
	'##contig=<ID=9,length=141213431>',
	'##contig=<ID=10,length=135534747>',
	'##contig=<ID=11,length=135006516>',
	'##contig=<ID=12,length=133851895>',
	'##contig=<ID=13,length=115169878>',
	'##contig=<ID=14,length=107349540>',
	'##contig=<ID=15,length=102531392>',
	'##contig=<ID=16,length=90354753>',
	'##contig=<ID=17,length=81195210>',
	'##contig=<ID=18,length=78077248>',
	'##contig=<ID=19,length=59128983>',
	'##contig=<ID=20,length=63025520>',
	'##contig=<ID=21,length=48129895>',
	'##contig=<ID=22,length=51304566>',
	'##contig=<ID=X,length=155270560>',
	'##contig=<ID=Y,length=59373566>',
	'##ALT=<ID=DEL,Description="Deletion relative to the reference">',
	'##ALT=<ID=DUP,Description="Region of elevated copy number relative to the reference">',
	'##INFO=<ID=REFLEN,Number=1,Type=Integer,Description="Number of REF positions included in this record">',
	'##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
	'##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
	'##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
	'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
	'##FORMAT=<ID=BF,Number=1,Type=Float,Description="quality BF value">',
	'##FORMAT=<ID=CN,Number=1,Type=String,Description="Copy number estimate">',
	'##FORMAT=<ID=RE,Number=1,Type=Integer,Description="Number of reads expected">',
	'##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Number of reads observed">',
	'##FORMAT=<ID=RR,Number=1,Type=Float,Description="Reads Ratio">'
	]

	return vcf_header





 #####################################
 ##			PROGRAMME CODE 			##
 #####################################

if __name__ == '__main__':

 	# args
	parser = argparse.ArgumentParser()
	parser.add_argument('--rawdata','-d', help = 'pathway to decon output directory eg: post-processing/results/cnv_svs/raw_data/')
	parser.add_argument('--pedfile','-p', help = 'pathway to pedfile eg: post-processing/results/ped/<runid>.ped')
	parser.add_argument('--outfile','-o', help = 'filename for output eg: <runid>_decon.vcf')
	args = parser.parse_args()

	# get table of all CNVs from DeCon output files.
	CNV_table = get_CNV_table(args.rawdata)
	# print(CNV_table)

	# get list of sampleIDs from PED files
	sampleID_list = get_sampleIDs(args.pedfile)
	# print(sampleID_list)

	# get vcf dict format of all data
	vcf_dict = get_vcf_dict(CNV_table, sampleID_list)
	#print(vcf_dict)

	# get export table
	export_list = get_export_list(vcf_dict,sampleID_list)

	# get runid information
	runID = CNV_table.at[1,'runID']
	#print(runID)

	# get hard-coded VCF header list
	vcf_header = get_vcf_header(runID)

	# export information to file
	with open(args.outfile,'w',newline='') as file:
		for line in vcf_header:
			file.write(line + '\n')
		for line in export_list:
			file.write(line + '\n')

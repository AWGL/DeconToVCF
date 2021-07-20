import unittest
import pandas as pd
import os
from copy import deepcopy
import csv
import datetime
from DeconToVCF import get_vcf_dict, get_CNV_table, get_sampleIDs

"""
test usage:

run test script from within a directory which contains the "testdata" folder.
This folder must also contain:
raw_data directory with multiple decon output files
.ped file with sampleIDs relevant to the decon output files

"""


class TestDeconToVCF(unittest.TestCase):

	def test_get_sampleIDs(self):
		sample_IDs = get_sampleIDs('testdata/test.ped')
		self.assertEqual(sample_IDs, ['16M04452', '17M16605', '17M17936', '17M18406', '18M00728', '18M01315', '18M03204', '18M04957', '18M05431', '18M07042', '18M09699', '18M13937', '18M14187', '18M19396', '18m11396', 'NTC'])

	def test_get_CNV_table(self):
		CNV_table = get_CNV_table('testdata/raw_data', False)

		CNV_table = CNV_table[CNV_table['Sample'] == '190705_M00766_0238_000000000-CD3NW_16M04452']

		self.assertEqual(CNV_table.at[1,'Genomic.ID'], 'chr19:11221328-11221447<DEL>')

	def test_get_vcf_dict_simple_case(self):
		CNV_table = get_CNV_table('testdata/raw_data', False)
		sample_IDs = get_sampleIDs('testdata/test.ped')
		vcf_dict = get_vcf_dict(CNV_table, sample_IDs)
		self.assertEqual(vcf_dict['chr19:11226770-11227674<DEL>']['18M13937'].split(':')[0], '0/1')


	def test_get_vcf_dict_dupdel_case(self):
		# edge case test to show same genomic location can be both Dup and Del
		CNV_table = get_CNV_table('testdata/raw_data', False)
		sample_IDs = get_sampleIDs('testdata/test.ped')
		vcf_dict = get_vcf_dict(CNV_table, sample_IDs)
		self.assertIn('chr19:11210899-11227674<DEL>',vcf_dict.keys())
		self.assertIn('chr19:11210899-11227674<DUP>',vcf_dict.keys())

	def test_get_vcf_dict_complex(self):
		CNV_table = get_CNV_table('testdata/raw_data2', False)
		sample_IDs = get_sampleIDs('testdata/test.ped')
		vcf_dict = get_vcf_dict(CNV_table, sample_IDs)
		self.assertEqual(vcf_dict['chr19:11210899-11227674<DEL>']['18M09699'].split(':')[0], '0/1')

		self.assertEqual(vcf_dict['chr19:11226770-11227674<DEL>']['18M13937'].split(':')[0], '0/1')

		self.assertEqual(vcf_dict['chr19:11226770-11227674<DEL>']['17M16605'].split(':')[0], './.')

		self.assertEqual(vcf_dict['chr19:11226870-11226970<DEL>']['18M13937'].split(':')[0], '1/1')

		self.assertEqual(vcf_dict['chr19:11226870-11226970<DEL>']['17M18406'].split(':')[0], '1/1')

		self.assertEqual(vcf_dict['chrX:11226870-11246970<DUP>']['17M18406'].split(':')[0], '0/1')

if __name__ == '__main__':
	unittest.main()
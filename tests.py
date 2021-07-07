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
		self.assertEqual(sample_IDs, ['18M13937','17M17936','16M04452','18M04957','18M01315','18m11396','18M03204','18M07042','17M18406','18M05431','17M16605','18M00728','18M09699','18M19396','18M14187'])

	def test_get_CNV_table(self):
		CNV_table = get_CNV_table('testdata/raw_data')
		self.assertEqual(CNV_table.at[1,'Genomic.ID'], 'chr19:11210899-11227674<DEL>')

	def test_get_vcf_dict_simple_case(self):
		CNV_table = get_CNV_table('testdata/raw_data')
		sample_IDs = get_sampleIDs('testdata/test.ped')
		vcf_dict = get_vcf_dict(CNV_table, sample_IDs)
		self.assertEqual(vcf_dict['chr19:11226770-11227674<DEL>']['18M13937'], '0/1:19.7:2:5757:3970:0.69')


	def test_get_vcf_dict_dupdel_case(self):
		# edge case test to show same genomic location can be both Dup and Del
		CNV_table = get_CNV_table('testdata/raw_data')
		sample_IDs = get_sampleIDs('testdata/test.ped')
		vcf_dict = get_vcf_dict(CNV_table, sample_IDs)
		self.assertIn('chr19:11210899-11227674<DEL>',vcf_dict.keys())
		self.assertIn('chr19:11210899-11227674<DUP>',vcf_dict.keys())

if __name__ == '__main__':
	unittest.main()
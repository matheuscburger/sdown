#!/usr/bin/env python
# -*- coding: utf8 -*-

__author__ = 'Matheus Carvalho Bürger'
__license__ = 'agpl'

import xlrd
import argparse
import sys
import re

def changeT(x):
	return re.sub(u'\u0442', 'T', x)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Get a tsv file from a xls file")
	parser.add_argument("-o", type=argparse.FileType('w'), default= sys.stdout, help="output file")
	parser.add_argument("-l", type=argparse.FileType('w'), default= sys.stderr, help="log file")
	args = parser.parse_args()
	out =  args.o
	log =  args.l
	filename = "/data/projects/bit_data/down/htqpcr/mRNA.imune/raw/dadosInflamacaoCompletos_thresholdAjustado.xlsx"
	myfile = xlrd.open_workbook(filename)
	mysheet = myfile.sheet_by_index(0)
	header = mysheet.row_values(0)
	intCols = ('Biological Group Name', 'Sample Name', 'Target Name', u'C\u0442')	# colunas de interesse
	idxIntCols = map(header.index, intCols)	# indice das colunas de interesse
	print >>out, "\t".join(map(changeT, intCols))
	for i in xrange(1, mysheet.nrows):
		if sum(map(lambda x: x in intCols,mysheet.row_values(i))) > 0:
			print >>log, "Header in the middle of file, line ", i
			continue
		try:
			print >>out, "\t".join([str(mysheet.row_values(i, j, j+1).pop()) for j in idxIntCols])
		except Exception, e:
			print >>log, e
			print >>log, "Linha ", i
			break



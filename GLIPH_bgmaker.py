#! /usr/bin/python
import sys
import os
from string import split
import fileinput
from TCRseq_import_func import tcr_import
import argparse

# WHAT TYPE OF 
parser = argparse.ArgumentParser(description='This parser provides options for importing the tabulated output of a TCR repertoire dataset as rendered by the IPOP pipeline, or by a number of other commercial platforms for upstream handling of the raw sequencing data (see options that follow). Generates a summary report on number of nucleotide CDR3 sequences, amino acid CDR3 motifs, V-J cassette combinations, aa/VJ clonotypes, statistics on the distributions of amino acid CDR3 per VJ and VJ per amino acid CDR3. Also generates tabulated output of the repertoire according to these identifiers. Generates a directory using the basename of the file to place these outputs. See below for options regarding the inclusion of \"productive\" vs \"non-productive\" sequences, and parsing of unresolvable V- and J-cassettes.', epilog='Version 1.0 written in Python2.7.5. DEPENDENCIES: Functions in TCRseq_H_JSD_func.py and TCRseq_import_func.py.')
parser.add_argument('-f','--infile',action='append',help='Input files. Tab-delimited text file, windows formatted, *.txt or *.tsv.')
parser.add_argument('-t','--filetype',default='IPOP',choices=['IPOP','IREP','IMSQ2','IMSQ3','AUTO','CUSTOM','PS'],help='This program is designed to handle several styles of tabulated input. All must be tab-delimited text formatted for windows. The \"AUTO\" option will attempt to autodetect based on the column names. The \"CUSTOM\" option will allow user input to determine key column identities, so brace yourself for a lot of questions.  Also, do not run the \"CUSTOM\" option in the background (by putting \"&\" after your command); this angers bash. Using the \"CUSTOM\" option will eventually make it possible to perform downstram analysis given incomplete information (such as amino-acid only), but the summary output has not yet been debugged.')
parser.add_argument('-o','--outfile',help='Name of output file.')
args = parser.parse_args()

infiles = args.infile

# output of tcr_import: dCDR3,dVJ,dNT,dTOT,dAAperVJ,dVJperAA,ct,num_unprod_aa,reads_unprod_aa,num_notrans_aa,reads_notrans_aa,num_unres_vj,num_filt_nt,reads_filt_nt
a=0
b=0
c=0

dTOTDATA = {}

output = open(args.outfile,'w')

for infile in infiles:
	instrip = str(infile.rstrip('.txt''.tsv''.productive'))

	repdata = tcr_import(infile,args.filetype,a,b,c)			# dCDR3,dVJ,dNT,dTOT,dAAperVJ,dVJperAA,ct,len(dUNPROD.keys()),sum(dUNPROD.values()),len(dNOTRANS.keys()),sum(dNOTRANS.values()),len(dUNRES.keys()),sum(dUNRES.values()),len(dNTFILT.keys()),sum(dNTFILT.values())

	#dCDR3 = repdata[0]
	#dVJ = repdata[1]
	#dNT = repdata[2]
	dTOT = repdata[3]

	#dAAperVJ = repdata[4]
	#dVJperAA = repdata[5]
	ct = repdata[6]
	#num_unprod_aa = repdata[7]
	#reads_unprod_aa = repdata[8]
	#num_notrans_aa = repdata[9]
	#reads_notrans_aa = repdata[10]
	#num_unres_vj = repdata[11]
	#reads_unres_vj = repdata[12]
	#num_filt_nt = repdata[13]
	#reads_filt_nt = repdata[14]
	print instrip
	print ct,'total reads in'
	#print len(dCDR3.keys()),'aaCDR3'
	#print len(dVJ.keys()),'VJ'
	#print len(dNT.keys()),'ntCDR3'
	print len(dTOT.keys()),'totCDR3'

	for totid in dTOT.keys():
		aa = totid[0]
		combo = split(str(totid[1]),'.')
		vcass = combo[0]
		jcass = combo[1]
		reads = int(dTOT[totid])
		output.write('>%(vcass)s,%(jcass)s,%(aa)s;%(vcass)s 300 0;;%(jcass)s 30 0;%(aa)s;;;;;;;;;;;;;;;\n' % vars() )
		# >TRBV6-5,TRBJ2-3,CAAGGGITDTQYF;TRBV6-5 300 0;;TRBJ2-3 30 0;CAAGGGITDTQYF;;;;;;;;;;;;;;;
		output.write('%(aa)s\n' % vars() )
		# CAAGGGITDTQYF
output.close()

print 'done'


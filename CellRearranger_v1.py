#! /usr/bin/python
import sys
import os
from string import split
import fileinput
#from TCRseq_import_func import tcr_import
import argparse


# infile1 = CellRanger clonotypes.csv file, takes the clonotypes and their frequencies
# infile2 = CellRanger consensus_annotations.csv file

sampledir = sys.argv[1]

infile1 = str(sampledir) + '/outs/clonotypes.csv'
infile2 =  str(sampledir) + '/outs/consensus_annotations.csv'
instrip1 = str(infile1.rstrip('.csv'))
instrip1 = str(infile1.rstrip('.csv'))

dFREQ = {}

k1 = 0
for line in fileinput.input([infile1]):
	if k1 != 0:
		llist = split(line,',')
		clonotype = llist[0]
		freq = llist[1]
		dFREQ[clonotype] = freq
	k1 = k1 + 1
fileinput.close()

dTRA = {}			# key=clonotype, value =[TRXinfo1,TRXinfo2,...]; TRXinfo=[aa,vcass,jcass]
dTRB = {}

k2 =0
for line in fileinput.input([infile2]):
	if k2 != 0:
		llist = split(line,',')
		clonotype = llist[0]
		vcass = llist[4]
		jcass = llist[6]
		aa = llist[10]
		info = aa,vcass,jcass
		if llist[3] == 'TRA':		
			try:
				dTRA[clonotype].append(info)
			except KeyError:
				dTRA[clonotype] = [info]
		elif llist[3] == 'TRB':
			try:
				dTRB[clonotype].append(info)
			except KeyError:
				dTRB[clonotype] = [info]
	k2 = k2+1
fileinput.close()

############# WRITE OUTPUT ######################

samplesplit = split(sampledir,'/')
samplename = samplesplit[len(samplesplit)-1]

outfile = str(sampledir) + '/' + str(samplename) + '_forGLIPH.tsv'
output = open(outfile,'w')

h1 = 'CDR3b'
h2 = 'TRBV'
h3 = 'TRBJ'
h4 = 'CDR3a'
h5 = 'TRAV'
h6 = 'TRAJ'
h7 = 'Patient'  # Clonotype
h8 = 'Count'	# Frequency, from CellRanger, which gives whole numbers of cells
output.write('%(h1)s\t%(h2)s\t%(h3)s\t%(h4)s\t%(h5)s\t%(h6)s\t%(h7)s\t%(h8)s\n ' % vars() )

# FOR EACH CLONOTYPE MAKE A LINE FOR EACH CHAIN PERMUTATION AND DISTRIBUTE THE FREQ OF THAT CLONOTYPE OVER THOSE PERMUTATIONS

for clonotype in dFREQ.keys():
	freq = dFREQ[clonotype]
	try:
		tra_list = dTRA[clonotype]
		a = len(tra_list)
	except KeyError:
		tra_list = [['','','']]
		a = 1	
	try:
		trb_list = dTRB[clonotype]
		b = len(trb_list)		
	except KeyError:
		trb_list =[['','','']]
		b = 1
	freq_per = float(freq)/float(a*b)
	for tra in tra_list:
		a_aa = tra[0]
		a_vcass = tra[1]
		a_jcass = tra[2]
		for trb in trb_list:
			b_aa = trb[0]
			b_vcass = trb[1]
			b_jcass = trb[2]
			output.write('%(b_aa)s\t%(b_vcass)s\t%(b_jcass)s\t%(a_aa)s\t%(a_vcass)s\t%(a_jcass)s\t%(clonotype)s\t%(freq_per)s\n' % vars() )
output.close()


#! /usr/bin/python
import numpy as np
import sys
import os
import random
from string import split

# -------- FUNCIONS ---------

def revcomp(seq):						# borrowed this from a blogspot fellow called 'CrazyHotTommy'...
	seq_dict = {'A':'T','T':'A','G':'C','C':'G'}
	return "".join([seq_dict[base] for base in reversed(seq)])

# -------- INPUT --------------

cassettefile = sys.argv[1]
filesource = sys.argv[2]		#cdr3.txt or cdr3.err.txt from PS pipeline

cassettefilesplit = split(cassettefile,'.')
refstrip = cassettefilesplit[0]

runfiles = []

if os.path.isdir(filesource) == True:
	for file in os.listdir(filesource):
		runfiles.append(filesource + '/' + file)
else:
	runfiles.append(filesource)

print runfiles

dirname = str(refstrip) + '_fraghists'
cmd1 = 'mkdir %(dirname)s' % vars()
os.system(cmd1)

# ------------ CASSETTE REFERENCE SEQUENCES -----------

dSEQ = {}								# Dictionary of sequences from the cassette reference file

refinput = open(cassettefile,'r')
k=0
for line in refinput.readlines():
	if k == 1:
		llist = split(line)
		cass = llist[1]			# this is based on my "TR[A,B].mm10.ref.txt" files
		seq = llist[7]
		dSEQ[cass] = seq
	k=1
refinput.close()

#print dSEQ.keys()

print 'Cassette reference loaded!'

# ------------ INPUT FILES ------------------

dFRAG = {}						# accumulate fragment lengths per cassette
dMAPPED = {}					# accumulate read lengths per cassette
dNOT = {}						# accumulate cassettes that weren't in the reference

for infile in runfiles:
	print infile
	
	infilesplit = split(infile,'.')
	instrip = infilesplit[0].strip(str(filesource) + '/')

	input = open(infile,'r')
	k = 0
	m = 0
	n = 0
	for line in input.readlines():
		if k != 0 or k == 0:
			llist = split(line)
			line_id = llist[0]
			if llist[2].find('V') != -1:
				name1 = llist[2]
				read1 = str(llist[7])
				name2 = llist[3]
				read2 = str(llist[8])
				bigread = str(llist[9])
			elif llist[2].find('J') != -1:
				name1 = llist[3]
				read1 = str(llist[8])
				name2 = llist[2]
				read2 = str(llist[7])
				bigread = str(llist[9])
			else:
				print llist[2],llist[3]
				pass
			try:															# cassette that is a V
				seq = dSEQ[name1]
				lenseq = len(seq)
				index1 = int(seq.find(read1[0:20]))			# matched from the 5' of where the mapped region on the read starts
				if index1 > 0:
					fraglen = lenseq-index1					# what fragment of the reference sequence showed up in that read
					readlen = len(read1)					# the size of the cassette-mapped read
					try:
						lenlist = dFRAG[name1]
						lenlist.append(fraglen)
						dFRAG[name1] = lenlist
						seqlist = dMAPPED[name1]
						seqlist.append(readlen)
						dMAPPED[name1] = seqlist
					except KeyError:
						lenlist = []
						lenlist.append(fraglen)
						dFRAG[name1] = lenlist
						seqlist = []
						seqlist.append(readlen)
						dMAPPED[name1] = seqlist
					m=m+1
				elif index1 == 0:
					fraglen = lenseq
					readlen = len(read1)
					try:
						lenlist = dFRAG[name1]
						lenlist.append(fraglen)
						dFRAG[name1] = lenlist
						seqlist = dMAPPED[name1]
						seqlist.append(readlen)
						dMAPPED[name1] = seqlist
					except KeyError:
						lenlist = []
						lenlist.append(fraglen)
						dFRAG[name1] = lenlist
						seqlist = []
						seqlist.append(readlen)
						dMAPPED[name1] = seqlist
					m = m+1
				else:
					n = n+1
					# print '1fail',line_id,name1,index1,read1,bigread
			except KeyError:
				try:
					notlist = dNOT[name1]
					notlist.append(line_id)
					dNOT[name1]=notlist
				except KeyError:
					notlist = []
					notlist.append(line_id)
					dNOT[name1]=notlist
			try:															# cassette that is a J
				seq = dSEQ[name2]
				lenseq = len(seq)
				index2 = int(seq.rfind(read2[len(read2)-20:len(read2)+1]))	# matched from the 3'
				if index2 > 0:
					fraglen = index2
					readlen = len(read2)
					try:
						lenlist = dFRAG[name2]
						lenlist.append(fraglen)
						dFRAG[name2] = lenlist
						seqlist = dMAPPED[name2]
						seqlist.append(readlen)
						dMAPPED[name2] = seqlist
					except KeyError:
						lenlist = []
						lenlist.append(fraglen)
						dFRAG[name2] = lenlist
						seqlist = []
						seqlist.append(readlen)
						dMAPPED[name2] = seqlist
					m = m + 1
				elif index2 == 0:
					fraglen = lenseq
					readlen = len(read2)
					try:
						lenlist = dFRAG[name2]
						lenlist.append(fraglen)
						dFRAG[name2] = lenlist
						seqlist = dMAPPED[name2]
						seqlist.append(readlen)
						dMAPPED[name2] = seqlist
					except KeyError:
						lenlist = []
						lenlist.append(fraglen)
						dFRAG[name2] = lenlist
						seqlist = []
						seqlist.append(readlen)
						dMAPPED[name2] = seqlist
					m = m + 1
				else:
					if read2.rfind(seq[len(seq)-20:len(seq)+1]) != -1:
						fraglen = lenseq
						readlen = len(read2)
						try:
							lenlist = dFRAG[name2]
							lenlist.append(fraglen)
							dFRAG[name2] = lenlist
							seqlist = dMAPPED[name2]
							seqlist.append(readlen)
							dMAPPED[name2] = seqlist
						except KeyError:
							lenlist = []
							lenlist.append(fraglen)
							dFRAG[name2] = lenlist
							seqlist = []
							seqlist.append(readlen)
							dMAPPED[name2] = seqlist
					else:
						# print '2fail',line_id,name2,index2,index2r,read2,bigread
						n = n +1
			except KeyError:
				try:
					notlist = dNOT[name2]
					notlist.append(line_id)
					dNOT[name2]=notlist
				except KeyError:
					notlist = []
					notlist.append(line_id)
					dNOT[name2]=notlist
		k = k+1
	input.close()

	print n,'mapfails'
	print m,'readsmapped'
	print dNOT.keys()

# ---------------------------------------------------------------

dDATA1 = {}
dDATA2 = {}

for cassette in dFRAG.keys():
	fraglengths = dFRAG[cassette]
	num1 = len(fraglengths)
	median1 = np.median(fraglengths)
	mean1 = np.mean(fraglengths)
	stdev1 = np.std(fraglengths)
	dDATA1[cassette] = num1,median1,mean1,stdev1
	print 'ref',cassette,num1,median1,mean1,stdev1

	seqlengths = dMAPPED[cassette]
	num2 = len(seqlengths)
	median2 = np.median(seqlengths)
	mean2 = np.mean(seqlengths)
	stdev2 = np.std(seqlengths)
	dDATA2[cassette] = num2,median2,mean2,stdev2
	print 'read',cassette,num2,median2,mean2,stdev2

	if cassette.find('/') != -1:
		stupidsplit = split(cassette,'/')
		cassettename = str(stupidsplit[0]) + '.' + str(stupidsplit[1])
	else:
		cassettename = cassette

	outfile1 = dirname + '/' + str(cassettename) + '_fraglengths.out'
	output1 = open(outfile1,'w')
	for x in range(0,len(fraglengths)):
		pt1 = fraglengths[x]
		pt2 = seqlengths[x]
		pt3 = pt1-pt2
		output1.write('%(pt1)s\t%(pt2)s\t%(pt3)s\n' % vars())
	output1.close()


outfile2 = dirname + '/' + str(refstrip) + '_fraglengths.summary.txt'

output2 = open(outfile2,'w')
h1 = 'CassetteName'
h2 = 'CassetteLength'
h3 = 'NumberOfReads'
h4 = 'MedianLenMappedRef'
h5 = 'MeanLenMappedRef'
h6 = 'STDdevOfLenMappedRef'
h7 = 'MedianLenMappedRead'
h8 = 'MeanLenMappedRead'
h9 = 'STDdevOfLenMappedRead'
output2.write('%(h1)s\t%(h2)s\t%(h3)s\t%(h4)s\t%(h5)s\t%(h6)s\t%(h7)s\t%(h8)s\t%(h9)s\n' % vars())

for cassette in dDATA1.keys():
	pt1 = cassette
	pt2 = len(dSEQ[cassette])
	datalist = dDATA1[cassette]
	seqlist = dDATA2[cassette]
	pt3 = datalist[0]
	pt4 = datalist[1]
	pt5 = datalist[2]
	pt6 = datalist[3]
	pt7 = seqlist[1]
	pt8 = seqlist[2]
	pt9 = seqlist[3]
	output2.write('%(pt1)s\t%(pt2)s\t%(pt3)s\t%(pt4)s\t%(pt5)s\t%(pt6)s\t%(pt7)s\t%(pt8)s\t%(pt9)s\n' % vars())
output2.close()

print infile,'Done'
#! /usr/bin/python
import numpy as np
import sys
import os
import random
from string import split
import string

# -------- FUNCTIONS ---------

def revcomp(seq):						# borrowed this from a blogspot fellow called 'CrazyHotTommy'...
	seq_dict = {'A':'T','T':'A','G':'C','C':'G'}
	return "".join([seq_dict[base] for base in reversed(seq)])

def hamdist(str1,str2):		#borrowed from activestate
	diffs = 0
	for ch1, ch2 in zip(str1,str2):
		if ch1 != ch2:
			diffs += 1
	return diffs


# -------- INPUT --------------

cassette_ref = sys.argv[1]		# annotated cassette reference file (e.g. TRA.mm10.ref.txt)
length_ref = sys.argv[2]		# *_fraglengths_summary.txt from JS_maplocation.py
fam_name = sys.argv[3]			# curate this yourself
consensus = sys.argv[4]			# curate this by aligning all the family members of fam_name and finding a consensus sequence that is upstream of where you think the reads end

refstrip = cassette_ref[0]

#dirname = str(refstrip) + '_fraghists'
#cmd1 = 'mkdir %(dirname)s' % vars()
#os.system(cmd1)

# ------------ CASSETTE REFERENCE SEQUENCES: INPUT AND TRIM TO START AT CONSENSUS -----------

dREFSEQ = {}

cassnames = []
sequences = []

refinput = open(cassette_ref,'r')
k=0
for line in refinput.readlines():
	llist = split(line)
	if k == 1 and llist[1].startswith(fam_name):
		name = llist[1]
		cassnames.append(name)		# this is based on my "TR[A,B].mm10.ref.txt" files
		seq = llist[7]
		seqtrim = seq[seq.rfind(str(consensus))+len(consensus):len(seq)+1]		# grab the sequence from the 3' end of the consensus to the end of each cassette
		sequences.append(seqtrim)
		dREFSEQ[name] = len(seqtrim)
	k=1
refinput.close()

print len(cassnames),len(sequences)

print dREFSEQ.keys()

dCASSLEN = {}

lengthinput = open(length_ref,'r')
k=0
for line in lengthinput.readlines():
	if k != 0:
		llist = split(line,'\t')
		cassname = llist[0]
		if k == 1 and cassname.startswith(fam_name):
	#		length_diff = float(llist[3]) - float(llist[6])	# this is based on my "TR[A,B].mm10.ref.txt" files
	#		dCASSLEN[cassname] = int(max(length_diff,5))
	#		dCASSLEN[cassname] = int(length_diff)
			dCASSLEN[cassname] = int(round(float(llist[3])))		# median length of the reference sequence that appeared in a read
	k=1
lengthinput.close()

print dCASSLEN.keys()
print dCASSLEN.values()


# ------------ TRIM SEQUENCES ------------------

dHAM = {}

c = 0
for i in range(0,len(sequences)):
	for j in range(0,len(sequences)):
		pair = i,j
		frag1 = sequences[i]
		frag2 = sequences[j]
		try:
#			a = min(len(frag1),len(frag2),len(frag1)-dCASSLEN[cassnames[i]],len(frag2)-dCASSLEN[cassnames[j]])
			pos1 = len(frag1)-dCASSLEN[cassnames[i]]
			pos2 = len(frag2)-dCASSLEN[cassnames[j]]
			a = max(pos1,pos2)
			b = min(dCASSLEN[cassnames[i]],dCASSLEN[cassnames[j]])
#			print cassnames[i],cassnames[j],len(frag1),len(frag2),dCASSLEN[cassnames[i]],dCASSLEN[cassnames[j]],pos1,pos2
#			print frag1[a:a+b]
#			print frag2[a:a+b]
			dHAM[pair] = hamdist(frag1[a:a+b],frag2[a:a+b])
#			dHAM[pair] = hamdist(frag1[len(frag1)-a:len(frag1)],frag2[len(frag2)-a:len(frag2)])
#			print pair,frag1,frag2
#			print dHAM[pair],a,cassnames[i],len(frag1)-a,cassnames[j],len(frag2)-a
		except KeyError:
			dHAM[pair]='N/A'
			pass

for pair in dHAM.keys():
	try:
		if int(dHAM[pair]) <= 3:
			print cassnames[pair[0]],cassnames[pair[1]],dHAM[pair]
	except ValueError:			# If "N/A"
		pass

outfile = str(fam_name) + '_hammingtable.tsv'

output = open(outfile,'w')
h1 = fam_name
output.write('%(h1)s' % vars())
for x in range(0,len(sequences)):
	h2 = cassnames[x]
	output.write('\t%(h2)s' % vars())
output.write('\n' % vars())

for y in range(0,len(sequences)):
	pt1 = cassnames[y]
	output.write('%(pt1)s' % vars())
	for z in range(0,len(sequences)):
		pair = y,z
		pt2 = dHAM[pair]
		output.write('\t%(pt2)s' % vars())
	output.write('\n' % vars())
output.close()

print fam_name,'Done'
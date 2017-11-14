#! /usr/bin/python
import numpy as np
import sys
import os
import random
from string import split

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

cassette_ref = sys.argv[1]
length_ref = sys.argv[2]
fam_name = sys.argv[3]
consensus = sys.argv[4]

refstrip = cassette_ref[0]

#dirname = str(refstrip) + '_fraghists'
#cmd1 = 'mkdir %(dirname)s' % vars()
#os.system(cmd1)

# ------------ CASSETTE REFERENCE SEQUENCES: INPUT AND TRIM TO START AT CONSENSUS -----------

dCASS = {}

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
		sequences.append(seq[seq.rfind(str(consensus))+len(consensus):len(seq)+1])
		dCASS[name] = 5
	k=1
refinput.close()

print len(cassnames),len(sequences)

#print dSEQ.keys()

lengthinput = open(length_ref,'r')
k=0
for line in lengthinput.readlines():
	llist = split(line)
	cassname = llist[0]
	if k == 1 and cassname.startswith(fam_name):
		length_diff = float(llist[3]) - float(llist[6])	# this is based on my "TR[A,B].mm10.ref.txt" files
		dCASS[cassname] = int(max(length_diff,5))
	k=1
lengthinput.close()

print dCASS.keys()
print dCASS.values()


# ------------ TRIM SEQUENCES ------------------

dHAM = {}

c = 0
for i in range(0,len(sequences)):
	for j in range(0,len(sequences)):
		pair = i,j
		frag1 = sequences[i]
		frag2 = sequences[j]
		a = min(len(frag1),len(frag2),len(frag1)-dCASS[cassnames[i]],len(frag2)-dCASS[cassnames[j]])
		dHAM[pair] = hamdist(frag1[0:a],frag2[0:a])
		print dHAM[pair],a,cassnames[i],len(frag1)-a,cassnames[j],len(frag2)-a

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
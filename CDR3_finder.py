#! /bin/python
import os
import sys
import csv
from string import split
import string
import random
import argparse
from itertools import izip
import fileinput


def choose_delimiter(filename):
	d = ''
	if filename[len(filename)-4:len(filename)] == '.tsv': #or filename[len(filename)-4:len(filename)] == '.txt':
		d = '\t'
	if filename[len(filename)-4:len(filename)] == '.csv':
		d = ','
	else:
		d = '\t'
	return d
	print d

def countif(rangelist,value):
	ct = 0
	for pt in rangelist:
		if pt == value:
			ct = ct + 1
	return ct

# ------------------------- INPUT ARGUMENT PARSER -------------------------------

parser = argparse.ArgumentParser(description='This parser provides options for querying a repertoire file (\"rfile\") for whether it contains CDR3s from a query file, and reproducing the repertoire file with information from the query file appended.  Currently set up for the query file to be .csv because that\'s what Cell Ranger produces.', epilog='Version 1')

# QUERY FILE
parser.add_argument('qfile',help='Query file name')
parser.add_argument('qcol',help='Column of query file containing CDR3s (0 indexed).')
# REPERTOIRE FILE
parser.add_argument('rfile',help='Repertoire file name')
parser.add_argument('rcol',help='Column of repertoire file containing CDR3s (0 indexed).')
# OPTIONAL ARGUMENTS
parser.add_argument('-qh','--qheader',action='store_true',default=True,help='Query file has header row.')
parser.add_argument('-qt','--qtrim',default=0,help='Number of residues to trim from the end of the query CDR3 sequences')
parser.add_argument('-rh','--rheader',action='store_true',default=True,help='Repertoire file has header row.')
parser.add_argument('-rt','--rtrim',default=0,help='Number of residues to trim from the end of the repertoire CDR3 sequences')
parser.add_argument('-c','--classes',action='append',help='If there is a way CDR3 types will be designated within the query file (i.e. TRA, TRB), and you only want to query some of them, or you just want to have them annotated as such in the output file, flag \"-c\" for each class.')
parser.add_argument('-dt','--tcrdelim',default=';',help='Delimiter between items if there are multiple CDR3s in a cell of the table. For now default is \";\" because that is what Cell Ranger uses.')
parser.add_argument('-dc','--classdelim',default=':',help='Delimiter between the class and the item if there are multiple CDR3s in a cell of the table. For now default is \":\" because that is what Cell Ranger uses.')
parser.add_argument('-xcol','--extracolumn',action='append',help='Are there other columns of data you would like to have associated with the hits in the output file? (0 indexed)')

args=parser.parse_args()

qcol = int(args.qcol)
rcol = int(args.rcol)
classes = args.classes

print 'Querying CDR3s from column',args.qcol,' of ',args.qfile,' in column ',args.rcol,' of ',args.rfile,'...'

dCLASS = {}

try:
	print len(classes),'classes:'
	for i in range(0,len(classes)):
		item = classes[i]
		dCLASS[i] = item
		print item
except TypeError:
	print classes,'classes'

qtcr_list = []			# the 			
qextra_list = []		# Will have an entry for every TCR; that entry will be redundant between TCRs associated with, say, the same clonotype

qextra_names = []		# just the header names associated with the extra data columns

qdelim = choose_delimiter(args.qfile)

k1 = 0
c = 0
t = 0
for line in fileinput.input([args.qfile]):
	llist = split(line,str(qdelim))						# Delimiter auto-detected here, can change to hard code
	if (args.qheader == True) and (k1 == 0):
		if len(args.extracolumn) > 0:
			for pt in args.extracolumn:
				try:
					qextra_names.append(llist[int(pt)])
				except IndexError:
					print pt,llist
	else:
		tcrinfo = llist[int(qcol)].rstrip(',')
		if len(args.extracolumn) > 0:
			extras = []
			for pt in args.extracolumn:
				extras.append(llist[int(pt)])
		tcrlist = split(tcrinfo,args.tcrdelim)

		for item in tcrlist:
			if args.classes == None:
				qtcr_list.append(item)
				t = t+1
			else:
				datum = split(item,args.classdelim)
				classname = datum[0]
				try:
					tcr = str(datum[1])
				except IndexError:
					print datum
				tcr = tcr[0:len(tcr)-int(args.qtrim)]
	#			print llist[0],tcr
				qtcr_list.append([classname,tcr])
				t = t+1
			qextra_list.append(extras)					# This looks weird but yes every time you append a new TCR, you append the extras list for that line
		c = c+1
	if k1 == 10:
		print "qtcr_list",qtcr_list
		print "qextra_list",qextra_list
	k1 = k1+1
fileinput.close()

print k1,' lines'
print c,' clonotypes'
print t,' TCRs'

#print qtcr_list

try:
	[q_classes,q_tcrs] = zip(*qtcr_list)
except ValueError:
	q_tcrs = qtcr_list
	q_classes = len(qtcr_list)*['ALL']

#print q_classes,q_tcrs
#print qextra_list
#q_extras = zip(*qextra_list)

for name in classes:
	num = countif(q_classes,name)
	print num,name

outfile = str(args.rfile).rstrip('.tsv''.txt''.csv') + '_' + str(args.qfile).rstrip('.tsv''.txt''.csv') + '.tsv'
output = open(outfile,'w')

rdelim = choose_delimiter(args.rfile)

k2 = 0
m = 0
for line in fileinput.input([args.rfile]):
	if (args.rheader == True) and (k2 == 0):
		newline = line.rstrip('\n''\r')
		output.write('%(newline)s' % vars())
		output.write('\tMatched_Num' % vars())
		output.write('\tMatched_Class' % vars())
		for extrahead in qextra_names:
			output.write('\t%(extrahead)s' % vars())
	else:
		newline = line.rstrip('\n''\r')
		output.write('%(newline)s' % vars())
		llist = split(line,str(rdelim))		
		data = llist[int(rcol)]
		r_tcr = data[0:len(data)-int(args.rtrim)]

		num = 0
		q_extras = []
		try:
			match_index = [ i for i, x in enumerate(q_tcrs) if x == r_tcr]
#			print match_index
			if k2 == 10:
				print "match_index",match_index
			num = len(match_index)
			matchtcrlist = []
			matchclasslist = []
			q_extradatalist = []
			for i in range(0,len(match_index)):
				matchtcrlist.append(q_tcrs[match_index[i]])
				matchclasslist.append(q_classes[match_index[i]])
				if len(args.extracolumn) > 0:
					q_extradatalist.append(qextra_list[match_index[i]])
#					for j in range(0,len(args.extracolumn)):
#						q_extradata.append(q_extras[j][match_index[i]])
			if len(matchclasslist) > 0:
				output.write('\t%(num)s\t%(matchclasslist)s' % vars())
			else:
				output.write('\t%(num)s' % vars())			
			q_extradata = zip(*q_extradatalist)
#			print len(q_extradata)
			for pt in q_extradata:
				data = []
				for datum in pt:
					if len(datum) > 0:
						data.append(datum)
				output.write('\t%(data)s' % vars()) 
			m =  m + 1


 		except ValueError:
			matchtcrlist = ''
			matchclasslist = ''
			num = ''
			output.write('\t%(num)s' % vars())
#		print 'q_extradata',q_extradata

	output.write('\n' % vars())
	k2 = k2 + 1
fileinput.close()

print k2,' lines in ',args.rfile
print m,args.rfile,'CDR3s found by',args.qfile,'query'

batmansays = ['BOOM!','KAPOW!','SHAZZAM!','WHAMMY!','VRONK!','SPLAT!','BANG!','WHAP!','ZOWIE!','SPLAT!','BAM!(just kidding)','CLANK!(just kidding)','AWK!(just kidding)']
print args.qfile,'in',args.rfile,'...'+str(random.choice(batmansays))
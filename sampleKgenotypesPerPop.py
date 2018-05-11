#! /usr/bin/env python

# Author: David Marques, davidalexander.marques [at] eawag.ch, (c) 09.11.2014
# Title: sampleKgenotypesPerPop.py
# Written in Python 2.6.6
# What it does: This script takes a VCF file and creates a new VCF file:
# (1) without missing data and
# (2) with a number k of genotypes per population
# by randomly picking a subsample of available genotypes per population
# which are then put into synthetic individuals (1001.pop1, 1002.pop1, 1001.pop2 etc.)
# If there are not enough genotypes in a population, the site will be removed from the dataset.
# Re-assembled and removed sites are reported in the verbose output (default, -v)
#
# 2014-11-11: Bugfix, simplified list comprehension in the data section
# 2017-08-10: Different k for each population
# 2018-05-11: New, providing a population file needed to assign individuals to populations
#             Added optional input from STDIN and output to STDOUT, fixed broken pipe bug
#             Added functionality for in- and / or output in gzipped format
#             Added functionality to either specify k for each population or for all populations at once

from sys import *
import argparse, re, random
from collections import defaultdict
from itertools import repeat
import gzip

# Fixes a broken pipe issue if STDOUT is piped into head
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)

parser = argparse.ArgumentParser(description='Subsample VCF file to k non-missing genotypes per population')

parser.add_argument('-i', '--input', dest='i', help="Input file in VCF format or enter '-' for STDIN [required]", required=True)
parser.add_argument('-o', '--output', dest='o', help="Output file name or enter '-' for STDOUT [required]", required=True)
parser.add_argument('-p', '--popfile', dest='p', help="Population file name [required], format: file with two TAB-separated columns, individual/tpopulation, each individual in the VCF file as one line", required=True)
parser.add_argument('-k', '--genotypes', dest='k', help="Number of genotypes to sample per population, either an integer (applies to all populations) or in the format: 'pop1:12,pop2:10,pop3:11' for different cut-offs in different populations [required]", required=True)
parser.add_argument('-s', '--seed', dest='s', type=int, help="Seed for random number generator [default=20]", default=20)
parser.add_argument('-v', '--verbose', action='store_true', help="If -v is specified, reporting on input, output and filtered sites will be displayed", default=False)

args = parser.parse_args()

# Print verbose information to screen (STDERR)
if args.verbose:
	print >> stderr, "Subsampling the VCF file with %s  genotypes per population" % (args.k)

# Read information from population file
inputP=open(args.p,'r')
indpop=defaultdict(str)
for Line in inputP:
	columns=Line.strip("\n").split("\t")
	indpop[columns[0]]=columns[1]
inputP.close()

# Transform args.k into dictionary with population names (index) and requested number of individuals (value)
if len(args.k.split(",")) == 1:
	tmp=[(x,args.k) for x in set(indpop.values())]
else:
	tmp=[x.split(":") for x in args.k.split(",")]
kperpop=defaultdict(str)
poplist=[x[0] for x in tmp]
for i in tmp:
	kperpop[i[0]]=int(i[1])

# Checks whether populations given in argument -k are identical to populations in popfile
poplist2=list(set(indpop.values()))
if sorted(poplist2)!=sorted(poplist):
	print >> stderr, "Error: population list entered not identical with population list in VCF file\nEntered: "+' '.join(poplist)+"\nVCF file: "+' '.join(poplist2)
	exit(0)

# Determine input and output options: files, gzipped-files or STDIN/STDOUT
if args.i == "-":
	inputF=stdin
else:
	# Decide whether to open gzipped file or not
	if args.i.endswith(".gz"):
		inputF=gzip.open(args.i,'r')
	else:
		inputF=open(args.i,'r')
if args.o != "-":
	# Decide whether to write gzipped file or not
	if args.o.endswith(".gz"):
		outputF=gzip.open(args.o, 'wb')
	else:
		outputF=open(args.o, 'w')

# seed for random number generator
random.seed(args.s)

# Count passing / failing sites for reporting in screen (STDERR)
sitefail=0
sitepass=0

for Line in inputF:
	# HEADER SECTION of VCF file
	if re.match('^#',Line): 
		if re.match('^##',Line): # writes header information into outfile
			# Decide whether to write into STDOUT or file
			if args.o =="-":
				stdout.write(Line)
			else:
				outputF.write(Line)
		if re.match('^##',Line) is None: # header with individuals / popinfo for parsing / changing
			header=Line.strip("\n").split("\t")
			indid=header[9:len(header)]  # header now contains all the individual IDs
			noperpop=defaultdict(str)
			indexpop=defaultdict(str)
			# Assign individuals in VCF file to populations in the popfile
			pop=[indpop[x] for x in indid]
			for i in poplist:
				noperpop[i]=pop.count(i)  # counts no. of individuals per population
				indexpop[i]=[j for j,n in enumerate(pop) if n==i]  # get index of individuals from population i
			if True in [noperpop[x]<kperpop[x] for x in noperpop]:
				print >> stderr, "Error: one or more populations have less than k genotypes"
				exit(0)
			# make new header: kperpop times each population name with a fake individual name
			populationlist=[item for sublist in [[x]*kperpop[x] for x in poplist] for item in sublist]
			individuallist=[item for sublist in [map(str,range(1001,1001+kperpop[x])) for x in poplist] for item in sublist]
			newheaderlist= [i+j+k for i,j,k in zip(individuallist,"."*len(individuallist),populationlist)]
			# Decide whether to write into STDOUT or file
			if args.o =="-":
				stdout.write('\t'.join(header[0:9])+'\t'+'\t'.join(newheaderlist)+"\n")
			else:
				outputF.write('\t'.join(header[0:9])+'\t'+'\t'.join(newheaderlist)+"\n")
	# DATA SECTION of VCF file
	else:
		columns=Line.strip("\n").split("\t")
		genotypecolumns=columns[9:len(columns)]
		tmp=[x.split(":") for x in genotypecolumns]
		genotypes=[x[0] for x in tmp]
		newgenotypecolumns=[]
		genotypesperpop=defaultdict(str)
		fail=0
		for i in poplist:
			popdata=[genotypes[j] for j in indexpop[i]]
			subsample=[]
			missing=["./.","."]
			popnonmissing=[j for j,n in enumerate(popdata) if n not in missing]
			popgenotypes=[genotypecolumns[j] for j in indexpop[i]]
			popgenotypesnonmissing=[popgenotypes[j] for j in popnonmissing]
			if len(popnonmissing) > kperpop[i]:
				newgenotypecolumns.extend(random.sample(popgenotypesnonmissing,kperpop[i]))
			elif len(popnonmissing) == kperpop[i]:
				newgenotypecolumns.extend(popgenotypesnonmissing)
			else:
				fail=1
				sitefail+=1
				break
		if fail==0:
			sitepass+=1
			# Decide whether to write into STDOUT or file
			if args.o =="-":
				stdout.write('\t'.join(columns[0:9])+'\t'+'\t'.join(newgenotypecolumns)+"\n")
			else:
				outputF.write('\t'.join(columns[0:9])+'\t'+'\t'.join(newgenotypecolumns)+"\n")

# Reporting to screen (STDERR)
if args.verbose:
	print >> stderr, "\tInput file: %d Total sites of which" %(sitefail+sitepass)
	print >> stderr, "\t%d sites were successfully re-assembled to a dataset without missing data and" %(sitepass)
	print >> stderr, "\t%d sites were excluded due to too few genotypes present." %(sitefail)

# Close input and output files if applicable
if args.i != "-":
	inputF.close()
if args.o != "-":
	outputF.close()
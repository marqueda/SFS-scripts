#! /usr/bin/env python

# Author: David Marques, davidalexander.marques [at] eawag.ch, (c) 09.11.2014
# Title: choosekgenotypes.py
# Written in Python 2.6.6
#
# What it does: This script takes a VCF file and creates a new VCF file:
# (1) without missing data and
# (2) with an equal number of genotypes per population
# by randomly drawing a subsample of available genotypes per population
# which are then put into synthetic individuals (1001.pop1, 1002.pop1, 1001.pop2 etc.)
# If there are not enough genotypes in a population, the site will be removed from the dataset.
# Re-assembled and removed sites are reported in the verbose output (default, -v)
#
# NOTE: Individuals need to be labeled as "individual.popcode" with both individual identifier and population code
# that only contain the characters A-z and 0-9, separated by a period "."
#
# Usage: choosekgenotypes.py -i <infile>.vcf -o <outfile>.vcf -k INT
# -k INT  This number gives the number of genotypes per population to randomly draw from the available genotypes.
# -s INT  Seed for random number generator [default=20]
#
# 2014-11-11: Bugfix, simplified list comprehension in the data section

from sys import *
import argparse, re, random
from collections import defaultdict
from itertools import repeat

parser = argparse.ArgumentParser(description='Subsample VCF file to k non-missing genotypes per population')

parser.add_argument('-i', '--input', dest='i', help="Input file in VCF format [required]", required=True)
parser.add_argument('-o', '--output', dest='o', help="Output file name [required]", required=True)
parser.add_argument('-k', '--genotypes', dest='k', type=int, help="Number of genotypes to choose per population [required]", required=True)
parser.add_argument('-s', '--seed', dest='s', type=int, help="Seed for random number generator [default=20]", default=20)
parser.add_argument('-v', '--verbose', action='store_true', help="If -v is specified, reporting on input, output and filtered sites will be displayed", default=False)

args = parser.parse_args()

if args.verbose:
    print "\nSubsampling the VCF file with %d genotypes per population" % (args.k)

inputF=open(args.i,'r')
outputF=open(args.o, 'w')

# seed for random number generator
random.seed(args.s)

sitefail=0
sitepass=0

for Line in inputF:
	# HEADER SECTION of VCF file
	if re.match('^#',Line): 
		if re.match('^##',Line): # writes header information into outfile
			outputF.write(Line)
		if re.match('^##',Line) is None: # header with individuals / popinfo for parsing / changing
			header=Line.strip("\n").split("\t")
			indid=header[9:len(header)]  # header now contains all the individual IDs
			pop=[]
			for i in range(0,len(indid)):
				sub=indid[i].split('.')[1]  # extracts population code from the string "individual.popcode[.etc]"
				pop.append(sub)
			poplist=list(set(pop))
			noperpop=defaultdict(str)
			indexpop=defaultdict(str)
			for i in poplist:
				noperpop[i]=pop.count(i)  # counts no. of individuals per population
				indexpop[i]=[j for j,n in enumerate(pop) if n==i]  # get index of individuals from population i
			if True in [noperpop[x]<=args.k for x in noperpop]:
				print "Error: one or more populations have less than k genotypes"
				#exit(0)
			# make new header: k times each population name with a fake individual name
			populationlist=[x for item in noperpop.keys() for x in repeat(item, args.k)]
			individuallist=map(str,range(1001,1000+args.k+1)*len(noperpop.keys()))
			newheaderlist=[i+j+k for i,j,k in zip(individuallist,"."*len(noperpop.keys())*args.k,populationlist)]
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
		missing=["./.","."]
			popnonmissing=[j for j,n in enumerate(popdata) if n not in missing]
			popgenotypes=[genotypecolumns[j] for j in indexpop[i]]
			popgenotypesnonmissing=[popgenotypes[j] for j in popnonmissing]
			if len(popnonmissing) > args.k:
				newgenotypecolumns.extend(random.sample(popgenotypesnonmissing,args.k))
			elif len(popnonmissing) == args.k:
				newgenotypecolumns.extend(popgenotypesnonmissing)
			else:
		fail=1
				sitefail+=1
				break
		if fail==0:
			sitepass+=1
			outputF.write('\t'.join(columns[0:9])+'\t'+'\t'.join(newgenotypecolumns)+"\n")

if args.verbose:
	print "\tInput file: %d Total sites of which" %(sitefail+sitepass)
	print "\t%d sites were successfully re-assembled to a dataset without missing data and" %(sitepass)
	print "\t%d sites were excluded due to too few genotypes present." %(sitefail)

inputF.close()
outputF.close()

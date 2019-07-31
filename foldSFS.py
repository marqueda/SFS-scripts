#! /usr/bin/env python

# Author: David Marques, davidalexander.marques [at] eawag.ch, (c) 24.05.2018
# Title: foldSFS.py
# Written in Python 2.7.6
# What it does: This script folds an SFS file in fastsimcoal2 format of arbitrary dimension
# CAUTION: Currently only takes SFS files with a single SFS replicate!

# EXPANSION IDEAS: Make suitable for SFS with multiple replicates

from sys import *
import argparse, re
import gzip
import numpy as np
from itertools import product

# Fixes a broken pipe issue if STDOUT is piped into head
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)

parser = argparse.ArgumentParser(description='Converts a VCF file to an SFS file of arbitrary dimension by counting alleles from genotypes')

parser.add_argument('-i', '--input', dest='i', help="Input unfolded SFS file [_DAF*,_jointDAF*,_DSFS*] or enter '-' for STDIN [required]", required=True)
parser.add_argument('-o', '--output', dest='o', help="Prefix for output folded SFS file or enter '-' for STDOUT [required]", required=True)

args = parser.parse_args()

# Folding function for derived allele frequency spectrum, given SFS and population counts (two numpy arrays)
def fold_sfs(daf,popcount):
	# Create SFS array output MAF
	maf=np.zeros(popcount+1,dtype='float')
	# Get list of all possible coordinates in the SFS (indices as in fastsimcoal2 MSFS format)
	idxs=list(product(*[range(0,x+1) for x in popcount]))
	# Loop over SFS indices that need to be folded (entries with AF < 50% of total alleles sampled)
	for j in [i for i, x in enumerate([sum(x)<(float(sum(popcount))/2) for x in idxs]) if x]:
		maf[idxs[j]]=daf[idxs[j]]+daf[tuple(popcount-idxs[j])]
	# Fill diagonals of the MAF with observed values, i.e. the entries in which the AF = 50% of total alleles sampled
	for j in [i for i, x in enumerate([sum(x)==(float(sum(popcount))/2) for x in idxs]) if x]:
		maf[idxs[j]]=float(daf[idxs[j]]+daf[tuple(popcount-idxs[j])])/2
	return maf

# Determine input options: normal, gzipped or STDIN
if args.i == "-":
	inputF=stdin
else:
	# Determine whether observed or expected SFS is input
	if re.search('obs',args.i):
		suffix="obs"
	else:
		suffix="txt"
	# Decide whether to open gzipped file or not
	if args.i.endswith(".gz"):
		inputF=gzip.open(args.i,'r')
	else:
		inputF=open(args.i,'r')

# Read infile into buffer
data=inputF.readlines()

# Determine whether data is in 1D, 2D or multi SFS format
if re.search('d',data[1]):
	if re.search('d',data[2]):
		# joint-SFS format
		# Read number of SFS replicates in the file
		nsfs=int(re.sub(r' obs.*\n','',data[0]))
		# Read population indices
		popcol=data[1].strip('\n').split('\t')[1].strip('0').strip('_')[1:]
		poprow=data[2].strip('\n').split('\t')[0].strip('0').strip('_')[1:]
		# Count number of alleles per population from column and row names
		popcoln=len(data[1].strip('\n').split('\t')[1:])-1
		# Get index of data starts for each SFS in the file
		sfsstarts=[i for i, item in enumerate(data) if re.search('d'+str(poprow)+'_0',item)]
		if len(sfsstarts)>1:
			# If there is >1 SFS replicate per file, the difference between start and end (minus row 0) is the no. of alleles for the population in the rows
			poprown=sfsstarts[1]-sfsstarts[0]-1
		else:
			# If there is 1 SFS replicate per file, the file length -3 is the no. of alleles for the population in the rows
			poprown=len(data)-3
		# index of matching SFS starts 
		popcount=np.array([poprown,popcoln],dtype='int')
		# Output MAF / folded 2D-SFS in fastsimcoal2 format
		# Decide whether to write into STDOUT or file
		if args.o =="-":
			# Output 2D-SFS to STDOUT: first header (no. of observations) and second header (TAB dX_0,dX_1,etc.) line, followed by data line(s) with row headers (dY_0,dY_1,etc.)
			stdout.write(str(nsfs)+" observations\n")
			stdout.write("\td"+popcol+"_"+("\td"+popcol+"_").join([str(x) for x in range(0,popcount[1]+1)])+"\n")
			# Loop over multiple SFS replicates in the file
			for n in range(0,nsfs):
				# Create SFS array for derived AF and fill it with entries
				daf=np.zeros(popcount+1,dtype='float')
				for i in range(0,popcount[0]+1):
					daf[i,]=[float(x) for x in data[i+2+(poprown+1)*n].strip('\n').split('\t')[1:]]
				# Fold DAF into MAF / folded 2D-SFS
				maf=fold_sfs(daf,popcount)
				for i in range(0,maf.shape[0]):
					# Convert numbers individually to necessary decimals (1.0 -> 1, 0.2 -> 0.2, 0.23545 -> 0.23545)
					outdatastrings=[format(str(round(x, len(str(x%1))-2) if x%1!=0 else int(x))) for x in maf[i,:]]
					stdout.write("d"+poprow+"_"+str(i)+"\t"+"\t".join(outdatastrings)+"\n")
		else:
			# Wrote 2D-SFS file: first header (no. of observations) and second header (TAB dX_0,dX_1,etc.) line, followed by data line(s) with row headers (dY_0,dY_1,etc.)
			outputF=open(args.o+"_jointMAFpop"+poprow+"_"+popcol+"."+suffix, 'w')
			outputF.write(str(nsfs)+" observations\n")
			outputF.write("\td"+popcol+"_"+("\td"+popcol+"_").join([str(x) for x in range(0,popcount[1]+1)])+"\n")
			# Loop over multiple SFS replicates in the file
			for n in range(0,nsfs):
				# Create SFS array for derived AF and fill it with entries
				daf=np.zeros(popcount+1,dtype='float')
				for i in range(0,popcount[0]+1):
					daf[i,]=[float(x) for x in data[i+2+(poprown+1)*n].strip('\n').split('\t')[1:]]
				# Fold DAF into MAF / folded 2D-SFS
				maf=fold_sfs(daf,popcount)
				for i in range(0,maf.shape[0]):
					# Convert numbers individually to necessary decimals (1.0 -> 1, 0.2 -> 0.2, 0.23545 -> 0.23545)
					outdatastrings=[format(str(round(x, len(str(x%1))-2) if x%1!=0 else int(x))) for x in maf[i,:]]
					outputF.write("d"+poprow+"_"+str(i)+"\t"+"\t".join(outdatastrings)+"\n")
	else:
		# 1D-SFS format
		# Read number of SFS replicates in the file
		nsfs=int(re.sub(r' obs.*\n','',data[0]))
		# Read number of alleles into numpy array
		popcount=np.array(len(data[1].strip('\n').split('\t'))-1,dtype='int').reshape(1)
		# Read population index
		popcol=data[1].strip('\n').split('\t')[0].strip('0').strip('_')[1:]
		# Output MAF / folded 1D-SFS in fastsimcoal2 format
		# Decide whether to write into STDOUT or into file
		if args.o =="-":
			# Output 1D-SFS to STDOUT: first header (no. of observations) and second header (dX_0,dX_1,etc.) line, followed by data line(s)
			stdout.write(str(nsfs)+" observations\n")
			stdout.write("d"+popcol+"_"+("\td"+popcol+"_").join([str(x) for x in range(0,popcount+1)])+"\n")
			# Loop over multiple SFS replicates in the file
			for n in range(0,nsfs):
				# Read entries from input DAF into numpy array
				daf=np.array([float(x) for x in data[n+2].strip('\n').split('\t')[0:]],dtype='float')
				# Fold DAF into MAF / folded 1D-SFS
				maf=fold_sfs(daf,popcount)
				# Convert numbers individually to necessary decimals (1.0 -> 1, 0.2 -> 0.2, 0.23545 -> 0.23545)
				outdatastrings=[format(str(round(x, len(str(x%1))-2) if x%1!=0 else int(x))) for x in maf.flatten()[range(0,popcount+1)]]
				stdout.write("\t".join(outdatastrings)+"\n")
		else:
			# Write 1D-SFS file: first header (no. of observations) and second header (dX_0,dX_1,etc.) line, followed by data line(s)
			outputF=open(args.o+"_MAFpop"+popcol+"."+suffix, 'w')
			outputF.write(str(nsfs)+" observations\n")
			outputF.write("d"+popcol+"_"+("\td"+popcol+"_").join([str(x) for x in range(0,popcount+1)])+"\n")
			# Loop over multiple SFS replicates in the file
			for n in range(0,nsfs):
				# Read entries from input DAF into numpy array
				daf=np.array([float(x) for x in data[n+2].strip('\n').split('\t')[0:]],dtype='float')
				# Fold DAF into MAF / folded 1D-SFS
				maf=fold_sfs(daf,popcount)
				# Convert numbers individually to necessary decimals (1.0 -> 1, 0.2 -> 0.2, 0.23545 -> 0.23545)
				outdatastrings=[format(str(round(x, len(str(x%1))-2) if x%1!=0 else int(x))) for x in maf.flatten()[range(0,popcount+1)]]
				outputF.write("\t".join(outdatastrings)+"\n")
else:
	# MSFS format
	# Read number of SFS replicates in the file
	nsfs=int(re.sub(r' obs.*\n','',data[0]))
	# Read number of populations into numpy array
	popcount=np.array([int(x) for x in data[1].strip('\n').split('\t')][1:],dtype='int')
	# Output MAF / folded multi-dimensional SFS (MSFS) in fastsimcoal2 format
	# Decide whether to write into STDOUT or into file
	if args.o =="-":
		# Output MSFS to STDOUT: first header (no. of observations) and second header (deme and sample sizes) line, followed by data line(s)
		stdout.write(str(nsfs)+" observations. No. of demes and samples are on next line\n")
		stdout.write(str(len(popcount))+"\t"+"\t".join([str(x) for x in popcount])+"\n")
		# Loop over multiple SFS replicates in the file
		for n in range(0,nsfs):
			# Read entries from input DAF into numpy array
			entries=[float(x) for x in data[n+2].strip('\n').split('\t')]
			daf=np.array(entries,dtype='float').reshape(popcount+1)
			# Fold DAF into MAF / folded MSFS
			maf=fold_sfs(daf,popcount)
			# Convert numbers individually to necessary decimals (1.0 -> 1, 0.2 -> 0.2, 0.23545 -> 0.23545)
			outdatastrings=[format(str(round(x, len(str(x%1))-2) if x%1!=0 else int(x))) for x in maf.flatten()]
			stdout.write("\t".join(outdatastrings)+"\n")
	else:
		# Write MSFS file: first (no. of observations) and second (deme and sample sizes) header line, followed by data line(s)
		outputF=open(args.o+"_MSFS."+suffix, 'w')
		outputF.write(str(nsfs)+" observations. No. of demes and samples are on next line\n")
		outputF.write(str(len(popcount))+"\t"+"\t".join([str(x) for x in popcount])+"\n")
		# Loop over multiple SFS replicates in the file
		for n in range(0,nsfs):
			# Read entries from input DAF into numpy array
			entries=[float(x) for x in data[n+2].strip('\n').split('\t')]
			daf=np.array(entries,dtype='float').reshape(popcount+1)
			# Fold DAF into MAF / folded MSFS
			maf=fold_sfs(daf,popcount)
			# Convert numbers individually to necessary decimals (1.0 -> 1, 0.2 -> 0.2, 0.23545 -> 0.23545)
			outdatastrings=[format(str(round(x, len(str(x%1))-2) if x%1!=0 else int(x))) for x in maf.flatten()]
			outputF.write("\t".join(outdatastrings)+"\n")

# Close input file if applicable
if args.i != "-":
	inputF.close()

# Close output file if applicable
if args.o != "-":
	outputF.close()
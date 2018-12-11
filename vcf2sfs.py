#! /usr/bin/env python

# Author: David Marques, davidalexander.marques [at] eawag.ch, (c) 24.05.2018
# Title: vcf2sfs.py
# Written in Python 2.7.6

# What it does: This script converts VCF files into site frequency spectra (SFS) files
# of the respective dimensions in fastsimcoal2 format (1D, 2D or multi-D)
# by simply counting alleles from genotypes for sites without missing data.

# IMPORTANT: Sites with missing data are discarded. Use my other script choosekgenotypes.py 
# first to subsample genotypes per population to get rid of missing data.

# 2018-09-07: Added functionality to output SFS for non-overlapping windows in the genome, with
#             (1) windows stored as replicates ("observations") in a single file and
#             (2) window coordinates stored in a separate BED file.
# 2018-09-07: Added functionality to output bootstrapped SFS either by
#             (1) bootstrapping all sites randomly, with replacement, or by
#             (2) bootstrapping blocks in the genome as specified by -w
#                 (also resulting in a BED file with block coordinates)
#             Both versions create a folder bootSFS or blockbootSFS with subfolders rep1,...,repN
#             containing the bootstrapped SFS
# 2018-12-11: Fixed several bugs in the -w option:
#             - if run with 'bp', site counting starts at site 0
#             - the output BED file now correctly displays start/end position corresponding
#               to fixed-size windows ('bp') or windows for sequenced sites ('sites')
#             Added new functionality: -f to read an input BED file with window definitions
#             WARNING: Windows in the BED file need to be sorted by chromosome / position as in the VCF file
#                      Windows in the BED file cannot be overlapping (=no sliding windows)!

# Load Python modules
from sys import *
import argparse, re, os
from collections import defaultdict
import gzip
import numpy as np
from itertools import product

# Fixes a broken pipe issue if STDOUT is piped into head
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)

# Read input arguments using the argparse module
parser = argparse.ArgumentParser(description='Converts a VCF file to an SFS file in fastsimcoal2 format by counting alleles from genotypes for sites without missing data')
parser.add_argument('-i', '--input', dest='i', help="Input file in VCF format or enter '-' for STDIN [required]", required=True)
parser.add_argument('-o', '--output', dest='o', help="Prefix for SFS filename or enter '-' for STDOUT [required]", required=True)
parser.add_argument('-p', '--popfile', dest='p', help="Population file name [required]. Format: file with two TAB-separated columns, individual\tpopulation, for each individual in the VCF file (in the same order), one per line", required=True)
parser.add_argument('-q', '--poporder', dest='q', help="Desired order of populations for the SFS outfile [optional, default: alphabetic order]. Format: pop3,pop1,pop2", required=False,default="XX")
parser.add_argument('-w', '--windowSFS', dest='w', help="Compute SFS in non-overlapping windows [optional, default: off]. Format: [chr/nochr],[sites/bp],size [chr/nochr]: indicate whether the windows should respect chromosome boundaries (=chr) or not (=nochr) [sites/bp]: give 'sites' or 'bp' to indicate whether the window size is specified in no. of sites in the VCF file or no. of base pairs along a chromosome. Also outputs a BED file with window coordinates.", required=False,default="XX")
parser.add_argument('-f', '--windowSFSBED', dest='f', help="Compute SFS in non-overlapping windows specificed in a BED file [optional, default: off]. Format: <BED file name>. WARNING: windows need to be non-overlapping and sorted by chromosome / position (same order as VCF file).", required=False,default="XX")
parser.add_argument('-b', '--bootstrapSFS', dest='b', help="Compute N bootstrap replicates of the observed SFS or block-bootstrap replicates with blocks specified with the Option -w [optional, default: off]. Format: N = number of replicates. Option incompatible with output to STDOUT (-o -)", required=False, default=0)
parser.add_argument('-v', '--verbose', action='store_true', help="If -v is specified, reporting on progress, number of processed and filtered sites will be displayed", required=False, default=False)
parser.add_argument('-r', '--reportafterNbp', dest='r', help="Reporting on progress after N basepairs. [optional, default: N = 1 Mbp]", required=False,default=1000000)
args = parser.parse_args()

# Check whether incompatible options are given and if, then abort.
if (args.b > 0) and (args.o == "-"):
	exit("ERROR: Bootstrap (-b) and output to STDOUT (-o -) are incompatible. Please specify an output file prefix with -o")

# Parse -w option for window-SFS if given
if args.w != "XX":
	winchr=args.w.split(',')[0]
	winmet=args.w.split(',')[1]
	winsiz=int(args.w.split(',')[2])
	coord=[]

# Parse -f option for window-SFS specified in BED file if given
if args.f != "XX":
	winbed=open(args.f,'r')
	wfi=0
	wfc=list()
	wfs=list()
	wfe=list()
	for Line in winbed:
		if not (re.match('^track',Line) or re.match('^browser',Line)):
			entry=Line.strip("\n").split("\t")
			wfc.append(entry[0])
			wfs.append(entry[1])
			wfe.append(entry[2])
	winbed.close()

# Parse -b option, number of bootstrap replicates, if given
if args.b > 0:
	args.b=int(args.b)

# Read population file and store information in a dictionary [indpop]
inputP=open(args.p,'r')
indpop=defaultdict(str)
for Line in inputP:
	columns=Line.strip("\n").split("\t")
	indpop[columns[0]]=columns[1]
inputP.close()

# Create and optionally sort list of populations [poplist], the latter if option -q is given
poplist=list(set(indpop.values()))
if args.q == "XX":
	poplist.sort(key=str.lower)
else:
	popsort=args.q.split(',')
	if len(set(poplist) & set(popsort)) == len(poplist):
		poplist=popsort
	else:
		exit("Aborted!\nPopulations given in the population order agrument (-q):\n\t"+args.q+"\ndo not match all the populations given in the population file (-p):\n\t"+",".join(poplist))

# Count number of individuals per population and store in array [popcount]
popcount=[indpop.values().count(k) for k in poplist]

# Create an array [sfs] containing two or more site-frequency spectra as one- or multidimensional numpy-arrays
# The first SFS (index 0) is a DUMMY SFS, but is important for calling the correct SFS (index > 0) later on
# Numpy arrays will contain allele counts, assuming diploid individuals in the VCF file
sfs=[]
sfs.append(np.zeros([x*2+1 for x in popcount],dtype='int'))
sfs.append(np.zeros([x*2+1 for x in popcount],dtype='int'))
sfsidx=1

# Function to add and switch to additional SFS in the array [sfs]
def addsfs():
	sfs.append(np.zeros([x*2+1 for x in popcount],dtype='int'))
	global sfsidx
	sfsidx+=1

# Parse -i option, normal input VCF file, gzipped input VCF file or STDIN
if args.i == "-":
	inputF=stdin
else:
	# Decide whether to open gzipped file or not
	if args.i.endswith(".gz"):
		inputF=gzip.open(args.i,'r')
	else:
		inputF=open(args.i,'r')

# Initialize variables before working across VCF file
# Count passing / failing sites and progress report intervals for reporting in screen (STDERR)
sitefail=0
sitepass=0
winstart=0
winsitecount=0
winbpcount=0
laststate='nowindow'
lastchr=""
reportcnt=1
reportintv=int(args.r)

# Working across VCF file
for Line in inputF:
	# HEADER SECTION of VCF file: parse header information on individuals / populations
	if re.match('^#',Line): 
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
			if True in [noperpop[x]!=indpop.values().count(x) for x in noperpop]:
				print >> stderr, "Error: individuals in population file and VCF file not identical."
				exit(0)
	# DATA SECTION of VCF file: fill SFS entries
	else:
		columns=Line.strip("\n").split("\t")
		# Initialize values for current and previous chromosme and position
		if lastchr=="":
			lastchr=columns[0]
			winend=columns[1]
		# Option Window-SFS: create new empty SFS after window boundary
		if args.w != "XX":
			if winmet=="sites":
				winstart=int(columns[1])-1
			winbpcount=int(columns[1])-winstart
			winsitecount+=1
			# Start new SFS and initialize variables at chromosome switch or window size exceeds limit
			if (winchr=="chr" and (lastchr != columns[0])) or (winmet=="sites" and (winsitecount > winsiz)) or (winmet=="bp" and (winbpcount > winsiz)):
				addsfs()
				winsitecount=1
				if winmet=="sites":
					coord.append(lastchr+"\t"+str(winstart)+"\t"+winend)
					winstart=int(columns[1])-1
				else:
					coord.append(lastchr+"\t"+str(winstart)+"\t"+str(winstart+winsiz))
					winstart=((int(columns[1])-1)/winsiz)*winsiz
		# Option Window-SFS-File: define start and end of windows, create new empty SFS after window boundary
		if args.f != "XX":
			# Check whether the new current window exists
			if (wfi<=(len(wfc)-1)):
				# Check whether the current site falls into the current window in the BED file
				if (columns[0]==wfc[wfi]) and (int(columns[1])>int(wfs[wfi])) and (int(columns[1])<=int(wfe[wfi])):
					# Empty the SFS if the last site did not fall into a window
					if laststate=='nowindow':
						# Empty last SFS window
						sfs[sfsidx]=np.zeros([x*2+1 for x in popcount],dtype='int')
					# Define current site as falling into the current window
					laststate='window'
				else:
					# If the previous site fell into the current window, switch to a new current window
					if laststate=='window':
						addsfs()
						wfi+=1
					else:
						# Empty last SFS window
						sfs[sfsidx]=np.zeros([x*2+1 for x in popcount],dtype='int')
					# Check whether the new current window exists
					if (wfi<=(len(wfc)-1)):
						# Check whether current site preceeds or exceeds new current window
						if (columns[0]==wfc[wfi]) and (int(columns[1])<=int(wfs[wfi])):
							# If the current site preceeds the new current window, empty last SFS window
							sfs[sfsidx]=np.zeros([x*2+1 for x in popcount],dtype='int')
							laststate='nowindow'
						elif (columns[0]==wfc[wfi]) and (int(columns[1])>int(wfs[wfi])) and (int(columns[1])<=int(wfe[wfi])):
							# If the current site matches the new current window, set laststate to active window
							laststate='window'
						else:
							# If the current site exceeds the new current window, switch to newer current window until
							# (a) current site preceeds or matches newer current window or (b) no more windows exist
							while (((columns[0]!=wfc[wfi]) or ((columns[0]==wfc[wfi]) and (int(columns[1])>int(wfe[wfi])))) and (wfi<=(len(wfc)-1))):
								addsfs()
								wfi+=1
								# Check whether the new current window exists
								if (wfi<=(len(wfc)-1)):
									# Check whether the current sites falls into the new current window in the BED file
									if (columns[0]==wfc[wfi]) and (int(columns[1])>int(wfs[wfi])) and (int(columns[1])<=int(wfe[wfi])):
										laststate='window'
									else:
										laststate='nowindow'
								else:
									laststate='nowindow'
					else:
						laststate='nowindow'
		# Sum distance (number of basepairs) and number of sites processed so far
		winend=columns[1]
		genotypecolumns=columns[9:len(columns)]
		tmp=[x.split(":") for x in genotypecolumns]
		genotypes=[x[0] for x in tmp]
		alleles=[x for s in [re.findall(r"[\d.]+",x) for x in genotypes] for x in s]
		# Check if the chromosome changed
		if lastchr != columns[0]:
			lastchr=columns[0]
			reportcnt=1
		# Check if there is missing data, if the site is an indel or a multiallelic SNP
		if ('.' in list(set(alleles))) or (len(columns[3])>1) or (len(columns[4].strip('<NON_REF>').strip(','))>1):
			sitefail+=1
		else:
			alleles=map(int,alleles)
			entry=[]
			for i in poplist:
				entry.append(sum([alleles[y] for y in [x for n in [[x*2,x*2+1] for x in indexpop.get(i)] for x in n]]))
			# Add entry to SFS
			sfs[sfsidx][tuple(entry)]+=1
			sitepass+=1
		if args.verbose and int(columns[1]) > reportcnt*reportintv:
			print >> stderr, "Progress: chromosome %s, position > %d. %d sites processed." %(columns[0],reportcnt*reportintv,sitefail+sitepass)
			reportcnt+=1

# Close input file if applicable
if args.i != "-":
	inputF.close()

# Determine whether to write 1D, 2D or multidimensional SFS format
if len(poplist) == 1:
	suffix="_DAFpop0"
	def printsfs(socket,sfs,popcount,sfsidx):
		# Write first Header line
		socket.write(str(sfsidx)+" observations\n")
		# Write second Header line
		socket.write("d0_"+"\td0_".join([str(x) for x in range(0,popcount[0]*2+1)])+"\n")
		# Write 1D-SFS
		for j in range(1,sfsidx+1):
			socket.write("\t".join([str(x) for x in sfs[j].flatten()])+"\n")
elif len(poplist) == 2:
	suffix="_jointDAFpop1_0"
	def printsfs(socket,sfs,popcount,sfsidx):
		# Write first Header line
		socket.write(str(sfsidx)+" observations\n")
		# Write second Header line
		socket.write("\td0_"+"\td0_".join([str(x) for x in range(0,popcount[0]*2+1)])+"\n")
		# Write 2D-SFS
		for j in range(1,sfsidx+1):
			for i in range(0,sfs[j].shape[1]):
				socket.write("d1_"+str(i)+"\t"+"\t".join([str(x) for x in sfs[j][:,i]])+"\n")
else:
	suffix="_DSFS"
	def printsfs(socket,sfs,popcount,sfsidx):
		# Write first Header line
		socket.write(str(sfsidx)+" observations. No. of demes and samples are on next line\n")
		# Write second Header line: demes, sample sizes
		socket.write(str(len(poplist))+"\t"+"\t".join([str(x*2) for x in popcount])+"\n")
		# Write multi-SFS
		for j in range(1,sfsidx+1):
			socket.write("\t".join([str(x) for x in sfs[j].flatten()])+"\n")

# Write window coordinate file (BED format) to disk, if applicable
if args.w != "XX":
	# Add entry of last window
	if winmet=="sites":
		coord.append(lastchr+"\t"+str(winstart)+"\t"+winend)
	else:
		coord.append(lastchr+"\t"+str(winstart)+"\t"+str(winstart+winsiz))
	coordF=open(args.o+suffix+".bed", 'w')
	coordF.write("\n".join(coord)+"\n")
	coordF.close()

# Write additional empty SFS for -f option if no sites were found
if args.f != "XX":
	if laststate=='window':
		addsfs()
		wfi+=1
	else:
		# Empty last SFS window
		sfs[sfsidx]=np.zeros([x*2+1 for x in popcount],dtype='int')
	while (wfi<=(len(wfc)-1)):
		addsfs()
		wfi+=1
	# Delete late SFS entry
	del sfs[sfsidx]
	sfsidx=sfsidx-1

# Decide whether to write SFS into STDOUT or file
if args.o =="-":
	printsfs(stdout,sfs,popcount,sfsidx)
else:
	outputF=open(args.o+suffix+".obs", 'w')
	printsfs(outputF,sfs,popcount,sfsidx)
	outputF.close()

# Optional: Bootstrap SFS or block-bootstrap SFS (if -w is specified)
if args.b > 0:
	if (args.w == "XX") and (args.f == "XX"):
		def boot_sfs(sfs,popcount,sfsidx): # SFS bootstrapping function
			# Create SFS array for bootstrapped
			bootSFS=[]
			bootSFS.append(np.zeros([x*2+1 for x in popcount],dtype='int'))
			bootSFS.append(np.zeros([x*2+1 for x in popcount],dtype='int'))
			# Get list of all possible coordinates in the SFS (indices as in fastsimcoal2 MSFS format)
			idxs=list(product(*[range(0,x+1) for x in [x*2 for x in popcount]]))
			# Resample sites from SFS
			boot=np.unique(np.random.choice(range(0,sfs[sfsidx].size),np.sum(sfs[sfsidx]),replace=True,p=sfs[sfsidx].ravel().astype(dtype='float')/np.sum(sfs[sfsidx])),return_counts=True)
			for j in range(0,len(boot[0])):
				bootSFS[1][idxs[boot[0][j]]]+=boot[1][j]
			return bootSFS
		if not os.path.exists("bootSFS"):
			os.makedirs("bootSFS")
		for i in range(0,args.b):
			os.makedirs("bootSFS/rep"+str(i+1))
			outputF=open("bootSFS/rep"+str(i+1)+"/"+args.o+suffix+".obs", 'w')
			printsfs(outputF,boot_sfs(sfs,popcount,sfsidx),popcount,1)
			outputF.close()
	else:
		def boot_wsfs(sfs,popcount,sfsidx): # SFS block-bootstrapping function
			# Create SFS array
			bootSFS=[]
			bootSFS.append(np.zeros([x*2+1 for x in popcount],dtype='int'))
			bootSFS.append(np.zeros([x*2+1 for x in popcount],dtype='int'))
			# Resample block-SFSs 
			for i in np.random.choice(range(1,sfsidx+1),sfsidx,replace=True):
				bootSFS[1]+=sfs[i]
			return bootSFS
		if not os.path.exists("blockbootSFS"):
			os.makedirs("blockbootSFS")
		for i in range(0,args.b):
			if not os.path.exists("blockbootSFS/rep"+str(i+1)):
				os.makedirs("blockbootSFS/rep"+str(i+1))
			outputF=open("blockbootSFS/rep"+str(i+1)+"/"+args.o+suffix+".obs", 'w')
			printsfs(outputF,boot_wsfs(sfs,popcount,sfsidx),popcount,1)
			outputF.close()

# Reporting to screen (STDERR)
if args.verbose:
	print >> stderr, "\tInput file: %d Total sites of which" %(sitefail+sitepass)
	print >> stderr, "\t%d sites without missing data were integrated into the SFS and" %(sitepass)
	print >> stderr, "\t%d sites were excluded because of missing genotypes, of being indels or multiallelic SNPs." %(sitefail)


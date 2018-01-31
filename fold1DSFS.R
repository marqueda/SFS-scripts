#! /usr/bin/env Rscript

# Author: David Marques, davidalexander.marques [at] eawag.ch
# Title: fold2DSFS.R

# Usage: fold1DSFS.R filename_MAFpopX.[obs/txt]
# Folds 1-dimensional SFS in fastsimcoal2 format

# Parse file name
args<-commandArgs(trailingOnly=T)
base=sub("_MAF.*","",args[1])
pops=strsplit(sub(".*pop","",args[1]),"[.]")[[1]][1]
suffix=substr(args[1], nchar(args[1])-3+1, nchar(args[1]))

# Read SFS
a=read.csv(args[1],skip=1,sep=" ")

# Function to fold 2D-SFS
derived1maf=function(derived_sfs1d){
 n1=ncol(derived_sfs1d)
 maf_1dsfs=matrix(0,ncol=n1)
 colnames(maf_1dsfs)=colnames(derived_sfs1d)
 maf_1dsfs[1:ceiling(n1/2)]=unlist(derived_sfs1d[1:ceiling(n1/2)])
 for(i in 0:(ceiling(n1/2)-1)){
   maf_1dsfs[1,i+1]=maf_1dsfs[1,i+1]+derived_sfs1d[1,n1-i]
 }
 maf_1dsfs
}

# Apply to read 2D-SFS
b=derived1maf(a)

# Round if observed SFS
if(grepl("obs",suffix)){b=round(b)}

# Print new, folded 1D-SFS
cat(paste("1 observations","\n",
paste(c("",colnames(b)),collapse="\t"),"\n",sep=""),file=paste(base,"folded_MAFpop",pops,".",suffix,sep=""))
suppressWarnings(write.table(b,paste(base,"folded_MAFpop",pops,".",suffix,sep=""),col.names=F,row.names=F,
quote=F,sep="\t",append=T))

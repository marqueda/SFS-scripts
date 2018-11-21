#! /usr/bin/env Rscript

# (c) David A. Marques, 2016
# Usage: SFStools.R -t [mto2D,mto1D,print1D,print2D] -i infileprefix [-z -s]
# 1) converts multidimensional SFS or 2D-SFS to marginal 2D- or 1D-SFS
# 2) visualizes (marginal) 2D and 1D SFS and compares observed and expected SFS (if both are available in a folder)

# WARNING: currently handles only SFS-files with a single SFS, not multiple SFSs

# Load libaries
library(optparse)

# Read input arguments
option_list = list(
  make_option(c("-t", "--tool"), type="character", default=NULL, 
              help="choose one of the following tools:\n\tmto2D : converts multidimensional SFS to marginal 2D-SFS\n\tmto1D : converts multidimensional SFS to marginal 1D-SFS\n\t2Dto1D : converts 2D-SFS to marginal 1D-SFS\n\tprint1D : visualizes 1D-SFS from files <prefix>_[M/D]AFpopX.[obs/txt]\n\tprint2D : prints 2D-SFS from files <prefix>_joint[M/D]AFpopX_Y.[obs/txt]", metavar="character"),
  make_option(c("-i", "--infile"), type="character", default=NULL, 
              help="infile prefix", metavar="character"),
  make_option(c("-z", "--nozeros"), action="store_true", default=F, 
              help="optional: removes monomorphic sites during conversion or printing", metavar="character"),
  make_option(c("-s", "--nosingletons"), action="store_true", default=F, 
              help="optional: removes singletons in each population during conversion or printing", metavar="character"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Check whether input parameters are provided for -t and -i
# if not provided, abort
if(is.null(opt$tool) | is.null(opt$infile)){
  stop(paste("Aborted. Pleae provide arguments.\nUsage: SFStools.R -t [mto2D,mto1D,print1D,print2D] -i infileprefix [-z -s]",sep=""))
}

# Option 1: conversion of multidimensional SFS to marginal SFS
if(opt$tool %in% c("mto2D","mto1D")){
  
  # Step 1.1: get list of all multidimensional SFS present in folder
  infiles=grep(paste(opt$infile,"_[MD]SFS.[ot][bx][st]",sep=""),list.files(),value=T)
  
  # Check 1.1: confirm presence of multidimensional SFS files, otherwise abort
  if(length(infiles)>0){
    
    # Step 1.2: Loop over each multidimensional SFS file
    for(j in infiles){
      
      # Step 1.3: Parsing multi-dimensional SFS
      {
        # Reading infile into vector
        infile<-scan(file=j,character(0), sep ="\t")
        # Determine: minor or derived SFS
        if(grepl("MSFS",j)){sfstype<-"M"}else{sfstype<-"D"}
        # Determine: observed or expected SFS
        if(grepl("obs",j)){ending<-".obs"}else{ending<-".txt"}
        # Get the number of SFS stored in the file
        nsfs<-as.integer(gsub(" obs.*","",infile[1]))
        # Get the number of populations from the file
        npop<-as.integer(infile[2])
        # Read population sizes from the file
        popsiz<-as.integer(infile[2+1:npop])
        
        # Update number of populations if ghost population present
        # -> updates npop and popsiz to the populations present
        # -> deletes the entry "0" for the ghost population in the infile for correct data parsing
        if(FALSE %in% popsiz>0){
          npop<-sum(popsiz>0)
          infile=infile[-c(which(popsiz==0)+2)]
          popsiz<-popsiz[popsiz>0]
        }
        
        # Issue warning of >1 SFS stored in file
        if(nsfs>1){
          stop("SFStools.R is not yet suited to work with several
               multi-dimensional SFS stored in a single file. Sorry!")
        }
        
        # Parse data into multidmensional array
        # Note: index = allele count + 1, i.e. arr[1,1,1] = entry (0,0,0)
        # IMPORTANT: fastsimcoal2 and dadi store multi-dimensional SFS in the format:
        #            (0,0,0),(0,0,1),(0,0,2),(0,1,0),(0,1,1),(0,1,2),...
        #            first iterating the last population!
        #            thus,
        #            a) population sizes are reversed while parsing the data into an array
        #            and then,
        #            b) the array is reversed to restore the population order as in the file
        arr<-array(as.numeric(infile[-c(1,2,2+0:npop)]),rev(popsiz+1))
        arr<-aperm(arr)

        # Remove monomorphic sites and singletons if needed
        # WARNING: monomorphic sites / singletons will be removed in the multidimensional SFS
        #          therefore, there will still be intries in the 
        #          (0,0), (0,1) and (1,0) categories in the marginal 2D-SFS
        #          and entries in the MAC=0 and MAC=1 categories in the marginal 1D SFS
        if(opt$nozeros){
          arr[1]=0
        }
        if(opt$nosingletons){
          arr[matrix(head(rep(c(2,rep(1,npop)),npop),npop*npop),nrow=npop,byrow=T)]=0
        }
      }
      
      # Option 1-A: convert multidimensional to marginal 1D-SFS
      if(opt$tool=="mto1D"){
        
        # Create matrix, one row per population and columns as in the population with most individuals
        marSFS<-matrix(NA,nrow=npop,ncol=max(popsiz+1))
        # Fill the matrix with the marginal SFS entries for each population
        # i.e. the sum of all entries across the mSFS for this 1D-category
        # e.g. all entries / SNPs where population 1 (size=10) has 9 major and 1 minor allele (-> category MAC=1)
        for(i in 1:npop){
          marSFS[i,1:dim(arr)[i]]<-apply(arr,i,function(x) sum(x))
        }
        colnames(marSFS)<-paste("d_",0:max(popsiz),sep="")
        
        # Write marginal 1D-SFS to separate files
        for(i in 0:(npop-1)){
          # Subset the matrix for each population and its number of entries
          outSFS=marSFS[i+1,]
          outSFS=outSFS[!is.na(outSFS)]
          write.table("1 observations",paste(opt$infile,"_",sfstype,"AFpop",i,ending,sep=""),col.names=F,row.names=F,quote=F)
          write.table(c(paste(names(outSFS),collapse=" "),paste(outSFS,collapse=" ")),
                      paste(opt$infile,"_",sfstype,"AFpop",i,ending,sep=""),col.names=F,row.names=F,append=T,quote=F)
        }
      }else
      # Option 1-B: convert multidimensional to marginal 2D-SFS
      {
        
        # Prepare list, store all pairwise marginal 2D-SFS in this list
        D2SFS<-list()
        # Get combination of all pairwise population comparisons
        pairind<-combn(npop,2)
        # Loop over pairwise combinations and compute marginal 2D-SFS,
        # i.e. the sum of all entries across the mSFS for this pairwise 2D-category
        # e.g. all entries / SNPs where population 1 (size=10) has 9 major and 1 minor allele
        #      and population 2 (size = 5) has 4 major and 1 minor allele (-> category (1,1))
        # IMPORTANT: result is a matrix with the population of higher number in rows and lower number in columns!
        #            --> this is the input format for fastsimcoal2: e.g. file_joint[MD]AFpop2_1.[txt/obs]
        for(i in 1:ncol(pairind)){
          D2SFS[[i]]<-apply(arr,rev(pairind[,i]),function(x) sum(x))
        }
  
        # Write all marginal pairwise 2D-SFSs to separate files
        for(i in 1:ncol(pairind)){
          write.table(paste("1 observations",sep=""),
                      paste(opt$infile,"_joint",sfstype,"AFpop",(pairind[,i][2]-1),"_",(pairind[,i][1]-1),ending,sep=""),
                      col.names=F,row.names=F,quote=F)
          write.table(paste(c("",paste(paste("d",rep((pairind[,i][1]-1),(popsiz[pairind[,i][1]]+1)),sep=""),
                                       0:popsiz[pairind[,i][1]],sep="_")),collapse="\t"),
                      paste(opt$infile,"_joint",sfstype,"AFpop",(pairind[,i][2]-1),"_",(pairind[,i][1]-1),ending,sep=""),
                      col.names=F,row.names=F,append=T,quote=F,sep="\t")
          write.table(cbind(paste(paste("d",rep((pairind[,i][2]-1),(popsiz[pairind[,i][2]]+1)),sep=""),
                                  0:popsiz[pairind[,i][2]],sep="_"),as.data.frame(D2SFS[i])),
                      paste(opt$infile,"_joint",sfstype,"AFpop",(pairind[,i][2]-1),"_",(pairind[,i][1]-1),ending,sep=""),
                      col.names=F,row.names=F,append=T,quote=F,sep="\t")
      }
    }
      
    }
  }else{
    stop(paste("Aborted. No multi-dimenstional SFS file (",opt$infile,"_[M/D]SFS.[obs/txt]) found in this directory!",sep=""))
  }
}else if (opt$tool %in% c("2Dto1D"))

# Option 2: conversion of 2D-SFS to 1D-SFS
{
  # Step 2.1: get list of all 2D-SFS present in folder
  infiles=grep(paste(opt$infile,"_joint[MD]AFpop",sep=""),list.files(),value=T)
  
  # Check 2.1: confirm presence of 2D-SFS files, otherwise abort
  if(length(infiles)>0){
    
    # Step 2.2: Loop over each 2D-SFS file
    for(j in infiles){
  
      # Step 2.3: Parsing 2D-SFS
      {
        # Determine: minor or derived SFS
        sfstype=NULL;if(grepl("MAF",j)){sfstype=c(sfstype,"M")}
        if(grepl("DAF",j)){sfstype=c(sfstype,"D")}
        # Determine: observed or expected SFS
        ending=NULL;if(grepl("obs",j)){ending=".obs"};if(grepl("txt",j)){ending=".txt"}
        # Get population names from filename
        popnames=rev(unlist(strsplit(gsub(ending,"",gsub(paste0(".*joint",sfstype,"AFpop"),"",j)),"_")))
        
        # Reading infile into vector
        if(sum(grepl("observation",scan(j,what = "character")))>0){
          infile<-as.matrix(read.table(j,skip=1,header=T))
        }else{
          infile<-as.matrix(read.table(j,skip=0,header=T))
        }
        
        # Reading the two population sizes
        npop<-2
        popsiz<-dim(infile)-1
        
        # Create multidmensional array for useability
        # Note: index = allele count + 1, i.e. arr[1,1] = entry (0,0)
        # IMPORTANT: fastsimcoal2 stores 2D-SFS in the following format:
        #            the higher order population in rows, lower oder in columns
        #            (0,0),(1,0),(2,0),(0,1),(1,1),(2,1)...
        #            thus, when read as a vector in R, the array needs to be transposed
        arr<-array(as.numeric(infile),popsiz+1)
        arr<-aperm(arr)
        
        # Remove monomorphic sites and singletons if needed
        # WARNING: monomorphic sites / singletons will be removed in the 2D-SFS
        #          therefore, there will still be intries in the 
        #          MAC=0 and MAC=1-categories in the marginal 1D SFS
        if(opt$nozeros){
          arr[1]=0
        }
        if(opt$nosingletons){
          arr[matrix(head(rep(c(2,rep(1,npop)),npop),npop*npop),nrow=npop,byrow=T)]=0
        }
      }

      # Step 2.4: convert 2D-SFS to marginal 1D-SFSs and print marginals to file
      for(i in 1:npop){
        marSFS<-matrix(apply(arr,i,function(x) sum(x)),nrow=1)
        colnames(marSFS)<-paste("d_",0:rev(popsiz)[i],sep="")
        # Recover population names for filename
        write.table("1 observation",paste(opt$infile,"_",sfstype,"AFpop",popnames[i],ending,sep=""),
                    col.names=F,row.names=F,quote=F)
        write.table(c(paste(colnames(marSFS),collapse=" "),paste(marSFS,collapse=" ")),
                    paste(opt$infile,"_",sfstype,"AFpop",popnames[i],ending,sep=""),
                    col.names=F,row.names=F,append=T,quote=F)
      }
    }
  }else{
    stop(paste("Aborted. No 2D-SFS files (",opt$infile,"_joint[M/D]AF_XpopY.[obs/txt]) found in this directory!",sep=""))
  }
}else if(opt$tool %in% c("print1D"))

# Option 3: Printing of 1D-SFS
{
  # Step 3.1: get list of all 1D-SFS present in folder
  infiles=grep(paste(opt$infile,"_[MD]AFpop",sep=""),list.files(),value=T)
  
  # Check 3.1: confirm presence of 1D-SFS files, otherwise abort
  if(length(infiles)>0){
    
    # Step 3.2: get population indices
    pop=unique(as.integer(unlist(lapply(strsplit(infiles,"pop"),function(x){
      lapply(strsplit(x[2],"[.]"),function(x){x[1]})}))))
    
    # Step 3.3: determine the presence/absence of minor and derived SFS
    sfstype=NULL;if(sum(grepl("MAF",infiles))>0){sfstype=c(sfstype,"M")};
    if(sum(grepl("DAF",infiles))>0){sfstype=c(sfstype,"D")}
    
    # Step 3.4: Loop over SFS-type (minor / derived)
    for(fold in sfstype){
      
      # Step 3.5: Print visualization to PDF for each SFS-type (minor / derived)
      pdf(paste(opt$infile,"_1Dmarginal",fold,"AFs.pdf",sep=""),height=5,width=4,useDingbats=F,pointsize=10)
      par(mfrow=c(1,1),mar=c(0,4,0,1),oma=c(2.5,0,1,0))
      obslist=list()
      type=NULL;pops=NULL;n=1
      
      # Step 3.6: Loop over single populations
      for(j in pop){
        
        # Step 3.7: Parse input files and visualize either observed, expected or both 1D-SFS
        # First: Load observed SFS for population j, if available
        if(sum(grepl("obs",grep(paste(fold,"AFpop",j,sep=""),infiles,value = T)))>0){
          obs<-as.matrix(read.table(paste(opt$infile,"_",fold,"AFpop",j,".obs",sep=""),skip=1,header=T))
          
          # Remove monomorphic sites and singletons for plotting if chosen as option
          # WARNING: these are informative in marginal 1D-SFS for the model fit
          if(opt$nozeros){
            obs[1]=0
          }
          if(opt$nosingletons){
            obs[2]=0
          }
          
          # Store observed SFS in list of marginal 1D-SFS
          obslist[[n]]=obs;n=n+1
          type=c(type,"obs")
          pops=c(pops,j)
        }
        
        # Second: Load expected SFS for population j, if available
        if(sum(grepl("txt",grep(paste(fold,"AFpop",j,sep=""),infiles,value = T)))>0){
          # Check if expected SFS has header line
          if(sum(grepl("observation",scan(paste(opt$infile,"_",fold,"AFpop",j,".txt",sep=""),what = "character")))>0){
            sim<-as.matrix(read.table(paste(opt$infile,"_",fold,"AFpop",j,".txt",sep=""),skip=1,header=T))
          }else{
            sim<-as.matrix(read.table(paste(opt$infile,"_",fold,"AFpop",j,".txt",sep=""),skip=0,header=T))
          }
          
          # Remove monomorphic sites and singletons for plotting if chosen as option
          # WARNING: these are informative in marginal 1D-SFS for the model fit
          if(opt$nozeros){
            sim[1]=0
            sim=sim/sum(sim)
          }
          if(opt$nosingletons){
            sim[2]=0
            sim=sim/sum(sim)
          }
          
          # Store expected SFS in list of marginal 1D-SFS
          obslist[[n]]=sim;n=n+1
          type=c(type,"sim")
          pops=c(pops,j)
        }
        
        # Option 3-A: Print observed, expected and observed - expected marginal SFS
        if((sum(grepl("txt",grep(paste(fold,"AFpop",j,sep=""),infiles,value = T)))>0)&
           (sum(grepl("obs",grep(paste(fold,"AFpop",j,sep=""),infiles,value = T)))>0)){

          # Multiply probabilities of expected SFS with number of observed SNPs to plot them on the same scale
          sim<-sum(obs)*sim
          obslist[[n-1]]=sim
          
          # Plot difference betewen observed and expected SFS
          par(fig=c(0,1,0,0.18))
          dif=sim/obs
          dif[is.na(dif)]=1
          dif[is.infinite(dif)]=1
          max(abs(dif))
          b=barplot((c(dif)-1)*100,names.arg = NA,col="grey",border=0,space = 0,
                    ylim=c(min(dif)-1,max(dif)-1)*100,ylab="",yaxt="n");box()
          axis(1,0:(length(dif)-1),at=b)
          axis(2)
          abline(h=0);mtext("% SNPs",2,line=3)
          
          par(fig=c(0,1,0.2,0.38),new=T)
          dif=sim-obs
          max(abs(dif))
          b=barplot(c(dif),names.arg = NA,col="grey",border=0,space = 0,
                  ylim=c(min(dif),max(dif)),ylab="");box()
          abline(h=0);mtext("No. of SNPs",2,at=grconvertY(0.6,"ndc","user"),line=3)
          xlimit=par("usr")[1:2]

          # Plot both observed and expected SFS
          par(fig=c(0,1,0.4,1),new=T)
          plot(range(b),c(0,max(c(obs,sim))),type="n",xlim=xlimit,xaxs="i",ylab="",xlab="",axes=F)
          lines(c(0,rep(b,each=2)+c(-0.5,0.5),max(b)+0.5),c(0,rep(obs,each=2),0),lwd=1,col="grey")
          lines(c(0,rep(b,each=2)+c(-0.5,0.5),max(b)+0.5),c(0,rep(sim,each=2),0),lwd=1,lty=2,col="blue")
          axis(2);box()
          legend("top",c("observed","expected","obs. - exp."),
                 title=paste("population",j),lty=c(1,2,NA),lwd=c(1,1,NA),col=c("grey","blue",1),
                 fill=c(NA,NA,"grey"),border=c(NA,NA,"grey"),bty="n")
          
        }else if((sum(grepl("obs",grep(paste(fold,"AFpop",j,sep=""),infiles,value = T)))>0))
          
        # Option 3-B: Print observed marginal SFS only
        {
          plot(0,0,ylim=c(0,max(c(obs))),xlim=c(0,length(obs)-1),type="n",axes=F,ylab="No. of SNPs",xlab="")
          lines(sort(c(0:(length(obs)-1)-0.5,0:(length(obs)-1)+0.5)),
                rep(obs,each=2),lwd=1,col="grey")
          axis(1,0:(length(obs)-1),lwd=0,lwd.ticks=1)
          axis(2);box()
          legend("top",c("observed"),title=paste("population",j),lty=c(1),lwd=1,col=c("grey"),bty="n")
        }else
          
        # Option 3-C: Print expected marginal SFS only
        {
          plot(0,0,ylim=c(0,max(c(sim))),xlim=c(0,length(sim)-1),type="n",axes=F,ylab="Proportion of sites",xlab="")
          lines(sort(c(0:(length(sim)-1)-0.5,0:(length(sim)-1)+0.5)),
                rep(sim,each=2),lwd=1,lty=2,col="blue")
          axis(1,0:(length(sim)-1),lwd=0,lwd.ticks=1)
          axis(2);box()
          legend("top",c("expected"),title=paste("population",j),lty=c(2),lwd=1,col=c("blue"),bty="n")
        }
      }
      
      # Step 3.7: Plot all observed / simulated on top of each other in different colors
      par(fig=c(0,1,0,1))
      
      plot(0,0,ylim=c(0,max(unlist(lapply(obslist,max)))),
           xlim=c(0,max(unlist(lapply(obslist,length)))-1),type="n",axes=F,
           ylab=c("No. of SNPs","Proportion of Sites")[1+1*(max(unlist(lapply(obslist,max)))<1)],xlab="")
      lcol<-colorRampPalette(c("#a5002688","#d7302788","#f46d4388","#fdae6188",
                               "#fee09088","#e0f3f888","#abd9e988","#74add188","#4575b488","#31369588"))
      for(j in 1:length(obslist)){
        lines(sort(c(0:(length(obslist[[j]])-1)-0.5,0:(length(obslist[[j]])-1)+0.5)),
              rep(obslist[[j]],each=2),lwd=1,lty=c(1,2)[1+1*(type[j]=="sim")],
              col=lcol(length(unique(pops)))[match(pops[j],unique(pops))])
      }
      legend("top",c(paste("population",0:(length(unique(pops))-1)),"observed","expected"),
             lty=c(rep(1,length(obslist)/2),1,2),col=c(lcol(length(unique(pops))),1,1),bty="n",lwd=1)
      axis(1,0:(max(unlist(lapply(obslist,length)))-1),lwd=0,lwd.ticks=1)
      axis(2);box()
      dev.off()
    }
  }else{
    stop(paste("No 1D-SFS files (",opt$infile,"_[M/D]AFpopX.[obs/txt]) found in this directory!",sep=""))
  }
}else if(opt$tool %in% c("print2D"))

# Option 4: Printing of 2D-SFS
{
  # Step 4.1: get list of all 2D-SFS present in folder
  infiles=grep(paste(opt$infile,"_joint[MD]AFpop",sep=""),list.files(),value=T)
  
  # Check 4.1: confirm presence of 2D-SFS files, otherwise abort
  if(length(infiles)>0){
    
    # Step 4.2: get indices of population comparisons
    pop=unique((unlist(lapply(strsplit(infiles,"pop"),function(x){lapply(strsplit(x[2],"[.]"),function(x){x[1]})}))))
    
    # Step 4.3: determine the presence/absence of minor and derived SFS
    sfstype=NULL;if(sum(grepl("MAF",infiles))>0){sfstype=c(sfstype,"M")};
    if(sum(grepl("DAF",infiles))>0){sfstype=c(sfstype,"D")}
    
    # Set some plot settings
    # library(RColorBrewer)
    # ramp=colorRampPalette(brewer.pal(11,"RdYlBu"))
    ramp99=c("#A50026","#AA0426","#AF0926","#B40E26","#B91326","#BE1826","#C31D26","#C82226","#CD2726","#D22C26","#D73127","#DA372A","#DD3D2D","#E04330","#E34A33","#E65035","#E95638","#EC5C3B","#EF633E","#F26941","#F46F44","#F57647","#F67C4A","#F7834D","#F88A50","#F89053","#F99756","#FA9E59","#FBA45C","#FCAB5F","#FDB163","#FDB668","#FDBB6D","#FDC072","#FDC577","#FDCA7B","#FDCF80","#FDD485","#FDD98A","#FDDE8F","#FEE293","#FEE598","#FEE89D","#FEECA2","#FEEFA7","#FEF2AB","#FEF5B0","#FEF8B5","#FEFBBA","#FFFEBE","#FBFDC4","#F8FCCA","#F5FBD0","#F2FAD6","#EFF8DC","#ECF7E1","#E8F6E7","#E5F5ED","#E2F3F3","#DEF2F7","#D9EFF6","#D4EDF4","#CEEAF3","#C9E7F1","#C3E5F0","#BEE2EE","#B9DFEC","#B3DDEB","#AEDAE9","#A8D7E8","#A3D2E5","#9DCEE3","#97C9E0","#92C5DE","#8CC0DB","#87BCD9","#81B7D6","#7BB3D4","#76AED1","#71A9CF","#6CA3CC","#679EC9","#6298C6","#5D92C3","#598DC0","#5487BD","#4F81BA","#4A7BB7","#4576B4","#436FB1","#4169AE","#3F63AB","#3D5CA7","#3B56A4","#394FA1","#37499E","#35429B","#333C98","#313695")
    ramp30=c("#A50026","#B61026","#C72126","#D83227","#E24731","#EC5C3B","#F47145","#F7874F","#FA9E59","#FDB365","#FDC476","#FDD586","#FEE496","#FEEEA6","#FEF9B6","#F9FCC8","#EEF8DC","#E4F4F0","#D5EDF4","#C2E4EF","#B0DBEA","#9DCEE3","#8ABFDA","#77B0D2","#679DC9","#568ABF","#4676B5","#3E61AA","#374B9F","#313695")
    
    # Step 4.4: Loop over SFS-type (minor / derived)
    for(fold in sfstype){
      
      # Step 4.5: Print PDF
      pdf(paste(opt$infile,"_2D",fold,"AFs.pdf",sep=""),height=5.3,width=6,useDingbats=F,pointsize=10)
      {
        par(mfrow=c(2,2),oma=c(1,1,0.5,1.5),mar=c(2.5,2.5,2.5,4),mgp=c(1.5,0.5,0))
        
        # Step 4.5: Loop over populations
        for(j in pop){
          
          # Load observed 2D-SFS if present
          if(file.exists(paste(opt$infile,"_joint",fold,"AFpop",j,".obs",sep=""))){
            obs<-as.matrix(read.table(paste(opt$infile,"_joint",fold,"AFpop",j,".obs",sep=""),skip=1,header=T))
            
            # Remove monomorphic sites and singletons for plotting if chosen as option
            # WARNING: these are informative in marginal 2D-SFS for the model fit
            if(opt$nozeros){
              obs[1,1]=0
            }
            if(opt$nosingletons){
              obs[2,1]=0
              obs[1,2]=0
            }
          }else{obs=NA}
          
          # Load expected 2D-SFS if present
          if(file.exists(paste(opt$infile,"_joint",fold,"AFpop",j,".txt",sep=""))){
            # Check if expected SFS has header line
            if(sum(grepl("observation",scan(paste(opt$infile,"_joint",fold,"AFpop",j,".txt",sep=""),what = "character")))>0){
              sim<-as.matrix(read.table(paste(opt$infile,"_joint",fold,"AFpop",j,".txt",sep=""),skip=1,header=T))
            }else{
              sim<-as.matrix(read.table(paste(opt$infile,"_joint",fold,"AFpop",j,".txt",sep=""),skip=0,header=T))
            }

            # Remove monomorphic sites and singletons for plotting if chosen as option
            # WARNING: these are informative in marginal 1D-SFS for the model fit
            if(opt$nozeros){
              sim[1,1]=0
              sim=sim/sum(sim)
            }
            if(opt$nosingletons){
              sim[1,2]=0
              sim[2,1]=0
              sim=sim/sum(sim)
            }
          }else{sim=NA}

          # If both obs and sim are present,
          # multiply probabilities of expected 2D-SFS with number of observed SNPs
          if(is.matrix(obs) & is.matrix(sim)){
            sim<-sum(obs,na.rm=T)*sim
          }
          
          # Get upper and lower bounds for logarithmic color scale and get breaks
          uplog<-log10(max(c(obs,sim),na.rm = T)+0.00001)
          br<-c(10^(seq(log10(0.5),uplog,length.out=100)))
          
          # Print observed SFS
          if(is.matrix(obs)){
            image(0:(dim(obs)[1]-1),0:(dim(obs)[2]-1),obs,col=rev(ramp99),zlim=0.5,breaks=br,axes=F,
                  xlab="",ylab="",main="observed SFS",xpd=NA);box()
            axis(2);axis(1)
          }else{plot(0,0,type="n",xlab="",ylab="",axes=F)}
          mtext(paste("population",unlist(strsplit(j,"[_]"))[2]),2,at=grconvertY(0.5,"ndc","user"),line=2)
          
          # Print expected SFS
          if(is.matrix(sim)){
            image(0:(dim(sim)[1]-1),0:(dim(sim)[2]-1),sim,col=rev(ramp99),zlim=0.5,breaks=br,axes=F,
                xlab="",ylab="",main="expected SFS"); box()
            axis(2);axis(1)
          }else{plot(0,0,type="n",xlab="",ylab="",axes=F)}
          # Print color code legend
          logax=c(1,2,5)*rep(10^(seq(0,floor(uplog))),each=3);
          logax=c(0.5,logax[logax<10^(uplog)],round(10^(uplog),0))
          axis(4,at=grconvertY(0,"npc","user")+diff(grconvertY(0:1,"npc","user"))*
                 (log10(logax)-log10(logax)[1])/diff(log10(logax[c(1,length(logax))])),
               round(logax,1),las=2,line=1.4)
          rasterImage(as.raster(matrix(ramp30, ncol=1)),
                      par('usr')[2]+0.5*diff(grconvertX(0:1, 'inches', 'user'))*par('cin')[2]*par('cex')*par('lheight'),
                      grconvertY(0,"npc","user"),
                      par('usr')[2]+1.5*diff(grconvertX(0:1, 'inches', 'user'))*par('cin')[2]*par('cex')*par('lheight'),
                      grconvertY(1,"npc","user"),xpd=NA)
          mtext(expression(N[loci]),3,at=par('usr')[2]+
                  0.5*diff(grconvertX(0:1,'inches','user'))*par('cin')[2]*par('cex')*par('lheight'),adj=0,line=0.5)

          # Print difference between observed and expected
          if(is.matrix(sim) & is.matrix(obs)){
            dif=sim-obs
            # difcols=colorRampPalette(c(brewer.pal(10,"RdBu")[1:5],"#FFFFFF",brewer.pal(10,"RdBu")[6:10]))
            difcols99=c("#67001F","#6E0220","#760421","#7D0722","#850923","#8D0C25","#940E26","#9C1127","#A41328","#AB162A","#B2192B","#B6202F","#BA2832","#BD2F36","#C13639","#C53E3D","#C84540","#CC4C43","#D05447","#D35B4A","#D7624F","#DA6954","#DD7059","#E0775F","#E37E64","#E6866A","#E98D6F","#EC9475","#EF9B7A","#F2A27F","#F4A886","#F5AD8D","#F6B394","#F7B89B","#F8BEA2","#F9C3A9","#FAC9B0","#FACEB7","#FBD4BE","#FCD9C5","#FDDDCB","#FDE1D1","#FDE5D6","#FDE8DC","#FDECE2","#FEF0E8","#FEF3ED","#FEF7F3","#FEFBF9","#FFFEFE","#FAFCFD","#F5F9FB","#F0F7FA","#ECF4F8","#E7F1F7","#E2EFF5","#DEECF4","#D9E9F2","#D4E7F1","#CFE4EF","#C9E1ED","#C2DDEB","#BCDAEA","#B6D7E8","#AFD4E6","#A9D0E4","#A2CDE2","#9CCAE0","#95C6DF","#8EC2DC","#86BDDA","#7EB8D7","#76B3D4","#6EAED1","#66A9CF","#5EA4CC","#569FC9","#4E9AC6","#4695C4","#4090C1","#3D8BBF","#3987BC","#3682BA","#337DB8","#2F79B5","#2C74B3","#2870B1","#256BAE","#2166AC","#1E61A5","#1B5C9E","#195696","#16518E","#134B87","#10467F","#0D4077","#0A3B70","#073568","#053061")
            difcols30=c("#67001F","#800823","#9A1027","#B31A2C","#BF3337","#CC4C43","#D86450","#E27C62","#EC9475","#F4AA89","#F8BDA0","#FBCFB8","#FDDFCE","#FDECE2","#FEF8F5","#F7FAFC","#E7F1F7","#D7E8F2","#C3DEEC","#AED3E6","#98C8DF","#7EB8D7","#63A7CE","#4896C4","#3986BC","#2D77B4","#2267AC","#185594","#0E427A","#053061")

            uplogd=log10(max(c(dif),na.rm = T)+0.00001)
            brd<-c(10^(seq(log10(0.5),uplog,length.out=50)))
            brd=c(-rev(brd),brd)
            image(0:(dim(dif)[1]-1),0:(dim(dif)[2]-1),dif,col=difcols99,
                  zlim=c(-1,1)*max(abs(dif)),breaks=brd,axes=F,
                  xlab="",ylab="",main="expected - observed SFS",xpd=NA)
            axis(1);axis(2);box()

            logax=c(1,3)*rep(10^(seq(0,floor(uplogd))),each=2)
            axis(4,at=grconvertY(0,"npc","user")+diff(grconvertY(0:1,"npc","user"))*
                   (c(-uplogd,-rev(log10(logax)),log10(logax)[-1],uplogd)+uplogd)/(2*uplogd),
                 c(-round(10^uplogd),-rev(logax[-1]),0,logax[-1],round(10^uplogd)),las=2,line=1.4)
            rasterImage(as.raster(matrix(rev(difcols30), ncol=1)),
                        par('usr')[2]+0.5*diff(grconvertX(0:1, 'inches', 'user'))*par('cin')[2]*par('cex')*par('lheight'),
                        grconvertY(0,"npc","user"),
                        par('usr')[2]+1.5*diff(grconvertX(0:1, 'inches', 'user'))*par('cin')[2]*par('cex')*par('lheight'),
                        grconvertY(1,"npc","user"),xpd=NA)
            mtext(expression(N[loci]),3,at=par('usr')[2]+
                    0.5*diff(grconvertX(0:1,'inches','user'))*par('cin')[2]*par('cex')*par('lheight'),adj=0,line=0.5)
            
            # Relative: Percentage
            dif=sim/obs
            # difcols=colorRampPalette(c(brewer.pal(10,"RdBu")[1:5],"#FFFFFF",brewer.pal(10,"RdBu")[6:10]))
            difcols99=c("#67001F","#6E0220","#760421","#7D0722","#850923","#8D0C25","#940E26","#9C1127","#A41328","#AB162A","#B2192B","#B6202F","#BA2832","#BD2F36","#C13639","#C53E3D","#C84540","#CC4C43","#D05447","#D35B4A","#D7624F","#DA6954","#DD7059","#E0775F","#E37E64","#E6866A","#E98D6F","#EC9475","#EF9B7A","#F2A27F","#F4A886","#F5AD8D","#F6B394","#F7B89B","#F8BEA2","#F9C3A9","#FAC9B0","#FACEB7","#FBD4BE","#FCD9C5","#FDDDCB","#FDE1D1","#FDE5D6","#FDE8DC","#FDECE2","#FEF0E8","#FEF3ED","#FEF7F3","#FEFBF9","#FFFEFE","#FAFCFD","#F5F9FB","#F0F7FA","#ECF4F8","#E7F1F7","#E2EFF5","#DEECF4","#D9E9F2","#D4E7F1","#CFE4EF","#C9E1ED","#C2DDEB","#BCDAEA","#B6D7E8","#AFD4E6","#A9D0E4","#A2CDE2","#9CCAE0","#95C6DF","#8EC2DC","#86BDDA","#7EB8D7","#76B3D4","#6EAED1","#66A9CF","#5EA4CC","#569FC9","#4E9AC6","#4695C4","#4090C1","#3D8BBF","#3987BC","#3682BA","#337DB8","#2F79B5","#2C74B3","#2870B1","#256BAE","#2166AC","#1E61A5","#1B5C9E","#195696","#16518E","#134B87","#10467F","#0D4077","#0A3B70","#073568","#053061")
            difcols30=c("#67001F","#800823","#9A1027","#B31A2C","#BF3337","#CC4C43","#D86450","#E27C62","#EC9475","#F4AA89","#F8BDA0","#FBCFB8","#FDDFCE","#FDECE2","#FEF8F5","#F7FAFC","#E7F1F7","#D7E8F2","#C3DEEC","#AED3E6","#98C8DF","#7EB8D7","#63A7CE","#4896C4","#3986BC","#2D77B4","#2267AC","#185594","#0E427A","#053061")
            
            up=max(c(dif)[!is.infinite(c(dif))],na.rm = T)
            lo=min(c(dif)[!is.infinite(c(dif))],na.rm = T)
            brd=c(seq(lo,1,length.out=50),seq(1,up,length.out=50))
            image(0:(dim(dif)[1]-1),0:(dim(dif)[2]-1),dif,col=difcols99,
                  zlim=c(-1,1)*max(abs(dif)),breaks=brd,axes=F,
                  xlab="",ylab="",main="expected / observed SFS",xpd=NA)
            axis(1);axis(2);box()
            
            ax=c(ceiling(seq(lo*10,10,length.out=3))[-3]/10,1,floor(seq(10,up*10,length.out=3))[-1]/10)
            axis(4,at=grconvertY(0,"npc","user")+diff(grconvertY(0:1,"npc","user"))*
                   seq(0,1,length.out=5),
                 ax*100,las=2,line=1.4)
            rasterImage(as.raster(matrix(rev(difcols30), ncol=1)),
                        par('usr')[2]+0.5*diff(grconvertX(0:1, 'inches', 'user'))*par('cin')[2]*par('cex')*par('lheight'),
                        grconvertY(0,"npc","user"),
                        par('usr')[2]+1.5*diff(grconvertX(0:1, 'inches', 'user'))*par('cin')[2]*par('cex')*par('lheight'),
                        grconvertY(1,"npc","user"),xpd=NA)
            mtext("%",3,at=par('usr')[2]+
                    0.5*diff(grconvertX(0:1,'inches','user'))*par('cin')[2]*par('cex')*par('lheight'),adj=0,line=0.5)
          }else{plot(0,0,type="n",xlab="",ylab="",axes=F)}
          
          mtext(paste("population",unlist(strsplit(j,"[_]"))[1]),1,at=grconvertX(0.5,"ndc","user"),line=2)
        }
      }
      dev.off()
    }
  }else{
    stop(paste("No two-dimenstional SFS files (",opt$infile,"_joint[M/D]AFX_Y.[obs/txt]) found in this directory!",sep=""))
  }
}else{
  stop(paste("Aborted. Incorrect tools specification: ",opt$tool,
             ". Provide one of the following: mto2D, mto1D, print1D, print2D",sep=""))
}


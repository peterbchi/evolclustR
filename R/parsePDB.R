# 10/24/18 still in progress

parsePDB<-function(file, subunit = NULL){
#  wholefile<-strsplit(readLines(file),split=" ")	
#  structure<-strsplit(readLines("2i0q.pdb")[1061:4811],split=" ")

  wholefile <- strsplit(readLines(file),split=" ")
  structure <- wholefile[grep("ATOM", sapply(wholefile,function(x) x[1]))]

  Calphas<-NA
  for(i in 1:length(structure)){
     structure[[i]]<-structure[[i]][structure[[i]]!=""]
     if(structure[[i]][3]=="CA"){
         Calphas<-rbind(Calphas,structure[[i]])
     }
  }

  Calphas<-Calphas[-1,]
  Calphas<-data.frame(Calphas,stringsAsFactors=FALSE)
  Calphas[,c(2,6,7,8,9,10,11)]<-lapply(Calphas[,c(2,6,7,8,9,10,11)], as.numeric)	# 4/11/19 in the future, make sure this hard coding isn't a bad idea
  colnames(Calphas)<-c("type","atompos","atom","AA","nn","pos","x","y","z","w","ww","www")

  # Keep only subunit of interest
  if(!is.null(subunit)){
    Calphas <- Calphas[Calphas$nn == subunit,]
  }
  
  # Remove duplicates, keep first one in all cases
  i <- 1
  while(i < dim(Calphas)[1] - 1){
    if(Calphas$pos[i] == Calphas$pos[i+1]){
      Calphas <- Calphas[-(i+1),]
    } else {
      i <- i+1
    }
  }

  # rename amino acids that had duplicates (need to shave off leading "A")
  Calphas$AA <- ifelse(nchar(Calphas$AA) == 3, Calphas$AA, substring(Calphas$AA,2))

  # 12/5/18 still need to rename stuff?


  # remove sites with two CA atoms (and two of everything, really)
#  dup.sites<-as.numeric(names(which(table(Calphas$pos)==2)))   	# duplicate sites
#  ind.dup<-which(is.element(Calphas$pos,dup.sites))		# index of duplicate sites

  # one of them even has 3
#  tri.sites<-as.numeric(names(which(table(Calphas$pos)==3)))   	# duplicate sites
#  ind.tri<-which(is.element(Calphas$pos,tri.sites))		# index of duplicate sites


#  Calphas<-Calphas[-c(ind.dup[c(2,4,6,8,10,12,14,16,18,20,22)],ind.tri[c(2,3)]),]


  return(Calphas)
}

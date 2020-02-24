#' An evolClustR Function
#'
#' Function to generate one substitution in the background scheme
#' @param Cdata.gam PDB file with sasa and gamma values concatenated (run calcsasa and calcgamma first), and column for subbed sites
#' @param ratio the relative ratio of sasa vs. distance probabilities
#' @param b.length branch length (expected number of substitutions)	# actually this should go outside this function 9/12/18
#' @param sd.n.subs standard deviation for number of substitutions	# don't need this either
#' @param exponent power to raise distance by in denominator
#' @keywords cluster
#' @export
#' @examples
#' generate.one.sub()

# The purpose of this function is to choose one substitution, 
# with relative probabilities based on SASA and distance combined;
# P(A or B) = P(A) + P(B) - P(A and B); assume independence


# 4/9/19 It runs but need to double-check that it's doing the right thing.
# 7/15/19: Why does it take in subs as a vector of 0s and 1s, but output subs as position? CHANGE THIS.

generate.sasa.subs <- function(Cdata, b.length, ratio, exponent = 6, subs = NULL){
  # 4 categories, alpha=0.8, beta=1, scaled so the mean=1
  scaled.dgam <- qgamma(c(0.125, 0.375, 0.625, 0.875), 0.8, 1)/mean(qgamma(c(0.125, 0.375, 0.625, 0.875), 0.8, 1))

  sasa.cats <- cut(Cdata$acc.ratio, summary(Cdata$acc.ratio)[c(1,2,3,5,6)],include.lowest=T) 
  sasa.ind <- as.numeric(sasa.cats)
  scale.factor <- scaled.dgam[sasa.ind]


  # first round
  n<-dim(Cdata)[1]
  sasa.prob<-scale.factor / sum(scale.factor)
  
  # do multinomial 2/7/19
  # only do this if not provided any starting substitutions
  if(is.null(subs)){
    subs <- as.vector(rmultinom(1,1,sasa.prob))  # these will be row indices (not position)
  }
  
  start <- sum(subs) + 1	# find how many subs exist, start 1 after that

  # fill in rest of subs
  for(i in start:b.length){
    
    # first find min distances from existing subs (update every time)
    dists <- matrix(NA,ncol=sum(subs),nrow=dim(Cdata)[1])
    j <- 0
    for(s in which(subs==1)){		# 7/19/18: changed index from i to s. Can't be i b/c permutation.  
      j <- j+1
      dists[,j]<-find.dist(Cdata,s)    
    }
  
    min.dist <- apply(dists,1,min)
    dist.prob <- ifelse(min.dist==0,0,1/min.dist^exponent)
    dist.prob <- dist.prob/sum(dist.prob)     # normalize to sum to 1 (do i need this? 4/9/19)

    # use SASA and distance to calculate new probability
    new.prob <- sasa.prob*ratio + dist.prob - sasa.prob*ratio*dist.prob
    new.prob[which(subs==1)] <- 0
    new.prob <- new.prob/sum(new.prob)
    subs <- subs + as.vector(rmultinom(1,1,new.prob)) # both are vectors
    
  }
    
  return(Cdata[which(subs==1),]$pos)  # return position number of substituted sites
}


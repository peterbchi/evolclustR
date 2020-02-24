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
generate.sasa.null <- function(Cdata, b.length, ratio, exponent = 6){
  # first round
  n<-dim(Cdata)[1]
  sasa.prob<-Cdata$scale.factor / sum(Cdata$scale.factor)
  
  # do multinomial 2/7/19
  subs <- as.vector(rmultinom(1,1,sasa.prob))  # these will be row indices (not position)
  
  # fill in rest of subs
  for(i in 2:b.length){
    
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
    
  return(subs)
}


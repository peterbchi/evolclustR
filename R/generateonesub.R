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

# add probabilities from each process and subtract the product
# scale before and after
# 
generate.one.sub <- function(Cdata, ratio, exponent = 6, curr_subs){
  if(is.na(curr_subs)[1]){
    
  } else {
    # Find distance from each existing sub
    dists <- matrix(NA,ncol=sum(first.round),nrow=dim(Cdata)[1])
    j <- 0
    for(s in which(first.round)){		# 7/19/18: changed index from i to s. Can't be i b/c permutation.
      j <- j+1
      dists[,j]<-find.dist(Cdata,s)    
    }
    
    min.dist <- apply(dists,1,min)
    temp.prob <- ifelse(min.dist==0,0,1/min.dist^exponent)
    
  }
}



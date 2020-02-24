#' An evolClustR Function
#'
#' Function to generate subs used in experimental simulations, clustered by 1/r^6
#' @param Cdata.gam PDB file with sasa and gamma values concatenated (run calcsasa and calcgamma first)
#' @param scaling scaling factor for relative probabilities of substitution
#' @param exp.n.subs total expected number of substitutions
#' @param sd.n.subs standard deviation for number of substitutions
#' @param exponent power to raise distance by in denominator
#' @keywords cluster
#' @export
#' @examples
#' generate.subs()


generate.subs <- function(Cdata, scaling, exp.n.subs, sd.n.subs, exponent = 6){
  # pick substitutions
  # first round
  n<-dim(Cdata)[1]
  prob<-Cdata$scale.factor * scaling
  
  # make sure to get at least one
  first.round <- 0
  while(sum(first.round)==0){
    first.round <- rbernoulli(n, prob)   # these are row indexes (NOT site position). Gets converted below.
  }
  
  Cdata <- Cdata %>% 
    mutate(subbed = first.round)
  
  
  # find distance from each sub
  # 3/16/18: fix the finddist.R file (7/10/18: not sure if this has been done yet or not)
  dists <- matrix(NA,ncol=sum(first.round),nrow=dim(Cdata)[1])
  j <- 0
  for(s in which(first.round)){		# 7/19/18: changed index from i to s. Can't be i b/c permutation.
    j <- j+1
    dists[,j]<-find.dist(Cdata,s)    
  }
  
  min.dist <- apply(dists,1,min)
  temp.prob <- ifelse(min.dist==0,0,1/min.dist^exponent)

  num.subs <- round(rnorm(n=1, mean=exp.n.subs, sd=sd.n.subs),0)
  while(sum(first.round) < num.subs){
    dists <- matrix(NA,ncol=sum(first.round),nrow=dim(Cdata)[1])
    j <- 0
    for(s in which(first.round)){
      j <- j+1
      dists[,j]<-find.dist(Cdata,s)    
    }
    
    min.dist <- apply(dists,1,min)  # recalculate distance on each iteration
    temp.prob <- ifelse(min.dist==0,0,1/min.dist^exponent)
    rel.prob <- temp.prob/sum(temp.prob)
    first.round[sample(seq(1:n),1,prob=rel.prob)] <- TRUE
  }

  subs <- Cdata[first.round,6] # this needs to be position number to match previous code

  return(subs)			# AGAIN NOTE, THIS IS POSITION NUMBER, NOT LINE NUMBER
}



#' An evolClustR Function
#'
#' Function to return fractions in each cluster for any given radius
#' @param Cdata PDB file parsed to contain only central Carbon atoms
#' @param num.neighbors number of neighbors for every site
#' @param subs substituted sites (by position number)
#' @param radius Radius size of cluster of interest
#' @param contact TRUE/FALSE whether permutation should control for contact number
#' @param similarity If contact is TRUE, how similar does contact number have to be
#' @keywords cluster
#' @export
#' @examples
#' permuteone()


permute.one <- function(Cdata, subs, radius, num.neighbors=NA, contact=FALSE, similarity=0){
  if(contact & is.null(num.neighbors)){
    stop("If contact=TRUE, num.neighbors must be provided (usually calculated from parent function permute or permute.test)")
  }
	  
  n.subs <- length(subs)

  if(contact){
    perm.output<-rep(NA, n.subs)

    # Pull out num.neighbors for subbed sites. 'sample' to randomize order
    sub.neighbors <- sample(num.neighbors[which(is.element(Cdata$pos, subs))])
  
    j <- 1 # need manual index since i is doing cray cray stuff

    for(i in sub.neighbors){
      choices <- which(num.neighbors<=i+similarity & num.neighbors >=i-similarity)

      # remove from consideration any sites that are already subbed
      choices <- choices[!is.element(choices, perm.output)]

      temp.sim <- similarity # kind of dumb since it won't get used if we don't enter the while loop, but I can't think of a better way
      while(length(choices) == 0){ 
	temp.sim <- temp.sim + 1
        choices <- which(num.neighbors<=i+temp.sim & num.neighbors >=i-temp.sim)
        choices <- choices[!is.element(choices, perm.output)]
        
        warning("At least one permutation attempt had zero choices; similarity increased as needed.
		Consider re-running with higher similarity value.")
      }

      perm.output[j] <- ifelse(length(choices)>1, sample(choices, 1), choices)
      j <- j+1
    }
    return(Cdata$pos[perm.output])
  } else {
    return(sample(Cdata$pos, n.subs))
  }
}



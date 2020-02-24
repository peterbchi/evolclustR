#' An evolClustR Function
#'
#' Function to do permutation test
#' @param Cdata PDB file parsed to contain only central Carbon atoms
#' @param subs substituted sites (by position number)
#' @param reps number of replicates to produce
#' @param radius Radius size of cluster of interest
#' @param contact TRUE/FALSE whether permutation should control for contact number
#' @param similarity If contact is TRUE, how similar does contact number have to be
#' @keywords cluster
#' @export
#' @examples
#' permute()


permute <- function(Cdata, subs, reps, radius, contact=FALSE, similarity=0){
  if(contact){
    num.neighbors <- count.neighbors(Cdata,radius)
  } else {
    num.neighbors <- NA
  }

  permuted.data <- replicate(reps, permute.one(Cdata, subs, radius, num.neighbors, contact, similarity))

  # ultimately, find pvalues and return
  # for now, let script file do that while testing, and just return 
  return(permuted.data)  # rows are iterations, columns are each sub

}



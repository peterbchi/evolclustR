#' An evolClustR Function
#'
#' Function to do permutation test
#' @param Cdata PDB file parsed to contain only central Carbon atoms
#' @param subs substituted sites (by position number)
#' @param reps number of replicates to produce
#' @param radius Radius size of cluster of interest
#' @param contact TRUE/FALSE whether permutation should control for contact number
#' @param similarity If contact is TRUE, how similar does contact number have to be
#' @param test Which test to perform, default to 6=95th percentile, subs only
#' @keywords cluster
#' @export
#' @examples
#' permute()


# 8/13/19 hasn't been tested that thoroughly yet
permute.test <- function(Cdata, subs, reps, radius, contact=FALSE, similarity=0, test=6, fraction=TRUE){
  if(contact){
    num.neighbors <- count.neighbors(Cdata,radius)
  } else {
    num.neighbors <- NA
  }

  permuted.data <- replicate(reps, permute.one(Cdata, subs, radius, num.neighbors, contact, similarity))

  if(test == 6 | test == 7){
    perm.fractions <- apply(permuted.data, 2, fractions_subsonly, Cdata = Cdata, radius = radius, fraction=fraction)
    data.fractions <- fractions_subsonly(Cdata, subs, radius, fraction=fraction)
  } else {
    perm.fractions <- apply(permuted.data, 2, fractions, Cdata = Cdata, radius = radius, fraction=fraction)
    data.fractions <- fractions(Cdata, subs, radius, fraction=fraction)
  }

  # statistic: mean
  if(test==1 | test == 7){
    perm.mean <- apply(perm.fractions, 2, mean)
    data.mean <- mean(data.fractions)
    return(sum(perm.mean >= data.mean)/length(perm.mean))
  }

  # statistic: sd
  if(test==2){
    perm.sd <- apply(perm.fractions, 2, sd)
    data.sd <- sd(data.fractions)
    return(sum(perm.sd >= data.sd)/length(perm.sd))
  }

  # ks-test
  if(test == 3){
    perm.all <- as.vector(perm.fractions)
    return(ks.test(data.fractions, perm.all, alternative="less")$p.value)
  }

  # tail.prop.test
  if(test == 4){
    perm.all <- as.vector(perm.fractions)
    return(tail.prop.test(data.fractions, perm.all, 0.05))
  }

  # 95th quantile test
  if(test == 5 | test == 6){
    perm.95th <- apply(perm.fractions, 2, quantile, probs = 0.95, type=7)
    data.95th <- quantile(data.fractions, 0.95, type=7)
    return(sum(perm.95th >= data.95th)/length(perm.95th))
  }

}



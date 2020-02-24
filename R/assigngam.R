#' An evolClustR Function
#'
#' Function to assign gamma rates based on accessibility ratio
#' @param Cdata.sasa output from calc.sasa
#' @param n.cats number of discrete gamma categories
#' @param alpha alpha parameter for gamma distribution
#' @param beta beta parameter for gamma distribution
#' @keywords cluster
#' @export
#' @examples
#' assign.gam()

assign.gam <- function(Cdata.sasa, n.cats, alpha, beta){
   medians <- seq(1:n.cats)/n.cats - 0.5/n.cats		# find the percentile of the midpoint of each category
   scaled.dgam <- qgamma(medians, 0.8, 1)/mean(qgamma(medians, 0.8, 1))

  Cdata.gam <- Cdata.sasa %>% 
    mutate(sasa.cats = cut(acc.ratio, summary(acc.ratio)[c(1,2,3,5,6)],include.lowest=T)) %>% 
    mutate(sasa.ind = as.numeric(sasa.cats)) %>% 
    mutate(scale.factor = scaled.dgam[sasa.ind])

    return(Cdata.gam)

}

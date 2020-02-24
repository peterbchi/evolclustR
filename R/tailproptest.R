' An evolClustR Function
#'
#' Function to do proportions test on tail
#' @param nulldist Proportions in the null distribution
#' @param testdist Proportions in the test distribution
#' @param tailprop 
#' @keywords cluster
#' @export
#' @examples
#' tailproptest()

tail.prop.test<-function(testdist, nulldist, tailprop){
  count <- sum(testdist > quantile(nulldist, (1-tailprop)))
  p.hat <- count/length(testdist)
  pooled <- (count + tailprop*length(nulldist)) / (length(nulldist) + length(testdist))

  test.stat <- (p.hat - tailprop)/sqrt(pooled*(1-pooled) * (1/length(nulldist) + 1/length(testdist)))

  pvalue <- pnorm(test.stat, lower.tail = FALSE)

  return(pvalue)
}



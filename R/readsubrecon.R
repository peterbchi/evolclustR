#' An evolClustR Function
#'
#' Function to read substitutions from SubRecon output
#' @param file.name file containing verbose output from SubRecon
#' @keywords cluster
#' @export
#' @examples
#' read.subs(file.name)

read.subs <- function(file.name, blank.as.change = FALSE){
  wholefile <- strsplit(readLines(file.name),split=" ")
  results <- wholefile[grep("Result", sapply(wholefile,function(x) x[1]))]

  output <- NULL
  for(i in 1:length(results)){
    tmp <- unlist(strsplit(results[[i]], "\t"))
    if(blank.as.change){
      if(length(tmp)<4){
        output <- c(output, i)
      } else{
        tmp2 <- unlist(strsplit(tmp[4],""))
        if(tmp2[1] != tmp2[2]){
          output <- c(output, i)
        }
      }
    } else {
      if(length(tmp)==4){
        tmp2 <- unlist(strsplit(tmp[4],""))
        if(tmp2[1] != tmp2[2]){
          output <- c(output, i)
        }
      }
    }
  }
  return(output)
}

#' An evolClustR Function
#'
#' Function to read node [B] sequence from SubRecon output
#' @param file.name file containing verbose output from SubRecon
#' @keywords cluster
#' @export
#' @examples
#' get.sequence(file.name)

get.sequence <- function(file.name){
  wholefile <- strsplit(readLines(file.name),split=" ")
  results <- wholefile[grep("Result", sapply(wholefile,function(x) x[1]))]

  output <- NULL
  for(i in 1:length(results)){
    tmp <- unlist(strsplit(results[[i]], "\t"))
    if(length(tmp)<4){
      output <- c(output, NA)
    } else{
      tmp2 <- unlist(strsplit(tmp[4],""))
      output <- c(output, tmp2[2])
    }
  }
  return(output)
}

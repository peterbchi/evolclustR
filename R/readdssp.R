#' An evolClustR Function
#'
#' Function to read output file from dssp
#' @param dssp.filename name of dssp file
#' @keywords cluster
#' @export
#' @examples
#' read.dssp()

read.dssp <- function(dssp.filename){
  dssp.temp <- readLines(dssp.filename)
  h.ind <- grep(" ACC ", dssp.temp)   # find index of header
  dssp.temp <- dssp.temp[-c(1:h.ind)] # get rid of everything header and above
  
  n <- length(dssp.temp)      
  residue <- rep(NA,n)
  subunit <- rep(NA,n)
  sasa <- rep(NA, n)
  for(i in 1:length(dssp.temp)){
      residue[i] <- as.numeric(
        paste(
          unlist(
            strsplit(
              dssp.temp[i], 
              "")
          )[c(8:10)],sep="",collapse=""
        )
      )
      
      subunit[i] <- unlist(strsplit(dssp.temp[i],""))[12]
      
      sasa[i] <- as.numeric(
      paste(
        unlist(
          strsplit(
            dssp.temp[i], 
            "")
          )[c(36:38)],sep="",collapse=""
        )
      )
  }
  
  return(data.frame(residue,subunit,sasa,stringsAsFactors = F))
}


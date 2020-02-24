#' An evolClustR Function
#'
#' Function to calculate sasa (ratio) from dssp output and max allowable accessibility
#' @param Cdata PDB file parsed to contain only central Carbon atoms
#' @param dssp.output output from read.dssp
#' @keywords cluster
#' @export
#' @examples
#' calc.sasa()

calc.sasa <- function(Cdata, dssp.output){

  residue_max_acc <- read.csv("acc.txt")

  # original values from Claudia.
  # Above was typed in manually based on this; keep this
  # for posterity.
  #  {'A': 129.0, 'R': 274.0, 'N': 195.0, 'D': 193.0, \
  #  'C': 167.0, 'Q': 225.0, 'E': 223.0, 'G': 104.0, \
  #  'H': 224.0, 'I': 197.0, 'L': 201.0, 'K': 236.0, \
  #  'M': 224.0, 'F': 240.0, 'P': 159.0, 'S': 155.0, \
  #  'T': 172.0, 'W': 285.0, 'Y': 263.0, 'V': 174.0}

  Cdata.new <- merge(Cdata, residue_max_acc)
  Cdata.new <- merge(Cdata.new,dssp.output, by.x=c("pos","nn"), by.y=c("residue","subunit"))
  Cdata.new <- Cdata.new[order(Cdata.new$pos),]

  Cdata.new$acc.ratio <- Cdata.new$sasa / Cdata.new$aa.acc

  return(Cdata.new)
}

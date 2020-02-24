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

  residue_max_acc <- data_frame(AA=c("ALA","ARG","ASN","ASP",
                                   "CYS","GLN","GLU","GLY",
                                   "HIS","ILE","LEU","LYS",
                                   "MET","PHE","PRO","SER",
                                   "THR","TRP","TYR","VAL"),
                                 aa.acc=c(129, 274, 195, 193,
                                       167, 225, 223, 104,
                                       224, 197, 201, 236,
                                       224, 240, 159, 155,
                                       172, 285, 263, 174))
  # original values from Claudia.
  # Above was typed in manually based on this; keep this
  # for posterity.
  #  {'A': 129.0, 'R': 274.0, 'N': 195.0, 'D': 193.0, \
  #  'C': 167.0, 'Q': 225.0, 'E': 223.0, 'G': 104.0, \
  #  'H': 224.0, 'I': 197.0, 'L': 201.0, 'K': 236.0, \
  #  'M': 224.0, 'F': 240.0, 'P': 159.0, 'S': 155.0, \
  #  'T': 172.0, 'W': 285.0, 'Y': 263.0, 'V': 174.0}


  # shave "A" off of length-4 AA names
  Cdata <- Cdata %>% 
    mutate(AA = ifelse(nchar(AA) == 3, AA, substring(AA,2)))

  Cdata.new <- Cdata %>% 
    left_join(dssp.output, by=c("pos"="residue", "nn"="subunit")) %>% 
    left_join(residue_max_acc, by="AA") 

  Cdata.new <- Cdata.new %>% 
    mutate(acc.ratio = sasa/aa.acc)

  return(Cdata.new)
}

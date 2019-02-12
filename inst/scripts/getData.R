#' getData.R
#'
#' @description
#'
#'
#' @export

getData <- function() {

  message("Reading Ensembl data...")
  martFile <- file.path("../data", "mart_export.txt")

  martDF <- read.csv(martFile, stringsAsFactors = FALSE)

  # Remove any sapien genes that do not have a HGNC symbol
  # or PDB entry
  martDF <- martDF[(("" != martDF$HGNC.symbol) & ("" != martDF$HGNC.transcript.name.ID)
                        & ("" != martDF$PDB.ENSP.mappings)),]

  return(martDF)

}

#[END]

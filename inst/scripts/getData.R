#' getData.R
#'
#' \code{getData} Retrieve the data from the mart_export.txt file
#' from BiomaRts webpage and download services. Guidelines
#' provided in the readMe. \code{getData} returns a
#' dataframe containing the biomaRt data.
#'
#' Details.
#'
#' @return martDF; a dataframe containing the biomaRt data
#' for homosapien genes
#'

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

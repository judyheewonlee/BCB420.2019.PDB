#' addpdbChain.R
#'
#' \code{addpdbChain} Add the PDB chains and the corresponding HGNC ID of transcripts
#' to the database.
#'
#' Details.
#'
#' @param myDB A list. The database which is being built.
#' @param martDF The dataframe containing BiomaRt data.
#' @return myDB with the pdbChains added to the database.
#'

addpdbChain <- function(myDB, martDF) {

  # Add the PDB chains and their corresponding
  # transcript IDs to the database
  pdbChains <- data.frame(ID = martDF$PDB.ENSP.mappings,
                          transcriptHGNC =
                            martDF$HGNC.transcript.name.ID,
                          stringsAsFactors = FALSE)
  pdbChains <- unique(data.table::as.data.table(pdbChains))

  myDB$pdbChains <- as.data.frame(pdbChains)

  return(myDB)
}

# [END]

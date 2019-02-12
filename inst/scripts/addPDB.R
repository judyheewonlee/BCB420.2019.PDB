#' addPDB.R
#'
#' \code{addPDB} Add the PDB IDs and the corresponding PDB chains in that
#' PDB ID to the database.
#'
#' Details.
#'
#' @param myDB A list. The database which is being built.
#' @param martDF The dataframe containing BiomaRt data.
#' @return myDB with the PDB IDs added to the database.
#'

addPDB <- function(myDB, martDF) {

  # Add the PDB IDs and their corresponding
  # PDB chains to the database

  pdbIDs <- data.frame(ID = martDF$PDB.ID,
                      chainID = martDF$PDB.ENSP.mappings,
                      stringsAsFactors = FALSE)

  pdbIDs <- unique(data.table::as.data.table(pdbIDs))

  myDB$pdbID <- as.data.frame(pdbIDs)
  return(myDB)

}

# [END]

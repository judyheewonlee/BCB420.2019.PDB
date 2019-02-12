#' addPDBdata.R
#'
#' \code{addPDBdata} Add the PDB data retrieved from PDB. This include
#' The sequence and Resolution of each PDB chain in the database.
#'
#' Details.
#'
#' @param myDB A list. The database which is being built.
#'
#' @return myDB with the PDB data added to the database.
#'

addPDBdata <- function(myDB) {

  message("Retrieving PDB data...")
  # Call fetchPDBXML to retrieve resolutions and sequences
  # from PDB
  myPDBData <- unique(fetchPDBXML(myDB))

  # Set Resolution and sequence
  myDB$pdbChains$Resolution <- NA
  myDB$pdbChains$Sequences <- NA

  # Create a vector of the resolutions and sequences of
  # each PDB chain and add it to the DB

  for (chain in unique(myDB$pdbChains$ID)) {
    pdbID <- toupper(gsub("\\..*","", chain))

    resolution <- unique(myPDBData[pdbID
                          == myPDBData$IDs,]$Resolution)

    chainSeq <- unique(myPDBData[sub('.*\\.', '', chain) == myPDBData$ChainIDs &
                          pdbID == myPDBData$IDs,]$Sequence)[1]

    myDB$pdbChains[chain == myDB$pdbChains$ID,]$Resolution <- resolution
    myDB$pdbChains[chain == myDB$pdbChains$ID,]$Sequences <- chainSeq
  }

  return(myDB)

}

# [END]


#' PDBdataset.R
#'
#' @description
#'
#' @export
#'
#'

PDBdataset <- function() {

  # Fetch Data from ensembl and retrieve a dataframe
  martDF <- getData()

  message("Building the database...")
  # Build the database
  myDB <- list()
  myDB <- addHGNC(myDB, martDF)
  myDB <- addTranscripts(myDB, martDF)

  message("Adding PDB chains...")
  myDB <- addpdbChain(myDB, martDF)

  message("Adding PDB IDS ...")
  myDB <- addPDB(myDB, martDF)
  myDB <- addPDBdata(myDB)

  message("Database successfully generated!")

  return(myDB)

}

# [END]


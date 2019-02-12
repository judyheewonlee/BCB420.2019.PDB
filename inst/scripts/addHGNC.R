#' addHGNC.R
#'
#' \code{addHGNC} Add the HGNC symbols to the data set along
#' with the corresponding gene description, stable ID,
#' stable version ID and gene name.
#'
#' Details.
#'
#' @param myDB A list. The database which is being built.
#' @param martDF The dataframe containing BiomaRt data.
#' @return myDB with the HGNC data added to the database.
#'

addHGNC <- function(myDB, martDF) {

  geneRows <- martDF[!duplicated(martDF$HGNC.symbol),]
  myDB$HGNC <- data.frame(ID = geneRows$HGNC.symbol,
                          description = geneRows$Gene.description,
                          stableID = geneRows$Gene.stable.ID,
                          version = geneRows$Gene.stable.ID.version,
                          name = geneRows$Gene.name,
                          stringsAsFactors = FALSE)

  return (myDB)
}

# [END]

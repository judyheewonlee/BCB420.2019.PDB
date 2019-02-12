#' addHGNC.R
#'
#' @description
#'
#' @param myDB
#'
#' @param geneList
#'
#' @export

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

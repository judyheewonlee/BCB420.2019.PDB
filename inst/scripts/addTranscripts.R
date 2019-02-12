#' addTranscripts.R
#'
#' @description
#'
#' @param myDB
#'
#' @param geneList
#'
#' @export


addTranscripts <- function(myDB, martDF) {

  # Retrieve rows with unduplicated transcript IDs
  geneRows <- martDF[!duplicated(martDF$Transcript.stable.ID),]
  myDB$Transcripts <- data.frame(ID = geneRows$HGNC.transcript.name.ID,
                              version = geneRows$Transcript.stable.ID.version,
                              hgncID = geneRows$HGNC.symbol,
                              stableID = geneRows$Transcript.stable.ID,
                              start = geneRows$Transcript.start..bp,
                              end = geneRows$Transcript.end..bp,
                              stringsAsFactors = FALSE)

  return(myDB)

}

# [END]


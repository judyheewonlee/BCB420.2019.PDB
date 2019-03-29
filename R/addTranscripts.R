#' addTranscripts.R
#'
#' \code{addTranscripts} Add the transcript data to the database.
#' This includes the HGNC ID for the transcript, the stable ID, the
#' stable ID with the version number, The corresponding HGNC symbols,
#' and the start and end gene coordinates.
#'
#' Details.
#'
#' @param myDB A list. The database which is being built.
#' @param martDF The dataframe containing BiomaRt data.
#' @return myDB with the Transcript data added to the database.
#'

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


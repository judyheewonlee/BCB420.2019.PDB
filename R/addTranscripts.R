



addTranscripts <- function(myDB, geneList) {

  # Retrieve rows with unduplicated transcript names
  geneRows <- geneList[!duplicated(geneList[,5]),]

  f <- function(geneRows) {
    x <- data.frame(
      ID = geneRows[[5]],
      hgncID = geneRows[[1]],
      transcriptName = geneRows[[4]],
      stringsAsFactors = FALSE
    )
  }

  myDB$Transcripts <- rbind(myDB$Transcripts,
                         do.call("rbind", apply(geneRows, 1, f)))

  return(myDB)

}

# [END]


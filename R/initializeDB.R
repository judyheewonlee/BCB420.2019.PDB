




initializeDB <- function(myDB) {

  myDB$hgnc <- data.frame(ID = character(),
                          stringsAsFactors = FALSE
                          )

  myDB$Transcripts <- data.frame(
                      ID = character(),
                      hgncID = character(),
                      transcriptName = character(),
                      stringsAsFactors = FALSE
                      )

  myDB$pdbChain <- data.frame(
                  ID = character(),
                  transcriptID = character(),
                  stringsAsFactors = FALSE
                  )

  myDB$pdb <- data.frame(
              ID = character(),
              chainID = character(),
              stringsAsFactors = FALSE
              )

  return (myDB)

}


# [END]

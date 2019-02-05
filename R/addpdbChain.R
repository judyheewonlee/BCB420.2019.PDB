



addpdbChain <- function(myDB, geneList) {

  # Add the PDB chains and their corresponding
  # transcript IDs to the database
  # Note that each transcript has a unique PDB chain
  # but not all PDB chains are exclusive to one transcript
  geneRows <- geneList[!duplicated(geneList[,5]),]

  f <- function(geneRows) {
    x <- data.frame(ID = geneRows[[3]],
                    transcriptID = geneRows[[5]],
                    stringsAsFactors = FALSE)
  }

  myDB$pdbChain <- rbind(myDB$pdbChain,
                         do.call("rbind", apply(geneRows, 1, f)))

  return(myDB)
}

# [END]

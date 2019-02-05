





addPDB <- function(myDB, geneList) {

  # Since each PDB ID can have multiple PDB chain ID's
  # but each PDB chain ID has a unique PDB ID, we remove
  # any duplicated PDB chain ID's and determine the
  # corresponding PDB ID for that chain
  geneRows <- geneList[!duplicated(geneList[,3]),]

  f <- function(geneRows) {
    x <- data.frame(
      ID = geneRows[[2]],
      chainID = geneRows[[3]],
      stringsAsFactors = FALSE
    )
  }

  myDB$pdb <- rbind(myDB$pdb,
                            do.call("rbind", apply(geneRows, 1, f)))

  return(myDB)

}

# [END]

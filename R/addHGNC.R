


addHGNC <- function(myDB, geneList) {

  geneRows <- geneList[!duplicated(geneList[,1]),]

  f <- function(geneRows) {
    x <- data.frame(
      ID = geneRows[[1]],
      stringsAsFactors = FALSE
    )
  }

  myDB$hgnc <- rbind(myDB$hgnc,
                      do.call("rbind", apply(geneRows, 1, f)))

  return (myDB)
}

# [END]

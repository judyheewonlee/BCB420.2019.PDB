#' fetchPDBXML.R
#'
#' \code{fetchPDBXML} Fetch the PDB data including Sequence and Resolution
#' from PDB's API web services and return a dataframe including
#' the PDB ID, chain identifier, Sequence and Resolution for each chain
#' in the database.
#'
#' Details.
#'
#' @param myDB A list. The database which is being built.
#'
#' @return myPDBData; a data frame with the PDB data.
#'

fetchPDBXML <- function(myDB) {

  # Generate a list of PDB chains, remove the
  # chain identifier
  chainList <- unique(myDB$pdbChains$ID)
  chainList <- gsub("\\..*", "", chainList)
  myPDBData <- data.table::data.table()

  # Code to split list referenced by
  # https://stackoverflow.com/questions/7060272/split-up-a-dataframe-by-number-of-rows
  # Split the chainList into sets of 1500 elements
  chunk <- 1500
  n <- length(chainList)
  r  <- rep(1:ceiling(n/chunk),each=chunk)[1:n]
  sptChainList <- split(chainList, r)

  #Fetch the XML data for each chain
  for (i in 1:length(sptChainList)) {

    url <- paste("https://www.rcsb.org/pdb/rest/customReport.xml?pdbids=",
                 paste(sptChainList[[i]], collapse = ","),
                 "&customReportColumns=structureId,resolution,sequence",
                 sep = "")
    myRefData <- readXML(url)
    myPDBData <- rbind(myPDBData, myRefData)
  }

  return (myPDBData)
}

#[END]

#' readXML.R
#'
#' @description Retrieve the XML data from the PDB web services
#' and generate a dataframe containing the resolution
#' and StructureID for a PDB structure of interest
#'
#' @param URL The URL to the XML data from the PDB webservices in the
#' form of a character vector
#'
#' @return The dataframe containing the protein structure ID and
#' the corresponding resolution
#'

readXML <- function(URL) {

  library(xml2)

  myXML <- read_xml(URL)
  xml_name(myXML)
  xml_children(myXML)


  myIDnodes <- xml_find_all(myXML, ".//dimEntity.structureId")
  myResNodes <- xml_find_all(myXML, ".//dimStructure.resolution")
  myChainNodes <- xml_find_all(myXML, ".//dimEntity.chainId")
  mySeqNodes <- xml_find_all(myXML, ".//dimEntity.sequence")

  myPDBData <- data.frame(IDs =
                          paste(tolower(as.character(xml_contents(myIDnodes))),
                                as.character(xml_contents(myChainNodes)),
                                sep = "."),
                          Resolution = as.character(xml_contents(myResNodes)),
                          ChainIDs = as.character(xml_contents(myChainNodes)),
                          Sequence = as.character(xml_contents(mySeqNodes)),
                          stringsAsFactors = FALSE)

  return (myPDBData)

}


# [END]

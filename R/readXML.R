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


  myXML <- xml2::read_xml(URL)
  xml2::xml_name(myXML)
  xml2::xml_children(myXML)


  myIDnodes <- xml2::xml_find_all(myXML, ".//dimEntity.structureId")
  myResNodes <- xml2::xml_find_all(myXML, ".//dimStructure.resolution")
  myChainNodes <- xml2::xml_find_all(myXML, ".//dimEntity.chainId")
  mySeqNodes <- xml2::xml_find_all(myXML, ".//dimEntity.sequence")

  myPDBData <- data.frame(IDs =
                          paste(tolower(
                                as.character(xml2::xml_contents(myIDnodes))),
                                as.character(xml2::xml_contents(myChainNodes)),
                                sep = "."),
                          Resolution = as.character(xml2::xml_contents(myResNodes)),
                          ChainIDs = as.character(xml2::xml_contents(myChainNodes)),
                          Sequence = as.character(xml2::xml_contents(mySeqNodes)),
                          stringsAsFactors = FALSE)

  return (myPDBData)

}


# [END]





getData <- function() {

  # Load and install required dependencies if required
  if (!requireNamespace("biocLite", quietly=TRUE)) {
    install.packages("biocLite")
  }

  if (!requireNamespace("biomaRt", quietly=TRUE)) {
    biocLite("biomaRt")
  }

  library(biomaRt)

  # Retrieve HGNC symbols from ensembl biomaRt dataset
  sapienData <- useEnsembl("ensembl", "hsapiens_gene_ensembl")
  attributes <- c("hgnc_symbol", "pdb", "sifts_import",
                   "hgnc_trans_name", "ensembl_transcript_id")

  geneList <- getBM(attributes = attributes, mart = sapienData)

  # Remove any sapien genes that do not have a HGNC symbol
  # or PDB entry
  geneList <- geneList[(("" != geneList[,1]) & ("" != geneList[,2])
                        & ("" != geneList[,3])),]

  # Create database
  myDB <- createDB(geneList)

  # Retrieve sequences
  seq <- getSequence(id = "ENST00000288986", type = "ensembl_transcript_id",
                     mart = sapienData, seqType = "peptide")
  show(seq)

}

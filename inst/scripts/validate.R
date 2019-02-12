#' validate.R
#'
#' @description This function provides users with validation on the
#' PDB and HGNC database generated from the BCB420.2019.PDB
#' package. It performs a MSA between the provided Transcript
#' and PDB chain sequences to determine if the two sequenes are
#' identical. It also verifies that for both the
#' Transcript and PDB chain map to the same gene coordinates.
#'
#' @param transcriptHGNC The stable transcript ID of a transcript
#'
#' @param chainID The ID of the PDB chain
#'
#' @param myDB The annotated database loaded in the environment
#'
#' @param showConsensus A boolean value determining if the user
#' would like to see the consensus sequence.
#' Set to TRUE automatically
#'
#' @param showPrint A boolean value determining if the user
#' would like to have a PDF file made of the MSA.
#' Set to TRUE automatically
#'
#' @export
#'

validate <- function(transcriptHGNC, chainID, myDB, gffData,
                     showConsensus = TRUE) {

  ### ============ Validate transcript and PDB chain sequences ============= ###

  # Retrieve the transcript ID to fetch the sequence from biomaRt
  transcriptID <- myDB$Transcripts[myDB$Transcripts$ID == transcriptHGNC,]$stableID

  # Retrieve transcript sequences using Ensembl
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  transcriptSeq <- biomaRt::getSequence(id = transcriptID, type = "ensembl_transcript_id",
                     mart = mart, seqType = "peptide")

  # Retrieve PDB chain sequences through the annotated database
  chainSeq <- myDB$pdbChains[chainID == myDB$pdbChains$ID,]$Sequence


  # Use multiple sequence alignment to validate sequences
  msaSeqs <- Biostrings::AAStringSet(c(transcriptSeq[[1]], chainSeq))
  names(msaSeqs) <- c(as.character(transcriptHGNC), as.character(chainID))
  myMsa <- msa::msa(msaSeqs, method = "ClustalW",
                    type = "protein")


  # Show the consensus sequence and generate printed MSA if
  # the user requests it
  if (showConsensus) {
    print(msaConsensusSequence(myMsa))
  }


  ### ================ Validate gene coordinates =============== ###

  coordinateData <- gffData[!is.na(gffData$transcript_id) & !is.na(gffData$start)
                     & !is.na(gffData$end),]
  coordinateData <- unique(gffData[grepl(transcriptID,
                                  gffData$transcript_id),][,c("transcript_id","start","end")])

  gffStart <- coordinateData$start[1]
  gffEnd <- coordinateData$end[1]

  pdbStart <- myDB$Transcripts[as.character(transcriptHGNC)
                               == myDB$Transcripts$ID,]$start
  pdbEnd <- myDB$Transcripts[as.character(transcriptHGNC)
                             == myDB$Transcripts$ID,]$end

  if (gffStart == pdbStart & gffEnd == pdbEnd) {
    message ("The coordinates are the same for the transcript and the
             PDB chain.")

  }

  else {
    message("The coordinates are not equal for the transcript and
            the PDB chain")
  }
}

# [END]


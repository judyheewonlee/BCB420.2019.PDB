## PDB Dataproject

###  PDB Dataproject for BCB420 2019

-----------------------------------------------

Version: 1.0  
Author: Judy Heewon Lee (heewon.lee@mail.utoronto.ca)  
Versions: 1.0 

-----------------------------------------------


This package describes the workflow and relation between the PDB database from the PDB database and how to map the PDB IDs to HGNC symbols. It takes into account the PDB ID along with the corresponding transcripts and HGNC symbls and HGNC transcripts. The biggest part of this code is representing the relations between PDB chains and how they map to specific HGNC transcripts and symbols.

----------------------------------------------
This package follows the structure and process 
suggested by Hadley Wickham in:

R Packages
http://r-pkgs.had.co.nz/

This package is also based off the rpt package developed by Professor Boris Steipe:

https://github.com/hyginn/rpt

-----------------------------------------------
Some useful keyboard shortcuts for package authoring:

Build and Reload Package:  'Cmd + Shift + B'  
Update Documentation:      'Cmd + Shift + D' or devtools::document()  
Check Package:             'Cmd + Shift + E'  

-----------------------------------------------

Load the package with:  
devtools::install_github("judyheewonlee/BCB420.2019.PDB")

# `BCB420.2019.PDB`

#### (PDB data annotatation of human genes)

&nbsp;

###### [Judy Heewon Lee](http://steipe.biochemistry.utoronto.ca/abc/students/index.php/User:Judy_Lee), BCB420 Project. Bioinformatics specialist and Cell systems biology major.

----

**If any of this information is ambiguous, inaccurate, outdated, or incomplete, please check the [most recent version](https://github.com/judyheewonlee/BCB420.2019.PDB) of the package on GitHub and if the problem has not already been addressed, please [file an issue](https://github.com/judyheewonlee/BCB420.2019.PDB/issues).**

----

## 1 About this package:

This package describes the workflow and relation between the PDB database from the [PDB](https://www.rcsb.org)  and how to map the PDB IDs to HGNC symbols. It takes into account the PDB ID (more specifically the PDB chain) along with the corresponding transcripts and HGNC symbols. The biggest part of this dataset is representing the relations between PDB chains and how they map to specific HGNC transcripts and symbols. Data is retrieved from both the [PDB](https://www.rcsb.org) and [BiomaRt Ensembl databases](http://useast.ensembl.org/index.html) to build the database and validate it.

The package serves dual duty, as an RStudio project, as well as an R package that can be installed. Package checks **pass without errors, warnings, or notes**.

&nbsp;

#### In this project ...

```text
--BCB420.2019.STRING/
|__.gitignore
|__.Rbuildignore
|__BCB420.2019.STRING.Rproj
|__DESCRIPTION
|__dev/
|__rptTwee.R
|__toBrowser.R               # display .md files in your browser
|__inst/
|__extdata/
|__ensp2sym.RData         # ENSP ID to HGNC symbol mapping tool
|__xSetEdges.tsv          # annotated example edges
|__img/
|__[...]                  # image sources for .md document
|__scripts/
|__recoverIDs.R           # utility to use biomaRt for ID mapping
|__LICENSE
|__NAMESPACE
|__R/
|__zzz.R
|__README.md                    # this file

```

&nbsp;

----

## 2 PDB, Ensembl and Gencode Data

The PDB is one of the largest datasources for annotations on biological proteins. The PDB provides access to 3D structure data for large biological molecules (proteins, DNA, and RNA). Data files contained in the [PDB archive](ftp://ftp.wwpdb.org) are free of all copyright restrictions and made fully and freely available for both non-commercial and commercial use.

The [Ensembl](http://useast.ensembl.org/info/genome/genebuild/genome_annotation.html) gene annotation database includes automatic annotation of genome-wide determination of transcripts. The database has several species with annotations such as humans, mice and zebrafish which can be manually curated as determined by the user.
 
[Gencode](https://www.gencodegenes.org/) is another database that identifies and classifies all gene features in the human and mouse genome with high accuracy. They base their information on biological evidence and release these annotations for users to use. The PDB often refers to the Gencode database. For the development of this workflow, we access Gencode in order to retrieve gene coordinates of PDB chains since the PDB uses Gencode to display gene coordinate information in their database. 

It should be noted that both Gencode and Ensembl carry almost identical data, but for the sake of removing ambiguity we will access Gencode and compare it with Ensembl data during verification of the dataset. More information on the correlation between Gencode and Ensembl can be found on their [FAQ](https://www.gencodegenes.org/pages/faq.html) and [publications](https://europepmc.org/articles/PMC4602055).

&nbsp;

#### 2.1 Ensembl Data semantics

Ensembl BiomaRt gene data has several possible annotations for each gene. Here there are several annotations that the PDB dataset will be using in order to map specific PDB chains to HGNC symbols:

1. **HGNC symbol**: HGNC symbol for a gene
2. **Gene Description**: A basic description of a gene
3. **Gene Stable ID**: The gene ID 
4. **Gene Stable ID Version**: The gene ID with the version number
5. **Gene name**: The name of the gene, commonly similar to the HGNC symbol
6. **Transcript Stable ID**: The transcript ID
7. **Transcript Stable ID Version**: The transcript ID with the version number
8. **Transcript Start**: The start coordinate of the transcript in a gene given by bp
9. **Transcript End**: The end coordinate of the transcript in a gene given by bp
10. **PDB-ENSP mappings**: PDB chain mappings to specific transcripts
11. **PDB IDs**: The associated PDB ID to a gene

For each gene, BiomaRt has several annotations.  Each gene has several transcripts, and each transcript has several possible PDB-ENSP mappings. 

**HGNC symbols, Transcript IDs, PDB-ENSP mappings and PDB IDs** are utilized in order to map the workflow of the dataset. **Transcript Start** and **Transcript End** are annotations used during **validation** of the dataset. Other attributes are annotations for the user to retrieve if interested.

#### PDB Data semantics

The PDB contains a very large amount of annotations for each protein structure. For the purpose of creating an annotated database showing the workflow between HGNC symbols and PDB structures, the following data was retrieved from the database:

1. **Resolution**: The resolution of the PDB structure
2. **Sequence**: The sequence of the PDB structure

These annotations are mainly utilized in order to **validate** the workflow.  Resolution is utilized as a way for the database to determine the best representative if PDB domains overlap one another. However, they have been added as annotations to the database for the user to retrieve and use.

#### Gencode Data semantics

Gencode

The Gencode dataset also includes several annotations in the human and mouse genome. In the workflow, we will be specifically retrieving the coordinates for each transcript selected from Ensembl.

1. **Transcript Stable ID Version**: The transcript ID with the version number
2. **Transcript Start**: The start coordinate of the transcript in a gene given by bp
3. **Transcript End**: The end coordinate of the transcript in a gene given by bp

The Gencode annotationa will be used in order to validate the data in the workflow by comparing these annotations with the same annotations provided by Ensembl.

&nbsp;

## 3 Data download and cleanup

To download the source data from BiomaRt ..:

1. To install the source code from BiomaRt follow the link [here](http://useast.ensembl.org/biomart/martview/9a56e6581b5b925d805764efae37fe5b?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.description|hsapiens_gene_ensembl.default.feature_page.transcript_start|hsapiens_gene_ensembl.default.feature_page.transcript_end|hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id_version|hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id|hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id_version|hsapiens_gene_ensembl.default.feature_page.hgnc_symbol|hsapiens_gene_ensembl.default.feature_page.pdb|hsapiens_gene_ensembl.default.feature_page.sifts_import|hsapiens_gene_ensembl.default.feature_page.hgnc_trans_name|hsapiens_gene_ensembl.default.feature_page.external_gene_name&FILTERS=&VISIBLEPANEL=attributepanel).

3. Download the following data file by clicking **Results** and selecting Export all files to **Compressed file (.gz)** and file type **CSV**: (Warning: large).

* `mart_export.txt.gz` (71.2 Mb)   All human genes with required annotations;

4. Uncompress the file and place it into a sister directory of your working directory which is called `data`. (It should be reachable with `file.path("..", "data")`). **Warning:**  `../data/mart_export.txt` is  1.73 GB.

To download the source data from Gencode:

1. Follow the link to the FTP site for Gencode [here](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/).

2. Select and download  `gencode.v27.basic.annotation.gff3.gz`. (28.4 Mb) (Warning: large)

3. Uncompress the file and place it into a sister directory of your working directory which is called `data`. (It should be reachable with `file.path("..", "data")`). **Warning:**  `../data/gencode.v29.basic.annotation.gff3` is 771 Mb.

&nbsp;

## 4 Mapping ENSEMBL IDs to HGNC symbols

STRING network nodes are Ensembl protein IDs. These can usually be mapped to HGNC symbols, but there might be ambiguities e.g. because alternatively spliced proteins might have different ENSP IDs that  map to the same HGNC symbol, or HGNC symbols have changed (they are frequently updated). To provide the best possible interpretation, we need to build a map of ENSP IDs to HGNC symbols. This requires care, because it is not guaranteed that all ENSP IDs can be mapped uniquely.** However, the usability of the dataset for annotation depends on the quality of this mapping.**

&nbsp;

#### Preparations: packages, functions, files

To begin, we need to make sure the required packages are installed:

**`xml2`** provides functions to read XML formatted data. This package will be used in order to retrieve sequence and resolution annotations of each transcript from the PDB RESTful API services.
&nbsp;

```R
if (!requireNamespace("xml2", quietly = TRUE)) {
install.packages("xml2")
}
```

**`rtracklayer`** provides functions to read `.gff3` files. This package will be used in order to read Gencode data from the file `gencode.v29.basic.annotation.gff3`.It is a Bioconductor package, and as such it needs to be loaded via the **`BiocManager`**,
&nbsp;

```R
if (! requireNamespace("BiocManager", quietly = TRUE)) {
install.packages("BiocManager")
}
if (!requireNamespace("rtracklayer", quietly = TRUE)) {
BiocManager::install("rtracklayer")
}

}
```

**`biomaRt`** biomaRt is a Bioconductor package that implements the RESTful API of biomart,
the annotation framwork for model organism genomes at the EBI. It is a Bioconductor package, and as such it needs to be loaded via the **`BiocManager`**,
&nbsp;

```R
if (! requireNamespace("BiocManager", quietly = TRUE)) {
install.packages("BiocManager")
}
if (! requireNamespace("biomaRt", quietly = TRUE)) {
BiocManager::install("biomaRt")
}
```

**`msa`** msa is a Bioconductor package that provides utilities to create multiple sequence alignments using provided sequences. It needs to be loaded via **`BiocManager`** like the `biomaRt` package. This package is used in order to validate the protein chain sequences provided by the PDB and the sequences of each transcript from Ensembl.
&nbsp;

```R
if (! requireNamespace("BiocManager", quietly = TRUE)) {
install.packages("BiocManager")
}
if (!requireNamespace("msa", quietly = TRUE)) {
BiocManager::install("msa", version = "3.8")
}
```

**`igraph`** is THE go-to package for everything graph related. We use it here to
compute some statistics on the STRING- and example graphs and plot.
&nbsp;

```R
if (! requireNamespace("igraph")) {
install.packages("igraph")
}
```

&nbsp;

Next we source a utility functions that can be run in order to map and create the PDB-HGNC database.
(Note: The functions have been separated into several different R files for organization)

&nbsp;

```R
source("inst/scripts/PDBdataset.R")
source("inst/scripts/getData.R")
source("inst/scripts/addHGNC.R")
source("inst/scripts/addTranscripts.R")
source("inst/scripts/addpdbChain.R")
source("inst/scripts/addPDB.R")
source("inst/scripts/addPDBData.R")
source("inst/scripts/fetchPDBXML.R")
source("inst/scripts/getGencode.R")
source("inst/scripts/validate.R")
```

&nbsp;

The user may simply use the function call in order to get the database:

```R
myDB <- PDBdataset()
```
Nonetheless, there will be a walkthrough of each function below showing the workflow of building the dataset.

&nbsp;

#### 4.1 Step one: which IDs do we have to map?

&nbsp;

```R
# Read the biomaRt Ensembl data from the data directory and create a dataframe

martFile <- file.path("../data", "mart_export.txt")
myDF <- read.csv(martFile)

# what does each row look like?
myDF[1,]
# Gene.stable.ID Transcript.stable.ID Gene.stable.ID.version Transcript.stable.ID.version  # Transcript.start..bp.
# 1 ENSG00000198888      ENST00000361390      ENSG00000198888.2                 
# ENST00000361390.2                  3307
# Transcript.end..bp. HGNC.symbol PDB.ENSP.mappings Protein.stable.ID 
# Protein.stable.ID.version
# 1                4262      MT-ND1            5xtc.s   ENSP00000354687        
# ENSP00000354687.2
# Gene.description
# 1 mitochondrially encoded NADH:ubiquinone oxidoreductase core subunit 1 [Source:HGNC 
# Symbol;Acc:HGNC:7455]
# Gene.start..bp. Gene.end..bp. Gene.name   HGNC.ID PDB.ID
# 1            3307          4262    MT-ND1 HGNC:7455   5XTD  

# Each gene has several annotations which are titled in the dataframe

# We want to remove any sapien genes that do not have a HGNC symbol
# or PDB entry or PDB chain mappings
myDF <- myDF[(("" != myDF$HGNC.symbol) & ("" != myDF$Transcript.stable.ID)
& ("" != myDF$PDB.ENSP.mappings)),]

# how many unique HGNC IDs do we have to map?
uHGNC <- unique(myDF$HGNC.symbol)  # 799 HGNC IDs need to be mapped

# how many unique transcript IDs do we have to map?
uTranscript <- unique(myDF$Transcript.stable.ID) # 1158 Transcript IDs need to be mapped

# how many unique PDB chain IDs do we have to map?
uPDBchains <- unique(myDF$PDB.ENSP.mappings) # 12006 

# how many unique PDB IDs do we have to map?
uPDB <- unique(myDF$PDB.ID) # 5722

```

&nbsp;

#### 4.2  Step two: mapping via biomaRt

We first fetch all unique HGNC symbols and create a dataframe entry for the database containing the HGNC symbols with the corresponding annotations: **HGNC symbol**, **Gene Description**, **Gene Stable ID**, **Gene Stable ID Version**, and **Gene name**. Every HGNC symbol can map to a large number of transcripts, as genes are constructed from many transcripts. Thus every HGNC symbol will have multiple transcripts but each transcript will only map to a single HGNC symbol. 

Once we generate the dataframe containing the HGNC symbols, we map the **transcript IDs** to several HGNC symbols. In this data frame, we will include **Transcript Stable ID**, **Transcript Stable ID Version**, **Transcript HGNC ID**, **Transcript Start**, and **Transcript End** as annotations. Every transcript can have multiple PDB chain IDs and every PDB chain can also have multiple transcripts. We must take this into account when building the database.

After we generate the transcript table, we map the PDB chains to each transcript ID. In the data frame containing the PDB chains, we will add the annotations: **PDB-ENSP mapping (chain ID)**, **Transcript IDs**, **Sequence** and **Resolution**. The **Sequence** and **Resolution** will be retrieved through the PDB API services discussed later. 

Lastly, we map each PDB-chain to a PDB ID to the general PDB structure. Oddly, some entries in the BiomaRt dataset have different corresponding chain ID's to the PDB ID. Most cases, the PDB ID corresponds to the PDB chain. For example: 

`5xtc.s` is a PDB chain ID where the suffix after the period is the chain identifier for the PDB protein ID. Sometimes these PDB chains may map to diferent PDB IDs. 

We will dicuss the problems and the solutions in retrieving this information and mapping it correctly in order to create a proper workflow.

&nbsp;

###### 4.2.1  Constructing the Dataset

&nbsp;

```R

# Since we already have read the data file for BiomaRt and 
# assigned it as myDF:
# nrow(myDF)  # 1368671 HGNC entries are present in the dataset ( with the entries without # PDB-ENSP mappings, HGNC symbols or transcript IDs removed )

# Note that the number of entries does not reflect the number of HGNC symbols. 
# There are many different combinations of transcripts and PDB chains that produce the high # level of entries in this dataset

```
&nbsp;



&nbsp;

**(1)** There might be more than one value returned. The ID appears more than
once in `tmp$ensembl_peptide_id`, with different mapped symbols.

```R
sum(duplicated(tmp$ensembl_peptide_id))  # Indeed: three duplicates!
```

&nbsp;

**(2)** There might be nothing returned for one ENSP ID. We have the ID in `uENSP`, but it does not appear in `tmp$ensembl_peptide_id`:

```R

sum(! (uENSP) %in% tmp$ensembl_peptide_id)  # 248
```
&nbsp;

**(3)** There might be no value returned: `NA`, or `""`. The ID appears in `tmp$ensembl_peptide_id`, but there is no symbol in `tmp$hgnc_symbol`.

```R
sum(is.na(ensp2sym$sym))  # 0
sum(ensp2sym$sym == "")   # 199 - note: empty strings for absent symbols.
```

&nbsp;

Let's fix the "duplicates" problem first. We can't have duplicates: if we encounter an ENSP ID, we need exactly one symbol assigned to it. What are these genes?

&nbsp;

```R

dupEnsp <- tmp$ensembl_peptide_id[duplicated(tmp$ensembl_peptide_id)]
tmp[tmp$ensembl_peptide_id %in% dupEnsp, ]

#                  ensp      sym
# 8668  ENSP00000344961  PLEKHG7
# 8669  ENSP00000344961 C12orf74
# 14086 ENSP00000380933  PLEKHG7
# 14087 ENSP00000380933 C12orf74
# 18419 ENSP00000480558   CCL3L3
# 18420 ENSP00000480558   CCL3L1

# ENSP00000380933 and ENSP00000344961 should both map to PLEKHG7
# CCL3L3 and CCL3L3 both have UniProt ID P16619, we map ENSP00000480558
# (arbitrarily) to CCL3L1

# validate target rows
tmp[tmp$hgnc_symbol %in% c("C12orf74", "CCL3L3"), ]

# remove target rows
tmp <- tmp[ ! (tmp$hgnc_symbol %in% c("C12orf74", "CCL3L3")), ]

# check result
any(duplicated(tmp$ensembl_peptide_id))   # now FALSE
```

&nbsp;

After this preliminary cleanup, defining the mapping tool is simple:

&nbsp;

```R
ensp2sym <- tmp$hgnc_symbol
names(ensp2sym) <- tmp$ensembl_peptide_id

head(ensp2sym)
# ENSP00000216487 ENSP00000075120 ENSP00000209884  
#        "RIN3"        "SLC2A3"        "KLHL20"   
#        
# ENSP00000046087 ENSP00000205214 ENSP00000167106
#      "ZPBP"         "AASDH"         "VASH1"

```

&nbsp;

###### 4.2.2  Cleanup and validation of `ensp2sym`

There are two types of IDs we need to process further: (1), those that were not returned at all from biomaRt, (2) those for which only an empty string was returned.

First, we add the symbols that were not returned by biomaRt to the map. They are present in uENSP, but not in ensp2sym$ensp:

&nbsp;

```R
sel <- ! (uENSP %in% names(ensp2sym))
x <- rep(NA, sum( sel))
names(x) <- uENSP[ sel ]

# confirm uniqueness
any(duplicated(c(names(x), names(ensp2sym))))  # FALSE

# concatenate the two vectors
ensp2sym <- c(ensp2sym, x)

# confirm
all(uENSP %in% names(ensp2sym))  # TRUE
```

&nbsp;

Next, we set the symbols for which only an empty string was returned to `NA`:

&nbsp;

```R
sel <- which(ensp2sym == "") # 199 elements
ensp2sym[head(sel)] # before ...
ensp2sym[sel] <- NA
ensp2sym[head(sel)] # ... after

# Do we still have all ENSP IDs accounted for?
all( uENSP %in% names(ensp2sym))  # TRUE

```

&nbsp;

###### 4.2.3  Additional symbols

A function for using biomaRt for more detailed mapping is in the file `inst/scripts/recoverIds.R`. We have loaded it previously, and use it on all elements of `ensp2sym` that are `NA`.

&nbsp;

```R

# How many NAs are there in "ensp2sym" column?
sum(is.na(ensp2sym))   # 447

# subset the ENSP IDs
unmappedENSP <- names(ensp2sym)[is.na(ensp2sym)]

# use our function recoverIDs() to try and map the unmapped ensp IDs
# to symboils via other cross-references
recoveredENSP <- recoverIDs(unmappedENSP)

# how many did we find
nrow(recoveredENSP)  # 11. Not much, but it's honest work.

# add the recovered symbols to ensp2sym
ensp2sym[recoveredENSP$ensp] <- recoveredENSP$sym

# validate:
sum(is.na(ensp2sym))  # 436 - 11 less than 447

```

&nbsp;

#### 4.4  Step four: outdated symbols

We now have each unique ENSP IDs represented once in our mapping table. But are these the correct symbols? Or did biomaRt return obsolete names for some? We need to compare the symbols to our reference data and try to fix any problems. Symbols that do not appear in the reference table will also be set to NA.

&nbsp;

```R
# are all symbols present in the reference table?
sel <- ( ! (ensp2sym %in% HGNC$sym)) & ( ! (is.na(ensp2sym)))
length(        ensp2sym[ sel ] )  # 137 unknown
length( unique(ensp2sym[ sel ]))  # they are all unique

# put these symbols in a new dataframe
unkSym <- data.frame(unk = ensp2sym[ sel ],
new = NA,
stringsAsFactors = FALSE)

# Inspect:
# several of these are formatted like "TNFSF12-TNFSF13" or "TMED7-TICAM2".
# This looks like biomaRt concatenated symbol names.
grep("TNFSF12", HGNC$sym) # 23984: TNFSF12
grep("TNFSF13", HGNC$sym) # 23985 23986: TNFSF13 and TNFSF13B
grep("TMED7",   HGNC$sym) # 23630: TMED7
grep("TICAM2",  HGNC$sym) # 23494: TICAM2

# It's not clear why this happened. We will take a conservative approach
# and not make assumptions which of the two symbols is the correct one,
# i.e. we will leave these symbols as NA


# grep() for the presence of the symbols in either HGNC$prev or
# HGNC$synonym. If either is found, that symbol replaces NA in unkSym$new
for (i in seq_len(nrow(unkSym))) {
iPrev <- grep(unkSym$unk[i], HGNC$prev)[1] # take No. 1 if there are several
if (length(iPrev) == 1) {
unkSym$new[i] <- HGNC$sym[iPrev]
} else {
iSynonym <- which(grep(unkSym$unk[i], HGNC$synonym))[1]
if (length(iSynonym) == 1) {
unkSym$new[i] <- HGNC$sym[iSynonym]
}
}
}

# How many did we find?
sum(! is.na(unkSym$new))  # 32

# We add the contents of unkSym$new back into ensp2sym. This way, the
# newly mapped symbols are updated, and the old symbols that did not
# map are set to NA.

ensp2sym[rownames(unkSym)] <- unkSym$new


```

#### 4.5 Final validation

Validation and statistics of our mapping tool:

```R

# do we now have all ENSP IDs mapped?
all(uENSP %in% names(ensp2sym))  # TRUE

# how many symbols did we find?
sum(! is.na(ensp2sym))  # 18845

# (in %)
sum(! is.na(ensp2sym)) * 100 / length(ensp2sym)  # 96.0 %

# are all symbols current in our reference table?
all(ensp2sym[! is.na(ensp2sym)] %in% HGNC$sym)  # TRUE

# Done.
# This concludes construction of our mapping tool.
# Save the map:

save(ensp2sym, file = file.path("inst", "extdata", "ensp2sym.RData"))

# From an RStudio project, the file can be loaded with
load(file = file.path("inst", "extdata", "ensp2sym.RData"))


```

&nbsp;

# 5 Annotating gene sets with STRING Data

Given our mapping tool, we can now annotate gene sets with STRING data. As a first example, we analyze the entire STRING graph. Next, we use high-confidence edges to analyze the network of our example gene set.


&nbsp;

```R

# Read the interaction graph data: this is a weighted graph defined as an
# edge list with gene a, gene b, confidence score (0, 999).

tmp <- readr::read_delim(file.path("../data", "9606.protein.links.v11.0.txt"),
delim = " ",
skip = 1,
col_names = c("a", "b", "score"))  # 11,759,454 rows

# do they all have the right tax id?
all(grepl("^9606\\.", tmp$a))  # TRUE
all(grepl("^9606\\.", tmp$b))  # TRUE
# remove "9606." prefix
tmp$a <- gsub("^9606\\.", "", tmp$a)
tmp$b <- gsub("^9606\\.", "", tmp$b)

# how are the scores distributed?

minScore <- 0
maxScore <- 1000
# we define breaks to lie just below the next full number
hist(tmp$score[(tmp$score >= minScore) & (tmp$score <= maxScore)],
xlim = c(minScore, maxScore),
breaks = c((seq(minScore, (maxScore - 25), by = 25) - 0.1), maxScore),
main = "STRING edge scores",
col = colorRampPalette(c("#FFFFFF","#8888A6","#FF6655"), bias = 2)(40),
xlab = "scores: (p * 1,000)",
ylab = "p",
xaxt = "n")
axis(1, at = seq(minScore, maxScore, by = 100))
abline(v = 900, lwd = 0.5)

```

![](./inst/img/score_hist_1.svg?sanitize=true "STRING score distribution")

We know that "channel 7 - databases" interactions are arbitrarily scored as _p_ = 0.9. This is clearly reflected in the scores distribution.

```R
# Zoom in

minScore <- 860
maxScore <- 1000
hist(tmp$score[(tmp$score >= minScore) & (tmp$score <= maxScore)],
xlim = c(minScore, maxScore),
breaks = c((seq(minScore, (maxScore - 4), by = 4) - 0.1), maxScore),
main = "STRING edge scores",
col = colorRampPalette(c("#FFFFFF","#8888A6","#FF6655"), bias = 1.2)(35),
xlab = "scores: (p * 1,000)",
ylab = "p",
xaxt = "n")
axis(1, at = seq(minScore, maxScore, by = 10))
abline(v = 900, lwd = 0.5)

```

![](./inst/img/score_hist_2.svg?sanitize=true "STRING score distribution (detail)")


```R

# Focus on the cutoff of scores at p == 0.9
sum(tmp$score >= 880 & tmp$score < 890) # 5,706
sum(tmp$score >= 890 & tmp$score < 900) # 5,666
sum(tmp$score >= 900 & tmp$score < 910) # 315,010
sum(tmp$score >= 910 & tmp$score < 920) # 83,756

# We shall restrict our dataset to high-confidence edges with p >= 0.9

tmp <- tmp[tmp$score >= 900, ]  # 648,304 rows of high-confidence edges

```

&nbsp;

Are these edges duplicated? I.e. are there (a, b) and (b, a) edges in the dataset? The common way to test for that is to created a composite string of the two elements, sorted. Thus if we have an edge betwween `"this"` and `"that"`, and an edge between `"that"` and `"this"`, these edges both get mapped to a key `"that:this"` - and the duplication is easy to recognize.

&nbsp;

```R

sPaste <- function(x, collapse = ":") {
return(paste(sort(x), collapse = collapse))
}
tmp$key <- apply(tmp[ , c("a", "b")], 1, sPaste)

length(tmp$key) # 648,304
length(unique(tmp$key)) # 324,152  ... one half of the edges are duplicates!

# We can remove those edges. And the keys.
tmp <- tmp[( ! duplicated(tmp$key)), c("a", "b", "score") ]
```

&nbsp;

Finally we map the ENSP IDs to HGNC symbols. Using our tool, this is a simple assignment:

&nbsp;

```R

tmp$a <- ensp2sym[tmp$a]
tmp$b <- ensp2sym[tmp$b]

# Validate:
# how many rows could not be mapped
any(grepl("ENSP", tmp$a))  # Nope
any(grepl("ENSP", tmp$b))  # None left here either
sum(is.na(tmp$a)) # 705
sum(is.na(tmp$b)) # 3501

# we remove edges in which either one or the other node is NA to
# create our final data:
STRINGedges <- tmp[( ! is.na(tmp$a)) & ( ! is.na(tmp$b)), ] # 319,997 edges

# Done.
# Save result
save(STRINGedges, file = file.path("..", "data", "STRINGedges.RData"))
# That's only 1.4 MB actually.

```

&nbsp;

#### 5 Network statistics

Simple characterization of network statistics:

&nbsp;

```R

# number of nodes
(N <- length(unique(c(STRINGedges$a, STRINGedges$b))))  # 12,196 genes

# coverage of human protein genes
N * 100 / sum(HGNC$type == "protein")  # 63.4 %

# number of edges
nrow(STRINGedges)   # 319,997

# any self-edges?
any(STRINGedges$a == STRINGedges$b) # yes
which(STRINGedges$a == STRINGedges$b)
STRINGedges[which(STRINGedges$a == STRINGedges$b), ]
#        a     b   score     # just one
#  1 ZBED6 ZBED6     940


# average number of interactions
nrow(STRINGedges) / N  # 26.2  ... that seems a lot - how is this distributed?

# degree distribution
deg <- table(c(STRINGedges$a, STRINGedges$b))
summary(as.numeric(deg))

hist(deg, breaks=50,
xlim = c(0, 1400),
col = "#3fafb388",
main = "STRING nodes degree distribution",
xlab = "degree (undirected graph)",
ylab = "Counts")
rug(deg, col = "#EE5544")

```

![](./inst/img/STRING_degrees_1.svg?sanitize=true "STRING network degree distribution")


## 6 Biological validation: network properties

For more detailed validation, we need to look at network properties 

&nbsp;

```R

sG <- igraph::graph_from_edgelist(matrix(c(STRINGedges$a,
STRINGedges$b),
ncol = 2,
byrow = FALSE),
directed = FALSE)

# degree distribution
dg <- igraph::degree(sG)

# is this a scale-free distribution? Plot log(rank) vs. log(frequency)
freqRank <- table(dg)
x <- log10(as.numeric(names(freqRank)) + 1)
y <- log10(as.numeric(freqRank))
plot(x, y,
type = "b",
pch = 21, bg = "#A5F5CC",
xlab = "log(Rank)", ylab = "log(frequency)",
main = "Zipf's law governing the STRING network")

# Regression line
ab <- lm(y ~ x)
abline(ab, col = "#FF000077", lwd = 0.7)

```

![](./inst/img/STRING_Zipf_plot_1.svg?sanitize=true "STRING score distribution (detail)")


```R
# What are the ten highest degree nodes?
x <- sort(dg, decreasing = TRUE)[1:10]
cat(sprintf("\t%d:\t%s\t(%s)\n", x, names(x), HGNC[names(x), "name"]))
# 1343:    RPS27A    (ribosomal protein S27a)
# 1339:    UBA52    (ubiquitin A-52 residue ribosomal protein fusion product 1)
# 1128:    UBC    (ubiquitin C)
# 1124:    UBB    (ubiquitin B)
# 918:    GNB1    (G protein subunit beta 1)
# 894:    GNGT1    (G protein subunit gamma transducin 1)
# 562:    APP    (amyloid beta precursor protein)
# 550:    CDC5L    (cell division cycle 5 like)
# 530:    GNG2    (G protein subunit gamma 2)
# 526:    RBX1    (ring-box 1)


```

&nbsp;

## 7 Annotation of the example gene set

To conclude, we annotate the example gene set, validate the annotation, and store the data in an edge-list format.

&nbsp;

```R

# The specification of the sample set is copy-paste from the 
# BCB420 resources project.

xSet <- c("AMBRA1", "ATG14", "ATP2A1", "ATP2A2", "ATP2A3", "BECN1", "BECN2",
"BIRC6", "BLOC1S1", "BLOC1S2", "BORCS5", "BORCS6", "BORCS7",
"BORCS8", "CACNA1A", "CALCOCO2", "CTTN", "DCTN1", "EPG5", "GABARAP",
"GABARAPL1", "GABARAPL2", "HDAC6", "HSPB8", "INPP5E", "IRGM",
"KXD1", "LAMP1", "LAMP2", "LAMP3", "LAMP5", "MAP1LC3A", "MAP1LC3B",
"MAP1LC3C", "MGRN1", "MYO1C", "MYO6", "NAPA", "NSF", "OPTN",
"OSBPL1A", "PI4K2A", "PIK3C3", "PLEKHM1", "PSEN1", "RAB20", "RAB21",
"RAB29", "RAB34", "RAB39A", "RAB7A", "RAB7B", "RPTOR", "RUBCN",
"RUBCNL", "SNAP29", "SNAP47", "SNAPIN", "SPG11", "STX17", "STX6",
"SYT7", "TARDBP", "TFEB", "TGM2", "TIFA", "TMEM175", "TOM1",
"TPCN1", "TPCN2", "TPPP", "TXNIP", "UVRAG", "VAMP3", "VAMP7",
"VAMP8", "VAPA", "VPS11", "VPS16", "VPS18", "VPS33A", "VPS39",
"VPS41", "VTI1B", "YKT6")

# which example genes are not among the known nodes?
x <- which( ! (xSet %in% c(STRINGedges$a, STRINGedges$b)))
cat(sprintf("\t%s\t(%s)\n", HGNC[xSet[x], "sym"], HGNC[xSet[x], "name"]))

# BECN2    (beclin 2)
# EPG5    (ectopic P-granules autophagy protein 5 homolog)
# LAMP3    (lysosomal associated membrane protein 3)
# LAMP5    (lysosomal associated membrane protein family member 5)
# PLEKHM1    (pleckstrin homology and RUN domain containing M1)
# RUBCNL    (rubicon like autophagy enhancer)
# TIFA    (TRAF interacting protein with forkhead associated domain)
# TMEM175    (transmembrane protein 175)
# TPCN1    (two pore segment channel 1)
# TPCN2    (two pore segment channel 2)

# That make sense - generally fewer interactions have been recorded for
# membrane proteins.


# For our annotation, we select edges for which both nodes are part of the
# example set:
sel <- (STRINGedges$a %in% xSet) & (STRINGedges$b %in% xSet)
xSetEdges <- STRINGedges[sel, c("a", "b")]
# Statistics:
nrow(xSetEdges)   # 206

# Save the annotated set

writeLines(c("a\tb",
sprintf("%s\t%s", xSetEdges$a, xSetEdges$b)),
con = "xSetEdges.tsv")

# The data set can be read back in again (in an RStudio session) with
myXset <- read.delim(file.path("inst", "extdata", "xSetEdges.tsv"),
stringsAsFactors = FALSE)

# From an installed package, the command would be:
myXset <- read.delim(system.file("extdata",
"xSetEdges.tsv",
package = "BCB420.2019.STRING"),
stringsAsFactors = FALSE)


# confirm
nrow(myXset) # 206
colnames(myXset) == c("a", "b") # TRUE TRUE

```

&nbsp;

#### 7.1 Biological validation: network properties

Explore some network properties of the exmple gene set.

&nbsp;

```R

# A graph ...
sXG <- igraph::graph_from_edgelist(matrix(c(xSetEdges$a,
xSetEdges$b),
ncol = 2,
byrow = FALSE),
directed = FALSE)

# degree distribution
dg <- igraph::degree(sXG)
hist(dg, col="#A5CCF5",
main = "Node degrees of example gene network",
xlab = "Degree", ylab = "Counts")

# scale free? log(rank) vs. log(frequency)
freqRank <- table(dg)
x <- log10(as.numeric(names(freqRank)) + 1)
y <- log10(as.numeric(freqRank))
plot(x, y,
type = "b",
pch = 21, bg = "#A5CCF5",
xlab = "log(Rank)", ylab = "log(frequency)",
main = "Zipf's law governing the example gene network")

# Regression line
ab <- lm(y ~ x)
abline(ab, col = "#FF000077", lwd = 0.7)

```

![](./inst/img/xGenes_Zipf_plot_1.svg?sanitize=true "xGenes degree distribution (log(#)/log(f))")


```R

# What are the ten highest degree nodes?
x <- sort(dg, decreasing = TRUE)[1:10]
cat(sprintf("\t%d:\t%s\t(%s)\n", x, names(x), HGNC[names(x), "name"]))

# 15:    VAMP8    (vesicle associated membrane protein 8)
# 15:    RAB7A    (RAB7A, member RAS oncogene family)
# 12:    PIK3C3    (phosphatidylinositol 3-kinase catalytic subunit type 3)
# 12:    GABARAP    (GABA type A receptor-associated protein)
# 12:    SNAP29    (synaptosome associated protein 29)
# 12:    STX17    (syntaxin 17)
# 11:    GABARAPL2    (GABA type A receptor associated protein like 2)
# 11:    BECN1    (beclin 1)
# 11:    GABARAPL1    (GABA type A receptor associated protein like 1)
# 10:    UVRAG    (UV radiation resistance associated)


# Plot the network
oPar <- par(mar= rep(0,4)) # Turn margins off
set.seed(112358)
plot(sXG,
layout = igraph::layout_with_fr(sXG),
vertex.color=heat.colors(max(igraph::degree(sXG)+1))[igraph::degree(sXG)+1],
vertex.size = 1.5 + (1.2 * igraph::degree(sXG)),
vertex.label.cex = 0.2 + (0.025 * igraph::degree(sXG)),
edge.width = 2,
vertex.label = igraph::V(sXG)$name,
vertex.label.family = "sans",
vertex.label.cex = 0.9)
set.seed(NULL)
par(oPar)

# we see several cliques (or near-cliques), possibly indicative of
# physical complexes.

```

![](./inst/img/xGenes_Network_1.svg?sanitize=true "xGenes functional interaction network")


&nbsp;

## 8 References

&nbsp;

Example code for biomaRt was taken taken from `BIN-PPI-Analysis.R` and example code for work with igraph was taken from `FND-MAT-Graphs_and_networks.R`, both in the [ABC-Units project](https://github.com/hyginn/ABC-units) (Steipe, 2016-1019). A preliminary version of a STRING import script was written as [starter code for the 2018 BCB BioHacks Hackathon](https://github.com/hyginn/ABC-units) at the UNiversity of Toronto (Steipe, 2018) - this script draws on the former.

&nbsp;

* Szklarczyk, D., Gable, A. L., Lyon, D., Junge, A., Wyder, S., Huerta-Cepas, J., Simonovic, M., Doncheva, N. T., Morris, J. H., Bork, P., Jensen, L. J., & von Mering, C. (2019). STRING v11: protein-protein association networks with increased coverage, supporting functional discovery in genome-wide experimental datasets. [_Nucleic acids research_, D1, D607-D613](https://academic.oup.com/nar/article/47/D1/D607/5198476).

* Huang, J. K., Carlin, D. E., Yu, M. K., Zhang, W., Kreisberg, J. F., Tamayo, P., & Ideker, T. (2018). Systematic Evaluation of Molecular Networks for Discovery of Disease Genes. _Cell systems_, 4, 484-495.e5.

&nbsp;

## 9 Acknowledgements

Thanks to Simon Kågedal's very useful [PubMed to APA reference tool](http://helgo.net/simon/pubmed/).

User `Potherca` [posted on Stack](https://stackoverflow.com/questions/13808020/include-an-svg-hosted-on-github-in-markdown) how to use the parameter `?sanitize=true` to display `.svg` images in github markdown.

&nbsp;

&nbsp;

<!-- [END] -->



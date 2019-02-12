# `BCB420.2019.PDB`

###  PDB Dataproject for BCB420 2019 (PDB data annotatation of human genes)

-----------------------------------------------

Version: 1.0  
Author: Judy Heewon Lee (heewon.lee@mail.utoronto.ca)  
Versions: 1.0 

-----------------------------------------------

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
--BCB420.2019.PDB/
|__.gitignore
|__.Rbuildignore
|__BCB420.2019.PDB.Rproj
|__DESCRIPTION              
|__inst/
    |__extdata/               # ENSP ID to HGNC symbol mapping tool
        |__xSetPDB.tsv          # annotated example edges
    |__img/
        |__[...]                  # image sources for .md document
    |__scripts/
        |__PDBdataset.R           # utilities for ID mapping, more details below
        |__getData.R
        |__addHGNC.R
        |__addTranscripts.R
        |__addpdbChain.R
        |__addPDB.R
        |__addPDBData.R
        |__fetchPDBXML.R
        |__validate.R  
        |__readXML.R
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

Here is a relational model showing the data and how it will be layed out in the dataset.

![](./inst/img/HGNC-PDB1.0.svg?sanitize=true "PDB-HGNC dataset annotation Relational Model")

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

**`data.table`** provides functions to manipulate and create `data.table` objects in R. This package is used in order to use the function `unique()` on data.table objects over data.frame objects to produce a significantly faster runtime.

&nbsp;
```R
if (!requireNamespace("data.table", quietly = TRUE)) {
    install.packages("data.table")
}
```

&nbsp;

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
source("inst/scripts/readXML.R")
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
martDFall <- read.csv(martFile)

# what does each row look like?
martDFall[1,]
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
martDF <- martDF[(("" != martDF$HGNC.symbol) & ("" != martDF$Transcript.stable.ID)
& ("" != martDF$PDB.ENSP.mappings)),]

# how many unique HGNC IDs do we have to map?
uHGNC <- unique(martDF$HGNC.symbol)  # 2554 HGNC IDs need to be mapped

# how many unique transcript IDs do we have to map?
uTranscript <- unique(martDF$Transcript.stable.ID) # 3783 Transcript IDs need to be mapped

# how many unique PDB chain IDs do we have to map?
uPDBchains <- unique(martDF$PDB.ENSP.mappings) # 37329 

# how many unique PDB IDs do we have to map?
uPDB <- unique(martDF$PDB.ID) # 16911

```

&nbsp;

#### 4.2  Step two: mapping and annotating via biomaRt

We first fetch all unique HGNC symbols and create a dataframe entry for the database containing the HGNC symbols with the corresponding annotations: **HGNC symbol**, **Gene Description**, **Gene Stable ID**, **Gene Stable ID Version**, and **Gene name**. Every HGNC symbol can map to a large number of transcripts, as genes are constructed from many transcripts. Thus every HGNC symbol will have multiple transcripts but each transcript will only map to a single HGNC symbol. 

Once we generate the dataframe containing the HGNC symbols, we map the **transcript IDs** to several HGNC symbols. In this data frame, we will include **Transcript HGNC ID**, **Transcript Stable ID Version**, **Transcript Stable ID**, **Transcript Start**, and **Transcript End** as annotations. Every transcript can have multiple PDB chain IDs and every PDB chain can also have multiple transcripts. We must take this into account when building the database.

After we generate the transcript table, we map the PDB chains to each transcript ID. In the data frame containing the PDB chains, we will add the annotations: **PDB-ENSP mapping (chain ID)**, **Transcript HGNC IDs**, **Sequences** and **Resolutions**. The **Sequence** and **Resolution** will be retrieved through the PDB API services discussed later. 

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

# To begin construction of the dataset, we will begin with creating the database. Below is # the general function to create the database:

PDBdataset <- function() {

# Fetch Data from ensembl and retrieve a dataframe
martDF <- getData()

message("Building the database...")
# Build the database
myDB <- list()
myDB <- addHGNC(myDB, martDF)
myDB <- addTranscripts(myDB, martDF)

message("Adding PDB chains...")
myDB <- addpdbChain(myDB, martDF)

message("Adding PDB IDS ...")
myDB <- addPDB(myDB, martDF)
myDB <- addPDBdata(myDB)

message("Database successfully generated!")

return(myDB)

}

# [END]

```

&nbsp;

Each of the functions in `PDBdataset` generates all the necessary dataframes and returns the final dataset. Each funtion will be discussed below. Note: The `getData()` function simply runs the code to read the `mart_export.txt` file discussed in **section 4.1**. It returns the mart data as a full dataframe.

&nbsp;

```R

getData <- function() {

message("Reading Ensembl data...")
martFile <- file.path("../data", "mart_export.txt")

martDF <- read.csv(martFile, stringsAsFactors = FALSE)

# Remove any sapien genes that do not have a HGNC symbol
# or PDB entry
martDF <- martDF[(("" != martDF$HGNC.symbol) & ("" != martDF$HGNC.transcript.name.ID)
& ("" != martDF$PDB.ENSP.mappings)),]

return(martDF)

}

```

&nbsp;

###### 4.2.2  Constructing the HGNC dataframe

&nbsp;

Each HGNC symbol will be placed as a unique entry into `myDB` and will contain annotations for each symbol. 

```R
addHGNC <- function(myDB, martDF) {

    # Here we want to remove any duplicated HGNC symbols since we want to map 
    # the unique HGNC symbols. 
    geneRows <- martDF[!duplicated(martDF$HGNC.symbol),]

    # We generate a dataframe containing the HGNC symbol (which will be used to map 
    # the symbols with the transcript HGNC symbols) and the corresponding annotations.
    myDB$HGNC <- data.frame(ID = geneRows$HGNC.symbol,
    description = geneRows$Gene.description,
    stableID = geneRows$Gene.stable.ID,
    version = geneRows$Gene.stable.ID.version,
    name = geneRows$Gene.name,
    stringsAsFactors = FALSE)

    # Finally we return the modified database.
    return (myDB)
}

```

###### 4.2.3  Constructing the Transcript HGNC ID dataframe

&nbsp;

For mapping transcripts, we come across the problem with duplications. Every HGNC symbol can have more than one transcript ID, but every transcript maps to a unique HGNC symbol. How do we prove this? We can retrieve each row with a unique transcript in the data frame retrieved from the BiomaRt data and get each HGNC symbol corresponding to each transcript. Then we can see if all HGNC symbols are contained in the set of unique HGNC symbols for the whole dataset from BiomaRt.

```R

# Retrieve the unique rows of transcript IDs
uTranscriptTable <- martDF[unique(martDF$Transcript.stable.ID),]

# Test to see if all unique transcript HGNC symbols are in uHGNC. uHGNC defined in section # 4.1 of the read me
all(uTranscriptTable$HGNC.symbol %in% uHGNC)

# [1] TRUE

```
&nbsp;

Since this is true, we can call the `addTranscript()` function which adds the dataframe for each unique transcript ID from the `martDF`

```R
addTranscripts <- function(myDB, martDF) {

    # Retrieve rows with unduplicated transcript IDs
    geneRows <- martDF[!duplicated(martDF$Transcript.stable.ID),]
    
    # Build the dataframe with annotations for each transcript
    myDB$Transcripts <- data.frame(ID = geneRows$HGNC.transcript.name.ID,
    version = geneRows$Transcript.stable.ID.version,
    hgncID = geneRows$HGNC.symbol,
    stableID = geneRows$Transcript.stable.ID,
    start = geneRows$Transcript.start..bp,
    end = geneRows$Transcript.end..bp,
    stringsAsFactors = FALSE)

    return(myDB)

}
```
&nbsp;

###### 4.2.4  Constructing the PDB Chain ID dataframe

&nbsp;

Each Transcript has several corresponding PDB chains and each PDB chain can be mapped to multiple transcripts. Therefore, each entry for the PDB chain must be added to the dataframe in order to map them to the correct transcripts. There may be duplicated entries due to annotations, but they are removed using the `unique()` function. The code for `addpdbChains()` is shown below:

&nbsp;

```R
addpdbChain <- function(myDB, martDF) {

# Add the PDB chains and their corresponding
# transcript IDs to the database
pdbChains <- data.frame(ID = martDF$PDB.ENSP.mappings,
transcriptHGNC =
martDF$HGNC.transcript.name.ID,
stringsAsFactors = FALSE)
pdbChains <- unique(data.table::as.data.table(pdbChains))

myDB$pdbChains <- as.data.frame(pdbChains)

return(myDB)

}
```
&nbsp;

By calling `unique` we remove redundancy of adding several of the same data into the dataframe. The PDB chain dataframe is further extended by adding **Sequence** and **Resolution** annotations retrieved from the PDB database.

In order to retrieve the data from PDB, the functions  `readXML()`, `fetchPDBXML()` and `addPDBdata()` have been developed and are shown below:

&nbsp;

```R

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

}

```
&nbsp;

`readXML()` is a function that takes a `URL` path to an XML file (in this case a XML path to PDB's API services) and retrieves defined nodes and converts the data into a dataframe. `fetchPDBXML()` is the function that calls `readXML()` and is also shown below:

&nbsp;

```R
fetchPDBXML <- function(myDB) {

# Generate a list of PDB chains, remove the
# chain identifier
chainList <- unique(myDB$pdbChains$ID)
chainList <- gsub("\\..*", "", chainList)
myPDBData <- data.frame()

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
```
&nbsp;

`fetchPDBXML()` takes the database being developed and produces a data frame which contains the **PDB ID**, **Resolution** and **Sequence** for each chain in `martDF`.  `addPDBdata` is the function that calls `fetchPDBXML()`  and is shown below:

&nbsp;

```R
addPDBdata <- function(myDB) {

message("Retrieving PDB data...")
# Call fetchPDBXML to retrieve resolutions and sequences
# from PDB
myPDBData <- unique(fetchPDBXML(myDB))

# Set Resolution and sequence
myDB$pdbChains$Resolution <- NA
myDB$pdbChains$Sequences <- NA

# Merge the two data frames using match

sel <- match(myDB$pdbChains$ID, myPDBData$IDs)
myDB$pdbChains$Resolution <- myPDBData[sel,]$Resolution
myDB$pdbChains$Sequences <- myPDBData[sel,]$Sequence

return(myDB)

}
```
Once these function calls are completed, the data table for Transcripts is completed with it's full annotations. 

&nbsp;

###### 4.2.5  Constructing the PDB ID dataframe

The PDB dataframe is quite simple since it simply links each PDB ID with the corresponding PDB ID according to BiomaRt's data. As earlier said, there is a problem where BiomaRt maps a PDB chain with a different PDB ID for several entries. It is evident that some PDB ID's are identical in structure but are given different PDB ID's depending on the environmental and conformational effects. This means that PDB chains with different PDB ID identifiers are most likely from a structurally identical protein. BiomaRt may have done this in order to associate PDB chains with higher resolution structures. Therefore, when constructing the PDB chain table and retrieving resolution data, the PDB chain was used to find the resolution of the structure and not the PDB ID. 

In order to construct the PDB dataframe, the function `addPDB()` is called.

```R
addPDB <- function(myDB, martDF) {

    # Add the PDB IDs and their corresponding
    # PDB chains to the database

    pdbIDs <- data.frame(ID = martDF$PDB.ID,
    chainID = martDF$PDB.ENSP.mappings,
    stringsAsFactors = FALSE)

    # Call the data.table package again in order to find unique entries quickly
    pdbIDs <- unique(data.table::as.data.table(pdbIDs))

    myDB$pdbID <- as.data.frame(pdbIDs)

    return(myDB)

}

```

&nbsp;

After these functions are called, the dataset is completed. We can now test for validation of the dataset using different function calls and tests.

&nbsp;

#### 4.3 Final validation or mapping

Validation and statistics of our mapping:

```R

# do we now have all HGNC symbols mapped?
all(uHGNC %in% myDB$HGNC$ID)  # TRUE

# do we now have all Transcript IDs mapped?
all(uTranscript %in% myDB$Transcripts$ID)  # TRUE

# do we now have all PDB chain IDs mapped?
all(uPDBchains %in% myDB$pdbChains$ID)  # TRUE

# do we now have all HGNC symbols mapped?
all(uHGNC %in% myDB$HGNC$ID)  # TRUE

# How many HGNC symbols did we miss from the original Ensembl Dataset?

length(unique(martDFall$HGNC.symbol)) # 2554
length(unique(myDB$HGNC$ID) #  2554, 100% of the HGNC symbols were mapped

# How many transcripts did we miss?

length(unqiue(martDFall$HGNC.Transcript.name.ID)) #3782
length(unique(myDB$Transcripts$ID)) # 3782, 100% of transcripts were mapped

# How many PDB chains did we miss?

length(unique(martDFall$PDB.ENSP.mappings) # 37329
length(unique(myDB$pdbChains$ID)) # 37329, 100% of PDB chains were mapped

# How many PDB IDs did we miss?

length(unique(martDF2$PDB.ID)) # 16911
length(unique(myDB$pdbID$ID)) # 16911, 100% of PDB ID's were mapped

# This data makes sense since the PDB is consistent with Ensembl data


# Done.
# This concludes construction of the database with mapping and annotations

save(myDB, file = file.path("inst", "extdata", "pdb2sym.RData"))

# From an RStudio project, the file can be loaded with
load(file = file.path("inst", "extdata", "pdb2sym.RData"))


```
&nbsp;

## 5 Biological Validation for Annotations: Multiple Sequence Alignments and Coordinate Checking

For more detailed validation, we need to look at the sequences of the PDB chains and their mapped transcripts and perform a multiple sequence alignment to see if they are identical or not. Then we also need to check the coordinates from Gencode for each PDB chain and the coordinates of each transcript retrieved from Ensembl. We will call the `validate()` function in order to do this.

&nbsp;

```R

# First retrieve the Gencode gene coordinate data for each PDB chain entry by reading the 'gencode.v29.basic.annotation.gff3` file

gffData <- gffFile <- file.path("../data", "gencode.v29.basic.annotation.gff3")
gffData <- rtracklayer::readGFF(gffFile)

validate <- function(transcriptHGNC, chainID, myDB, gffData,
showConsensus = TRUE, showPrint = TRUE) {

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

gffStart <- coordinateData$start
gffEnd <- coordinateData$end

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
```
&nbsp;

Lets test some PDB chains and transcript genes.


```R
# Assume we would like to validate some transcripts for MT-ND2, a HGNC symbol gene

transcripts <- myDB$Transcripts[myDB$Transcripts$hgncID == "MT-ND2",]

# ID           version hgncID        stableID start  end
# 2 MT-ND2-201 ENST00000361453.3 MT-ND2 ENST00000361453  4470 5511

# Let looks at the transcript and it's PDB chains
transcriptHGNC <- transcripts$ID

chains <- myDB$pdbChains[myDB$pdbChains$transcriptHGNC == transcriptHGNC,]

# ID transcriptHGNC
# 3 5xtc.i     MT-ND2-201
# 4 5xtd.i     MT-ND2-201

# Let's use "5xtc.i" as our PDB chain we want to validate
chainID <- chains$ID[1]

validate(transcriptHGNC, chainID, myDB, gffData)

# Console Output
This is the Consensus Sequence
?NPLAQPVIYSTIFAGTLITALSSHWFFTWVGLEMNMLAFIPVLTKKMNPRSTEAAIKYFLTQATASMILLMAILFNNMLSGQWTMTNTTNQYSSLMIMMAMAMKLGMAPFHFWVPEVTQGTPLTSGLLLLTWQKLAPISIMYQISPSLNV?LLLTLSILSIMAGSWGGLNQTQLRKILAYSSITHMGWMMAVLPYNPNMTILNLTIYIILTTTAFLLLNLNSSTTTLLLSRTWNKLTWLTPLIPSTLLSLGGLPPLTGFLPKWAIIEEFTKNNSLIIPTIMATITLLNLYFYLRLIYSTSITLLPMSNNVKMKWQFEHTKPTPFLPTLIALTTLLLPISPFMLMIL
The coordinates are the same for the transcript and the
PDB chain.

# The Consensus sequence shows the conserved residues in a sequence. It seems that there are two
# residues that are different comparing the transcript recieved from Ensembl vs the chain sequence
# from the PDB. This may be due to the fact that we are compraing protein sequences which are
# translated from Amino acid sequences which may cause slight variability.

# Lets look at the other chain matching the transcript

chainID <- chains$ID[2]

validate(transcriptHGNC, chainID, myDB, gffData)

# Console Output
This is the Consensus Sequence
?NPLAQPVIYSTIFAGTLITALSSHWFFTWVGLEMNMLAFIPVLTKKMNPRSTEAAIKYFLTQATASMILLMAILFNNMLSGQWTMTNTTNQYSSLMIMMAMAMKLGMAPFHFWVPEVTQGTPLTSGLLLLTWQKLAPISIMYQISPSLNV?LLLTLSILSIMAGSWGGLNQTQLRKILAYSSITHMGWMMAVLPYNPNMTILNLTIYIILTTTAFLLLNLNSSTTTLLLSRTWNKLTWLTPLIPSTLLSLGGLPPLTGFLPKWAIIEEFTKNNSLIIPTIMATITLLNLYFYLRLIYSTSITLLPMSNNVKMKWQFEHTKPTPFLPTLIALTTLLLPISPFMLMIL
The coordinates are the same for the transcript and the
PDB chain.

# These chains are the exact same sequence. This is valid even though the PDB chains may come from
# different PDB structures, usually those structures are identical but placed in different
# environments. They may also be structure which change conformation due to binding. Nonetheless, both
# chains are mapped to the HGNC transcript symbol and are both valid.

# You may perform validation over any of the genes in the dataset. 
```

&nbsp;

## 6 Annotation of the example gene set

To conclude, we annotate the example gene set, validate the annotation, and store the data in an edge-list format.

&nbsp;

```R

# The specification of the sample set is copy-paste from the 
# Professor Boris Steipe's String database project

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

# W
x <- which((xSet %in% myDB$HGNC$ID))
 xSetPDB <- myDB$Transcripts[match(xSet[x], myDB$Transcripts$hgncID), c("hgncID", "ID")]
 sel <- match(xSetPDB$ID, myDB$pdbChains$transcriptHGNC)
 
 xSetPDB <- cbind(xSetPDB, chainID = 
 myDB$pdbChains[sel, "ID"], Sequence = myDB$pdbChains[sel, "Sequences"],
 Resolution = myDB$pdbChains[sel, "Resolution"])
 sel <- match(xSetPDB$chainID, myDB$pdbID$chainID)
 
 xSetPDB <- cbind(xSetPDB, PDBID = myDB$pdbID[sel, "ID"])

# Save the annotated set

writeLines(c("hgncID\t chainID\t PDBID\t Sequence\t Resolution",
             sprintf("%s\t%s\t%s\t%s\t%s\t%s\n", xSetPDB$hgncID, xSetPDB$ID, xSetPDB$chainID, xSetPDB$PDBID, xSetPDB$Sequence, xSetPDB$Resolution)),
             con = "inst/extdata/xSetPDB.tsv")

# The data set can be read back in again (in an RStudio session) with
myXset <- read.delim(file.path("inst", "extdata", "xSetPDB.tsv"),
stringsAsFactors = FALSE)

# From an installed package, the command would be:
myXset <- read.delim(system.file("extdata",
"xSetPDB.tsv",
package = "BCB420.2019.PDB"),
stringsAsFactors = FALSE)

```

&nbsp;

## 7 References

&nbsp;

Sample code for writing annotating the sample gene set was taken from Professor Boris Steipe's STRING database annotation R package linked [here](https://github.com/hyginn/BCB420.2019.STRING). 

&nbsp;

*  Frankish A. et al. (2018) GENCODE reference annotation for the human and mouse genomes. Nucleic Acids Research, 47: 766-773.

*  H.M. Berman, J. Westbrook, Z. Feng, G. Gilliland, T.N. Bhat, H. Weissig, I.N. Shindyalov, P.E. Bourne.
     (2000) The Protein Data Bank Nucleic Acids Research, 28: 235-242.

&nbsp;

## 8 Acknowledgements

User `Ben Bolker` [posted on Stack](https://stackoverflow.com/questions/7060272/split-up-a-dataframe-by-number-of-rows) how to  split a dataframe by a distinct amount of rows.

&nbsp;

&nbsp;

<!-- [END] -->



% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fetchPDBXML.R
\name{fetchPDBXML}
\alias{fetchPDBXML}
\title{fetchPDBXML.R}
\usage{
fetchPDBXML(myDB)
}
\arguments{
\item{myDB}{A list. The database which is being built.}
}
\value{
myPDBData; a data frame with the PDB data.
}
\description{
\code{fetchPDBXML} Fetch the PDB data including Sequence and Resolution
from PDB's API web services and return a dataframe including
the PDB ID, chain identifier, Sequence and Resolution for each chain
in the database.
}
\details{
Details.
}

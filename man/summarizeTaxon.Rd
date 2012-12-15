\name{summarizeTaxon}
\alias{summarizeTaxon}
\title{Generate A Statistical Summary Of One Genus}
\description{
  Returns the mean, median, and quartiles for each column of one genus. The summary can be generated on the genus as a whole or by species. 
}
\usage{
summarizeTaxon(taxonname, byspecies=FALSE, columns=NA, paramnames=NA, paramvalues=NA)
}
\arguments{
  \item{taxonname}{is either a character literal containing the name of the genus or a data frame with data from the genus. If a character literal is provided, the optional parameters paramnames and paramvalues may also be used.}

  \item{byspecies}{is a boolean value indicating whether the analysis should be performed by species or aggregated on all species across the taxon. If TRUE, a summary for each species will be shown. If FALSE or if no value is provided, one summary aggregating across all species in the genus will be shown.}
    
  \item{columns}{is the vector of columns to include in the summary. All columns must be numeric.}
    
  \item{paramnames}{is a vector of additional, optional parameter names that will be used in the URI. See http://api.neotomadb.org/doc/resources/sampledata for a list of valid parameter names. If no vector is provided, no additional parameters will be used.}

  \item{paramvalues}{is a vector of additional, optional parameter values that correspond to the parameter names in the paramnames vector, which will be used to construct the URI. Must be the same length as paramnames. If a paramvalues vector is included, a paramnames vector must also be included.}
}
\examples{
summarizeTaxon("Mammuthus")
mammuthus <- getSampleData(paramnames="taxonname", paramvalues="Mammuthus")
summarizeTaxon(mammuthus, byspecies=TRUE, columns=c("SiteLatitudeNorth"))
summarizeTaxon("Microtus", columns=c("SiteLatitudeNorth"))
summarizeTaxon("Mammuthus", paramnames=c("ageold", "ageyoung"), 
    paramvalues=c(50000, 0), columns=c("SiteLatitudeNorth"))
}
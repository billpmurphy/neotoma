\name{getTaxon}
\alias{getTaxon}
\title{Transform One Genus Into An R Data Frame}
\description{
  Returns a data frame of one genus based on the genus name and optional additional parameters. 
}
\usage{
getTaxon(taxonname, paramnames=NA, paramvalues=NA)
}
\arguments{
  \item{taxonname}{is a character literal containing the name of the genus to be returned as a data frame.}
    
  \item{paramnames}{is a vector of additional, optional parameter names that will be used in the URI. See http://api.neotomadb.org/doc/resources/sampledata for a list of valid parameter names. If no vector is provided, no additional parameters will be used.}

  \item{paramvalues}{is a vector of additional, optional parameter values that correspond to the parameter names in the paramnames vector, which will be used to construct the URI. Must be the same length as paramnames. If a paramvalues vector is included, a paramnames vector must also be included.}
}
\examples{
getTaxon("Microtus")
getTaxon("Microtus", paramnames="ageyoung", paramvalues=1000)
}
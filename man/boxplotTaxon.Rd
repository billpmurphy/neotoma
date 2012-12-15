\name{boxplotTaxon}
\alias{boxplotTaxon}
\title{Create A Boxplot Of One Column In A Genus}
\description{
  Creates a boxplot in a new window of one numerical column from a genus. The boxplot can aggregate over all the values of one column for the genus, or 
}
\usage{
boxplotTaxon(taxonname, column, byspecies=FALSE, paramnames=NA, paramvalues=NA)
}
\arguments{
  \item{taxonname}{is either a character literal containing the name of the genus or a data frame with data from the genus. If a character literal is provided, the optional parameters paramnames and paramvalues may also be used.}
  
  \item{column}{is a character literal with name of the column to be plotted.}
  
  \item{byspecies}{is a boolean value determining whether one boxplot or a series of boxplots comparing each species within the genus will be generated. If TRUE, a boxplot for each species will be shown. If FALSE or if no value is provided, one boxplot aggregating across all species in the genus will be created.}
    
  \item{paramnames}{is a vector of additional, optional parameter names that will be used in the URI. See http://api.neotomadb.org/doc/resources/sampledata for a list of valid parameter names. If no vector is provided, no additional parameters will be used.}

  \item{paramvalues}{is a vector of additional, optional parameter values that correspond to the parameter names in the paramnames vector, which will be used to construct the URI. Must be the same length as paramnames. If a paramvalues vector is included, a paramnames vector must also be included.}
}
\examples{
boxplotTaxon("Microtus", column="SiteLongitudeEast")
microtusdata <- getTaxon("Microtus")
boxplotTaxon(microtusdata, column="SiteLongitudeEast")
boxplotTaxon(microtusdata, column="SiteLongitudeEast", paramnames="ageold", paramvalues=1000)
}
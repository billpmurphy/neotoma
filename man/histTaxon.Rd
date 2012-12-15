\name{histTaxon}
\alias{histTaxon}
\title{Histogram Of One Column From A Genus}
\description{
  Creates a histogram in a new window of one column from a genus. The number and size of bins is adjustable.  
}
\usage{
histTaxon(taxonname, column, species=NA, breaks="Sturges", paramnames=NA, paramvalues=NA)
}
\arguments{
  \item{taxonname}{is either a character literal containing the name of the genus or a data frame with data from the genus. If a character literal is provided, the optional parameters paramnames and paramvalues may also be used.}
  
  \item{column}{is a character literal with name of the column to be plotted.}
  
  \item{species}{is a character literal with the name of the species to be selected; if no value is entered, the histogram will aggregate across all species in the genus.}
  
  \item{breaks}{one of:
* a vector giving the breakpoints between histogram cells,
* a single number giving the number of cells for the histogram,
* a character string naming an algorithm to compute the number of cells
* a function to compute the number of cells.
In the last three cases the number is a suggestion only.}
      
  \item{paramnames}{is a vector of additional, optional parameter names that will be used in the URI. See http://api.neotomadb.org/doc/resources/sampledata for a list of valid parameter names. If no vector is provided, no additional parameters will be used.}
          
  \item{paramvalues}{is a vector of additional, optional parameter values that correspond to the parameter names in the paramnames vector, which will be used to construct the URI. Must be the same length as paramnames. If a paramvalues vector is included, a paramnames vector must also be included.}
}
\examples{
histTaxon("Microtus", column="SiteAltitude")
microtusdata <- getTaxon("Microtus")
histTaxon(microtusdata, column="SiteAltitude", species="Microtus longicaudus/M. montanus")
histTaxon(microtusdata, column="SiteAltitude", breaks=20)
}
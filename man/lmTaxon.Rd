\name{lmTaxon}
\alias{lmTaxon}
\title{Linear Regression On Two Numeric Columns in a Genus}
\description{
  Creates a scatterplot in a new window of two numeric columns, displays the correlation coeficcient, and plots the linear regression line.  
}
\usage{
lmTaxon(taxonname, x, y, paramnames=NA, paramvalues=NA)
}
\arguments{
  \item{taxonname}{is either a character literal containing the name of the genus or a data frame with data from the genus. If a character literal is provided, the optional parameters paramnames and paramvalues may also be used.}

  \item{x}{is the name of the column that will be the x column of the scatterplot.}

  \item{y}{is the name of the column that will be the y column of the scatterplot.}

  \item{paramnames}{is a vector of additional, optional parameter names that will be used in the URI. See http://api.neotomadb.org/doc/resources/sampledata for a list of valid parameter names. If no vector is provided, no additional parameters will be used.}

  \item{paramvalues}{is a vector of additional, optional parameter values that correspond to the parameter names in the paramnames vector, which will be used to construct the URI. Must be the same length as paramnames. If a paramvalues vector is included, a paramnames vector must also be included.}
}
\examples{
microtusdata <- getTaxon("Microtus")
lmTaxon(microtusdata, "SiteAltitude", "SiteLongitudeEast")
lmTaxon(microtusdata, "SiteAltitude", "SiteLongitudeEast", paramnames="ageyoung", paramvalues=1000) 
}
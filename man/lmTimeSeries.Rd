\name{lmTimeSeries}
\alias{lmTimeSeries}
\title{Linear Time Series Model of One Column in a Genus}
\description{
  Creates a scatterplot in a new window of the timeseries (using the SampleAge column if available, and the mean of SampleAgeYounger and SampleAgeOlder if not), displays the correlation coeficcient, and plots the linear regression line. 
}
\usage{
lmTimeSeries(taxonname, column, paramnames=NA, paramvalues=NA)
}
\arguments{
  \item{taxonname}{is either a character literal containing the name of the genus or a data frame with data from the genus. If a character literal is provided, the optional parameters paramnames and paramvalues may also be used.}
    
  \item{column}{is the name of the column that will be the y column of the scatterplot.}
  
  \item{paramnames}{is a vector of additional, optional parameter names that will be used in the URI. See http://api.neotomadb.org/doc/resources/sampledata for a list of valid parameter names. If no vector is provided, no additional parameters will be used.}

  \item{paramvalues}{is a vector of additional, optional parameter values that correspond to the parameter names in the paramnames vector, which will be used to construct the URI. Must be the same length as paramnames. If a paramvalues vector is included, a paramnames vector must also be included.}
}
\examples{
microtusdata <- getTaxon("Microtus")
lmTimeSeries(microtusdata, "SiteLatitudeNorth")
lmTimeSeries("Microtus", "SiteLatitudeNorth", paramnames="ageyoung", paramvalues=1000)    
}

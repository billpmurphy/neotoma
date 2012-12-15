\name{getSampleData}
\alias{getSampleData}
\title{Transform One SampleData Taxon Into An R Data Frame}
\description{
  Returns a dataframe based on user-entered parameters or a URI. 
}
\usage{
getSampleData(paramnames=NA, paramvalues=NA, uri=NA)
}
\arguments{
  \item{paramnames}{is a vector of parameter names that will be used in the URI. At least one parameter must be given. At least one parameter must be "taxonname" or "taxonids". See http://api.neotomadb.org/doc/resources/sampledata for a list of valid parameter names.}

  \item{paramvalues}{is a vector of parameter values that correspond to the parameter names in the paramnames vector, which will be used to construct the URI.}

  \item{uri}{is the URI to transform into a data frame. This parameter can be used instead of paramnames and paramvalues to select data from NeotomaDB.}
}
\examples{
getSampleData(paramnames="taxonname", paramvalues="Mammuthus")
getSampleData(paramnames=c("taxonname", "ageold", "ageyoung"), c("Mammuthus", 50000, 0))
#getSampleData(uri="http://api.neotomadb.org/data/sampledata?taxonname=Smilodon%")
}
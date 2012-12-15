library("rjson")
require("rjson")
createURI <- function(paramnames=NA, paramvalues=NA) {
    require("rjson")
    if(is.na(paramnames[1]))
        stop("Parameter name not provided")
    if(is.na(paramvalues[1])) 
        stop("Parameter value not provided")
    if(!("taxonname" %in% paramnames) & !("taxonids" %in% paramnames)) 
        stop("taxonname or taxonids required.")
    if(length(paramnames) != length(paramvalues)) 
        stop("Lists of parameter names and values are not the same length.") 

    #create URI
    completeURI <- "http://api.neotomadb.org/data/sampledata?"
    for (i in 1:length(paramnames)) {
        if (i > 1) 
            completeURI <- paste0(completeURI, "&")
        completeURI <- paste0(completeURI, paramnames[i], "=", paramvalues[i])
        if (paramnames[i] == "taxonname")
            completeURI <- paste0(completeURI, "%")
    }
    completeURI
}

getSampleData <- function(paramnames=NA, paramvalues=NA, uri=NA) {
    require("rjson")
    if(is.na(uri) & (is.na(paramnames[1]) | is.na(paramvalues[1])))
        stop("URI or parameter names/values not provided.")

    if (is.na(uri))
        dataURI <- createURI(paramnames=paramnames, paramvalues=paramvalues)
    else
        dataURI <- uri
    
    #retrieve JSON data    
    json <- fromJSON(file=dataURI) #get data, convert from JSON to list of lists
    if (json$success != 1) stop("Error retrieving data from api.neotomadb.org.")
    
    #turn json$data into a data frame
    databasecolumns <- names(json$data[[1]]) #get column names
    ##replace NULL with NA
    json$data <- lapply(json$data, function(x) lapply(x, function(x) ifelse(is.null(x), NA, x))) 
    ##convert to frame
    dataframe <- as.data.frame(t(do.call(cbind, lapply(json$data, unlist, use.names=FALSE))))
    ##apply column names to frame
    colnames(dataframe) <- databasecolumns 
    ##type conversions 
    dataframe <- transform(dataframe,
        SiteLongitudeWest = as.numeric(SiteLongitudeWest),
        SiteLatitudeSouth = as.numeric(SiteLatitudeSouth),
        TaxonName = as.character(TaxonName),
        VariableElement = as.character(VariableElement),
        Value = as.numeric(Value),
        VariableContext = as.character(VariableContext),
        TaxaGroup = as.character(TaxaGroup),
        SampleAgeYounger = as.numeric(SampleAgeYounger),
        SampleAgeOlder = as.numeric(SampleAgeOlder),
        SiteLongitudeEast = as.numeric(SiteLongitudeEast),
        SiteAltitude = as.numeric(SiteAltitude),
        VariableUnits = as.numeric(VariableUnits),
        DatasetID = as.numeric(DatasetID),
        SampleAge = as.numeric(SampleAge),
        SiteLatitudeNorth = as.numeric(SiteLatitudeNorth))
    dataframe
}

getTaxon <- function(taxonname, paramnames=NA, paramvalues=NA) {
    require("rjson")
    if(is.na(paramnames[1]) | is.na(paramvalues[1]))
        getSampleData(paramnames=c("taxonname"), paramvalues=c(taxonname))
    else
        getSampleData(paramnames=c("taxonname", paramnames), paramvalues=c(taxonname, paramvalues))
}

summarizeTaxon <- function(taxonname, byspecies=FALSE, columns=NA, paramnames=NA, paramvalues=NA) {
    require("rjson")
    if(!is.character(taxonname) & !is.data.frame(taxonname)) 
        stop("Type error: Provide a string or data frame.")
    if(!is.logical(byspecies)) 
        stop("Type error: Provide a TRUE or FALSE value for bySpecies.")
    
    #get data
    if(is.character(taxonname))
        dataframe <- getTaxon(taxonname, paramnames=paramnames, paramvalues=paramvalues)
    if(is.data.frame(taxonname))
        dataframe <- taxonname
    
    if(is.na(columns))
        columns <- names(dataframe[,sapply(dataframe, is.numeric) & 
            sapply(dataframe, function(x) !all(is.na(x)))])
    
    #error if no numeric columns or named columns
    if (length(columns) == 0)
        stop("Invalid columns selected")  
        
    if (!byspecies)    
        summary(dataframe[,columns])
    else
        aggregate(dataframe[,columns], by=list(dataframe$TaxonName), FUN=summary)
}

lmTaxon <- function(taxonname, x, y, paramnames=NA, paramvalues=NA) {
    require("rjson")
    if(!is.character(taxonname) & !is.data.frame(taxonname)) 
        stop("Provide a string or data frame for taxonname.")
    
    #get data
    if(is.character(taxonname))
        dataframe <- getTaxon(taxonname, paramnames=paramnames, paramvalues=paramvalues)
    if(is.data.frame(taxonname))
        dataframe <- taxonname

    if (!(x %in% names(dataframe)))
        stop("Invalid x column selected")
    if(!(x %in% names(dataframe[,sapply(dataframe, is.numeric)])))
        stop("Non-numeric x column selected.")
    if(!(x %in% names(dataframe[sapply(dataframe, function(z) !all(is.na(z)))])))
        stop("Null x column selected.")
    if (!(y %in% names(dataframe)))
        stop("Invalid y column selected")
    if(!(y %in% names(dataframe[,sapply(dataframe, is.numeric)])))
        stop("Non-numeric y column selected.")
    if(!(y %in% names(dataframe[sapply(dataframe, function(z) !all(is.na(z)))])))
        stop("Null y column selected.")
    
    genusname = strsplit(na.omit(dataframe$TaxonName)[1], " ")[[1]][1]
    linearmodel <- lm(dataframe[,x] ~ dataframe[,y])
    dev.new()
    plot(x=dataframe[,x], y=dataframe[,y], type="p", xlab=x, ylab=y)
    abline(reg=linearmodel, col="red")
    title(main=paste("Scatterplot of", genusname, x, "and", y), col.main="black", font.main=4)
    mtext(paste("R-squared =", round(summary(linearmodel)$r.squared, digits=6)))
}

lmTimeSeries <- function(taxonname, column, paramnames=NA, paramvalues=NA) {
    require("rjson")
    if(!is.character(taxonname) & !is.data.frame(taxonname)) 
        stop("Provide a string or data frame for taxonname.")
    
    #get data
    if(is.character(taxonname))
        dataframe <- getTaxon(taxonname, paramnames=paramnames, paramvalues=paramvalues)
    if(is.data.frame(taxonname))
        dataframe <- taxonname
    
    if (length(column) > 1 | !is.character(column))
        stop("Invalid column parameter.")    
    if (!(column %in% names(dataframe)))
        stop("Invalid column selected.")
    if(!(column %in% names(dataframe[,sapply(dataframe, is.numeric)])))
        stop("Non-numeric column selected.")
    if(!(column %in% names(dataframe[sapply(dataframe, function(z) !all(is.na(z)))])))
        stop("Null column selected.")

    genusname = strsplit(na.omit(dataframe$TaxonName)[1], " ")[[1]][1]
    timeseries <- ifelse(!is.na(dataframe$SampleAge), dataframe$SampleAge,
        (dataframe$SampleAgeYounger + dataframe$SampleAgeOlder)/2)
    linearmodel <- lm(timeseries ~ dataframe[,column])
    dev.new()
    plot(x=timeseries, y=dataframe[,column], type="p", xlab="Years Ago", ylab=column)
    abline(reg=linearmodel, col="red")
    title(main=paste("Time Series Scatterplot of",genusname,column), 
        col.main="black", font.main=4)
    mtext(paste("R-squared =", round(summary(linearmodel)$r.squared, digits=6)))    
}

histTaxon <-function(taxonname, column, species=NA, breaks="Sturges", paramnames=NA, paramvalues=NA) {
    require("rjson")
    if(!is.character(taxonname) & !is.data.frame(taxonname)) 
        stop("Provide a string or data frame for taxonname.")
    
    #get data
    if(is.character(taxonname))
        dataframe <- getTaxon(taxonname, paramnames=paramnames, paramvalues=paramvalues)
    if(is.data.frame(taxonname))
        dataframe <- taxonname

    if (!(column %in% names(dataframe)))
        stop("Invalid column selected")
    if(!(column %in% names(dataframe[,sapply(dataframe, is.numeric)])))
        stop("Non-numeric column selected.")
    if(!(column %in% names(dataframe[sapply(dataframe, function(x) !all(is.na(x)))])))
        stop("Null column selected.")
    if(!is.na(species) & !(species %in% unique(dataframe$TaxonName)))
        stop("Invalid species")

    if(is.na(species)) {
        genusname = strsplit(na.omit(dataframe$TaxonName)[1], " ")[[1]][1]
        dev.new()
        hist(dataframe[,column], xlab=column, main="", breaks=breaks)
        title(main=paste(column,"of",genusname), col.main="black", font.main=4)
    }
    else {
        dev.new()
        hist(dataframe[dataframe$TaxonName == species,column], xlab=column, main="", breaks=breaks)
        title(main=paste(species, column), col.main="black", font.main=4)
    }
}

boxplotTaxon <- function(taxonname, column, byspecies=FALSE, paramnames=NA, paramvalues=NA) {
    require("rjson")
    if(!is.character(taxonname) & !is.data.frame(taxonname)) 
        stop("Provide a string or data frame for taxonname.")
    if(!is.logical(byspecies)) 
        stop("Provide a TRUE or FALSE value for bySpecies.")
            
    #get data
    if(is.character(taxonname))
        dataframe <- getTaxon(taxonname, paramnames=paramnames, paramvalues=paramvalues)
    if(is.data.frame(taxonname))
        dataframe <- taxonname
        
    if (!(column %in% names(dataframe)))
        stop("Invalid column selected")
    if(!(column %in% names(dataframe[,sapply(dataframe, is.numeric)])))
        stop("Non-numeric column selected.")
    if(!(column %in% names(dataframe[sapply(dataframe, function(x) !all(is.na(x)))])))
        stop("Null column selected.")
      
    #create plot
    genusname = strsplit(na.omit(dataframe$TaxonName)[1], " ")[[1]][1]
    if(!byspecies) { 
        dev.new()
        boxplot(dataframe[,column], col="lightblue", ylab=column)
        title(main=paste("Genus",genusname), col.main="black", font.main=4)
    }
    else {
        specieslist <- list()
        for (i in unique(dataframe$TaxonName)) 
            specieslist [[i]] <- dataframe[dataframe$TaxonName==i, column]
        dev.new()
        boxplot(specieslist, col="lightblue", ylab=column, 
            names=1:length(unique(dataframe$TaxonName)))
        title(main=paste("Genus",genusname, "by Species"), col.main="black", font.main=4)        
        
        #put legend in new window
        leg <- character()
        leg[1] <- "Legend:"
        for(i in 1:length(unique(dataframe$TaxonName))) 
            leg[i+1]<- paste0(i,": ",dataframe$TaxonName[i])
        dev.new()
        plot.new()
        legend(x="center", legend=leg, cex=1)

    }
}




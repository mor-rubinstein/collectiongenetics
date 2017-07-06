##------------------------------##
####   Classes                ####
##------------------------------##

#' A S4 class for a matrix of genotype calls.
#'
#' A matrix of genotype calls, with each row representing a marker, and each
#' column representing an individual sample. The genotype call is in the form
#' of a character.
#'
#' @slot Calls matrix: The matrix of genotype calls
#' @slot Samples character: Individual samples IDs of the
#' @slot Markers character: Markers IDs
#'
setClass(Class = "GenotypeCalls", representation(Calls = "matrix", Samples = "character",
    Markers = "character"))

##-----------------------------------------------------------------------------
#' A S4 class for sample information
#'
#' A data frame which row names are the sample IDs. It can include any sample
#' information available. Usually it will have a column with the cultivar name.
#'
#' @slot SInfo data.frame: The table of sample info
#' @slot Samples character: Sample IDs. identical to the row names of SInfo.
#'
setClass(Class="SampleInfo",
         representation(SInfo="data.frame", Samples = "character")
)

##-----------------------------------------------------------------------------
#' A S4 class for marker information
#'
#' A data frame which row names are the markers IDs. It can include any marker
#' information available, such as the two alleles, genome positions if
#' available, and the quality ranks given by the Fluidigm team.
#'
#' @slot MInfo data.frame: The table of marker info
#' @slot Markers character: Marker IDs. identical to the row names of MInfo.
#'
setClass(Class="MarkerInfo",
         representation(MInfo="data.frame",Markers = "character")
)

##-----------------------------------------------------------------------------
#' A data frame for documenting the filtering process
#'
#' A data frame with the name of each filtering step, and the numbers of
#' markers and samples following that step. The first row in the table shows
#' the numbers of markers and samples in the raw, unfiltered data.
#'
#' @slot filter_table data.frame: The table of filtration information
#'
setClass(Class="filterTable",
         representation(filter_table="data.frame")
)

##-----------------------------------------------------------------------------
#' A S4 class for a Genotype database
#'
#' A class that typically represents a genotyping assay. It contains (by
#' inharitance) a GenotypeCalls object, which is the genotype call matrix, a
#' SampleInfo and a MarkerInfo objects, which hold the information about the
#' samples and the markers in the assay, and the filterTable, which documents
#' the cleaning of the assay data.
#'
#' @slot MarkerFilter integer: The indeces of markers that were not masked by
#' the filtration process.
#' @slot SampleFilter integer: The indeces of samples that were not masked by
#' the filtration process.
setClass(Class="Genotype",
         representation(MarkerFilter= "numeric",
                        SampleFilter= "numeric"
         ),
         contains=c("GenotypeCalls","SampleInfo","MarkerInfo", "filterTable")
)

##-----------------------------------------------------------------------------
##------------------------------##
####   Show methods           ####
##------------------------------##

setMethod("show","GenotypeCalls",function(object){
  if (SampleNo(object)==0)
    cat("No genotype calls matrix is available \n")
  else{
    cat("Genotype calls object:\n")
    cat("Number of sample: ", SampleNo(object),"\n")
    cat("Samples: ", head(SampleNames(object)),"\n")
    cat("Number of markers: ", MarkerNo(object),"\n")
    cat("Markers: ", head(MarkerNames(object)),"\n")
    head(object@Calls)
    cat("\n")
  }
})
setMethod("show","MarkerInfo",function(object){
  if (nrow(object@MInfo)==0)
    cat("No markers information is available \n")
  else{
    cat("Markers object:\n")
    cat("Number of markers: ", nrow(object@MInfo),"\n")
    cat("Markers: ", head(rownmames(object@MInfo)),"\n")
    cat("\nMarkers' information: \n")
    for (Cat in colnames(object@MInfo)){
      cat(Cat,": ",head(unique(object@MInfo[,Cat])),"\n")
    }
    cat("\n")
  }
})
setMethod("show","SampleInfo",function(object){
  if (nrow(object@SInfo)==0)
    cat("No sample information is available \n")
  else{
    cat("Samples object:\n")
    cat("Number of sample: ", nrow(object@SInfo),"\n")
    cat("Samples: ", head(rownames(object@SInfo)),"\n")
    cat("\nSamples' information: \n")
    for (Cat in colnames(object@SInfo)){
      cat(Cat,": ",head(unique(object@SInfo[,Cat])),"\n")
    }
    cat("\n")
  }
})
setMethod("show","filterTable",function(object) {
  if (nrow(object@filter_table)<2)
    cat("No Filtration has been made.\n")
  else{
    cat("Filtration:\n")
    cat(object@filter_table)
    cat("\n")
  }
})
setMethod("show","Genotype",function(object){
  if (SampleNo(object)==0)
    cat("No genotype calls matrix is available \n")
  else{
    cat("Genotype calls object:\n")
    cat("Number of sample: ", SampleNo(object),"\n")
    cat("Samples: ", head(SampleNames(object)),"\n")
    cat("Number of markers: ", MarkerNo(object),"\n")
    cat("Markers: ", head(MarkerNames(object)),"\n")
    head(object@Calls)
    cat("\n")
  }
  if (nrow(object@SInfo)==0)
    cat("No samples information is available \n")
  else{
    cat("\nSamples' information: \n")
    for (Cat in colnames(object@SInfo)){
      cat(Cat,": ",head(unique(object@SInfo[object@SampleFilter,Cat])),"\n")
    }
    cat("\n")
  }
  if (nrow(object@MInfo)==0)
    cat("No markers information is available \n")
  else{
    cat("\nMarkers' information: \n")
    for (Cat in colnames(object@MInfo)){
      cat(Cat,": ",head(unique(object@MInfo[object@MarkerFilter,Cat])),"\n")
    }
    cat("\n")
  }

})

##------------------------------##
####   Basic functions        ####
##------------------------------##

#' Marker information
#'
#' @param object. The object to get the markers information from.
#' @param Masked logical. Whether to avoid getting the information from masked
#' markers too
#' @return \code{MarkerInfo}
#' @examples
#' data("Mango_raw")
#' MI = getMarkerInfo(MxS)
#' head(MI)
setGeneric("getMarkerInfo",function(object,Masked=TRUE) {
  standardGeneric("getMarkerInfo")
  })
setMethod("getMarkerInfo","MarkerInfo",function(object) {
  object@MInfo
})
setMethod("getMarkerInfo","Genotype",function(object,Masked=TRUE) {
  if(Masked)
    callNextMethod()[object@MarkerFilter, ,drop=FALSE]
  else
    callNextMethod()
})

##-----------------------------------------------------------------------------

#' The number of markers
#'
#' @param object. The object to count the markers from
#' @param Masked logical. Whether to ignore masked markers
#' @return \code{numeric} the number of markers
#' @examples
#' data("Mango_raw")
#' MarkerNo(MxS)
#' MarkerNo(MxS, Masked=F)

setGeneric("MarkerNo",function(object, Masked=TRUE)  {
  standardGeneric("MarkerNo")
  })
setMethod("MarkerNo","GenotypeCalls",function(object,Masked=TRUE) {
  length(object@Markers)
})
setMethod("MarkerNo","MarkerInfo",function(object,Masked=TRUE) {
  length(object@Markers)
})
setMethod("MarkerNo","Genotype",function(object,Masked=TRUE) {
  if(Masked)
    length(object@MarkerFilter)
  else
    nrow(object@Calls)
})

##-----------------------------------------------------------------------------

#' The marker names
#'
#' @param object. The object to get the markers from
#' @param Masked logical. Whether to seclude masked markers
#' @return \code{character} the marker names
#' @examples
#' data("Mango_raw")
#' MarkerNames(MxS)
#' MarkerNames(MxS, Masked=F)

setGeneric("MarkerNames",function(object, Masked=TRUE) {
  standardGeneric("MarkerNames")
  })
setMethod("MarkerNames","GenotypeCalls",function(object,Masked=TRUE) {
  object@Markers
})
setMethod("MarkerNames","MarkerInfo",function(object,Masked=TRUE) {
  object@Markers
})
setMethod("MarkerNames","Genotype",function(object,Masked=TRUE) {
  if(Masked)
    object@Markers[object@MarkerFilter]
  else
    object@Markers
})

##-----------------------------------------------------------------------------

#' Reorder the markers
#'
#' @param object \code{Genotype}. The object to reorder the markers in
#' @param Order \code{numeric} or \code{character}. A numeric vector with the
#' order of markers, or a character vector with the marker names ordered as
#' desired.
#' @return The object with the markers re-ordered as specified.
#' @examples
#' data("Mango_raw")
#' a = MarkerNames(MxS)
#' b = MarkerNames(SetMarkerOrder(MxS, Order=sample(MarkerNo(MxS))))
#' head(a)
#' head(b)
#' setequal(a,b)

setGeneric("SetMarkerOrder",function(object,Order=NULL) {
  standardGeneric("SetMarkerOrder")
  })
setMethod("SetMarkerOrder","Genotype", function(object, Order=NULL) {
  Masked = (MarkerNo(object)>= MarkerNo(object, Masked=F) | MarkerNo(object, Masked=F) != length(Order))
  if (length(Order) != MarkerNo(object, Masked=Masked))
    stop("Order must have the same length as the number of markers.")
  m = mode(Order)
  if (m == "numeric") {
    if (!setequal(Order, 1:MarkerNo(object, Masked=Masked)))
      stop("Order must be a set of numbers between 1 and the number of markers.")
  } else {
    if (m == "character") {
      if (!setequal(Order, MarkerNames(object, Masked=Masked)))
        stop("Order must contain all and only the names of the marker.")
      Order = order(MarkerNames(object, Masked=Masked))[rank(Order)]
    } else stop("Order must be an integer or the set of marker names.")
  }
  if (Masked) {
    o = 1:length(object@Markers)
    o[object@MarkerFilter] = o[object@MarkerFilter][Order]
  } else o = Order
  object@MInfo = object@MInfo[o,]
  object@Markers = object@Markers[o]
  object@Calls = object@Calls[o, ]
  object
})

##-----------------------------------------------------------------------------
#' Change the marker names
#'
#' @param object \code{Genotype}. The object to change the markers names in.
#' @param marker.name code{character}. A character vector with the new marker
#' names.
#' @return The object with the markers having their new names.

setGeneric("SetMarkerName",function(object,marker.name=character(), Masked=F) {
  standardGeneric("SetMarkerName")
  })
setMethod("SetMarkerName","Genotype",function(object, marker.name=character(), Masked=F) {
  if (length(marker.name) != MarkerNo(object, Masked=Masked))
    stop ("The number of new marker names must be the same as the number of markers in the genotype object.")
  if (Masked) {
    object@Markers[object@MarkerFilter] = marker.name
    rownames(object@MInfo)[object@MarkerFilter] = marker.name
  } else {
    object@Markers = marker.name
    rownames(object@MInfo) = marker.name
  }
  object
})

##-----------------------------------------------------------------------------

#' Sample information
#'
#' @param object. The object to get the samples information from.
#' @param Masked logical. Whether to avoid getting the information from masked
#'  samples
#' @return \code{SampleInfo}
#' @examples
#' data("Mango_raw")
#' SI = getSampleInfo(MxS)
#' head(SI)
setGeneric("getSampleInfo",function(object,Masked=TRUE) {
  standardGeneric("getSampleInfo")
})
setMethod("getSampleInfo","SampleInfo",function(object,Masked=TRUE) {
  object@SInfo
})
setMethod("getSampleInfo","Genotype",function(object,Masked=TRUE) {
  if(Masked)
    callNextMethod()[object@SampleFilter,, drop=FALSE]
  else
    callNextMethod()
})
##-----------------------------------------------------------------------------

#' The number of samples
#'
#' @param object. The object to count the samples from
#' @param Masked logical. Whether to ignore masked samples
#' @return \code{numeric} the number of samples
#' @examples
#' data("Mango_raw")
#' SampleNo(MxS)
#' SampleNo(MxS, Masked=F)

setGeneric("SampleNo",function(object, Masked=TRUE) {
  standardGeneric("SampleNo")
})
setMethod("SampleNo","SampleInfo",function(object,Masked=TRUE) {
  length(object@Samples)
})
setMethod("SampleNo","Genotype",function(object,Masked=TRUE) {
  if(Masked)
    length(object@SampleFilter)
  else
    ncol(object@Calls)
})
setMethod("SampleNo","GenotypeCalls",function(object,Masked=TRUE) {
  length(object@Samples)
})

##-----------------------------------------------------------------------------

#' The sample names
#'
#' @param object. The object to get the samples from
#' @param Masked logical. Whether to seclude masked samples
#' @return \code{character} the sample names
#' @examples
#' data("Mango_raw")
#' SampleNames(MxS)
#' SampleNames(MxS, Masked=F)

setGeneric("SampleNames",function(object, Masked=TRUE) {
  standardGeneric("SampleNames")
})
setMethod("SampleNames","GenotypeCalls",function(object,Masked=TRUE){
  object@Samples
})
setMethod("SampleNames","SampleInfo",function(object,Masked=TRUE) {
  object@Samples
})
setMethod("SampleNames","Genotype",function(object,Masked=TRUE) {
  if(Masked)
    object@Samples[object@SampleFilter]
  else
    object@Samples
})
##-----------------------------------------------------------------------------

#' Reorder the samples
#'
#' @param object \code{Genotype}. The object to reorder the sample in
#' @param Order \code{numeric} or \code{character}. A numeric vector with the
#' order of samples, or a character vector with the sample names ordered as
#' desired.
#' @return The object with the samples re-ordered as specified.
#' @examples
#' data("Mango_raw")
#' a = SampleNames(MxS)
#' b = SampleNames(SetSampleOrder(MxS, Order=sample(SampleNo(MxS))))
#' head(a)
#' head(b)
#' setequal(a,b)

setGeneric("SetSampleOrder",function(object,Order=NULL) {
  standardGeneric("SetSampleOrder")
})
setMethod("SetSampleOrder","Genotype", function(object, Order=NULL) {
  if (length(Order) != SampleNo(object))
    stop("Order must have the same length as the number of samples.")
  s = mode(Order)
  if (s == "numeric") {
    if (!setequal(Order, 1:SampleNo(object)))
      stop("Order must be a set of numbers between 1 and the number of samples.")
  } else {
    if (s == "character") {
      if (!setequal(Order, SampleNames(object)))
        stop("Order must contain all and only the names of the sample.")
      Order = order(SampleNames(object))[rank(Order)]
    } else stop("Order must be an integer or the set of sample names.")
  }
  #   object@SampleFilter = object@SampleFilter[Order]
  o = 1:length(object@Samples)
  o[object@SampleFilter] = o[object@SampleFilter][Order]
  object@SInfo = object@SInfo[o,]
  object@Samples = object@Samples[o]
  object@Calls = object@Calls[, o]
  object
})

##-----------------------------------------------------------------------------

#' Get a matrix of genotype calls
#'
#' @param object. The object to get the genotype calls from
#' @param MaskedMarker logical. Whether to seclude masked markers
#' @param MaskedSample logical. Whether to seclude masked samples
#' @return \code{matrix} of genotype calls:
#' \describe{
#' \item{XX}{homozygous to the 1st allele}
#' \item{XY}{heterozygous}
#' \item{YY}{homozygous to the 2nd allele}
#' \item{NA}{genotype call not available}
#' }
#' @examples
#' data("Mango_raw")
#' a = getCalls(MxS)
#' SampleNames(MxS, Masked=F)

setGeneric("getCalls",function(object,MaskedMarker=TRUE,MaskedSample=TRUE) {
  standardGeneric("getCalls")
  })
setMethod("getCalls","GenotypeCalls",function(object,MaskedMarker=TRUE) {
  M = object@Calls
  rownames(M) = object@Markers
  colnames(M) = object@Samples
  M
})
setMethod("getCalls","Genotype",function(object,MaskedMarker=TRUE,MaskedSample=TRUE){
  if (MaskedMarker & MaskedSample)
    callNextMethod()[object@MarkerFilter,object@SampleFilter]
  else {
    if (MaskedMarker)
      callNextMethod()[object@MarkerFilter, ]
    else {
      if (MaskedSample)
        callNextMethod()[, object@SampleFilter]
      else
        callNextMethod()
    }
  }
})

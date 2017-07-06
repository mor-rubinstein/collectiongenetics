#' Import a table of genotype calls.
#'
#' Import a table of genotype calls and create a \code{Genotype} object.
#'
#' @param file character: the name of the file to read the data from.
#' @return a \code{Genotype} object
#'
#' @examples
#' MxS = read.genotype("Fluidigm_Output.csv")
setGeneric("read.genotype", function(file, object) {
  standardGeneric("read.genotype")
})
setMethod("read.genotype", signature(file="character"),
          function(file) {
  read.genotype(file, new("Genotype"))
})
setMethod("read.genotype", signature(file="character", object="GenotypeCalls"),
          function(file, object) {
  M = read.delim(file, stringsAsFactors = FALSE, na.strings = c("No Call",
        "invalid", "Invalid", "NA"), check.names = FALSE)
    if (length(M) == 1)
        M = read.delim(file, stringsAsFactors = FALSE, sep = ",", na.strings =
                         c("No Call", "invalid", "Invalid", "NA"),
                       check.names = FALSE)
    CallCols = apply(M, 2, function(x) any(c("XX", "YY", "XY", NA) %in% x))
    if (length(grep("^X", colnames(M[, CallCols]))) > 0)
        colnames(M[, CallCols]) <- gsub("^X", "", colnames(M[, CallCols]))
    object@Calls = as.matrix(M[, CallCols])
    object@Samples = colnames(M[, CallCols])
    object@Markers = as.character(M[, 1])
    object
})
setMethod("read.genotype",signature(file="character", object="Genotype"),
          function(file, object) {
  object<-callNextMethod()
  object = ResetMarkerFilter(object)
  object = ResetSampleFilter(object)
  object
})

##-------------------------------------------------------------------------------
#' Import the sample data.
#'
#' Import a table of sample data and create a \code{SampleInfo} object.
#'
#' @param object the object (typically of class \code{Genotype}) to add the
#' sample data into.
#' @param file character: the name of the file to read the sample data from.
#'    The 1st column in the file should contain the sample IDs.
#' @return the object containing the sample data.
#'
setGeneric("read.samples",function(object,file) {standardGeneric("read.samples")})
setMethod("read.samples","SampleInfo",function(object,file=character()){
  M = read.delim(file,row.names=1,stringsAsFactors=FALSE,check.names=FALSE,quote="")
  if (length(M) == 0)
    M = read.delim(file,row.names=1,stringsAsFactors=FALSE,check.names=FALSE, sep=",")
  if (SampleNo(object)==0) object@Samples<-rownames(M)
  n = sapply(M, function(x) sum(SampleNames(object) %in% x))
  i = which.max(n)
  if (n[i] > sum(SampleNames(object) %in% row.names(M))) row.names(M) = M[[i]]
  M = M[row.names(M) %in% object@Samples, ,drop=FALSE]
  x = setdiff(object@Samples, row.names(M))
  if (length(x)) {
    warning("Not all the samples in the Genotype object appear in the sample data")
    m = M[1:length(object@Samples), ,drop=FALSE]
    row.names(m)[(nrow(M)+1):nrow(m)] = x
    M = m
  }
  o = order(row.names(M))[rank(object@Samples)]
  M = M[o, ,drop=FALSE]
  object@SInfo= M
  object
})
setMethod("read.samples","Genotype",function(object,file=character()){
  object<-callNextMethod()
  if (SampleNo(object)!=nrow(getSampleInfo(object)))
    stop("Number of row of sample info table must be equal to the number of columns of genotype table  \n","Number of samples:\t",SampleNo(object),"\nSample information table size:\t",nrow(getSampleInfo(object))) else
      if(!all(SampleNames(object)%in%rownames(object@SInfo))) stop(cat(setdiff(SampleNames(object),rownames(object@SInfo)),sep="\n")) else
        object@SInfo = object@SInfo[SampleNames(object), , drop=FALSE]
      object = ResetSampleFilter(object)
      object
})

##-------------------------------------------------------------------------------
#' Import the marker data.
#'
#' Import a table of marker data and create a \code{MarkerInfo} object.
#'
#' @param object the object (typically of class \code{Genotype}) to add the
#' marker data into.
#' @param file character: the name of the file to read the marker data from.
#'    The 1st column in the file should contain the marker IDs.
#' @return the object containing the marker data.
#'
setGeneric("read.markers",function(object,file) {standardGeneric("read.markers")})
setMethod("read.markers","MarkerInfo",function(object,file=character()){
  M=read.delim(file,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)
  if (length(M) ==0)
    M=read.delim(file,row.names=1,stringsAsFactors=FALSE,check.names=FALSE, sep=",")
  if (length(object@Markers)==0) object@Markers<-rownames(M)
  object@MInfo=M
  object
})
setMethod("read.markers","Genotype",function(object,file=character()){
  object <- callNextMethod()
  if (MarkerNo(object)!=nrow(getMarkerInfo(object)))
    stop("Number of row of marker info table must be equal to the number of rows of genotype table  \n","Number of markers:\t",MarkerNo(object),"\nMarker information table size:\t",nrow(getMarkerInfo(object))) else
      if (!all(MarkerNames(object)%in%rownames(object@MInfo))) stop(cat(setdiff(MarkerNames(object),rownames(object@MInfo)),seop="\n")) else
        object@MInfo = object@MInfo[MarkerNames(object), ,drop=FALSE]
      object = ResetMarkerFilter(object)
      object
})

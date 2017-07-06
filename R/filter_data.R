#' Filtering the markers
#'
#' Filter (or 'mask') a specified set of markers from a \code{Genotype} object.
#' Update the filter table and the marker info.
#'
#' @param object \code{Genotype}. The object to be filtered.
#' @param Filter character, integer or logical. The markers to be filtered.
#'  This could be a vector of marker names, marker indeces, or a logical vector
#'  with TRUE values for markers to filter and FALSE for markers to keep.
#' @param print.n.markers logical. Should the number of remaining markers be
#' printed on the screen?
#' @param step_name character. The name of the filtration step. This will
#' in the step name column in the filter table.
#' @return The filtered \code{Genotype} object. The filter table holds a new
#' line with the number of markers left after this step, and the
#' 'Filtered_at_step' column in the MarkerInfo of the filtered markers, have
#' the step name.
#'
#' @examples
#' data("Mango_raw")
#' MarkerNo(MxS)
#' MxS = SetMarkerFilter(MxS, Filter=c("Contig1987_1443", "Contig2006_1970",
#'     "mango_c15408_840"), step_name="arbitrary filtering")
#' getMarkerInfo(MxS, Masked=F)["Contig1987_1443", ]
setGeneric("SetMarkerFilter",function(object,Filter=NULL,print.n.markers=T,
                                      step_name="marker filter") {
  standardGeneric("SetMarkerFilter")
  })
setMethod("SetMarkerFilter","Genotype",
          function(object,Filter=NULL, print.n.markers=T, step_name=
                     "Filter markers") {
  cls = class(Filter)
  object@MarkerFilter <- switch(cls,
                                NULL= stop("Filter is mandatory"),
                                integer = Filter[Filter%in%(1:MarkerNo(object,Masked=FALSE))],
                                logical = object@MarkerFilter[!Filter],
                                character = object@MarkerFilter[!MarkerNames(object)%in%Filter],
                                stop("Filter must be either characater or logical")
  )
  if (print.n.markers) cat("Number of remaining markers: ",length(object@MarkerFilter),"\n")
  MarkerInfo = getMarkerInfo(object, Masked=F)
  if (!("Filtered_at_step" %in% names(MarkerInfo))) {
    if (length(MarkerInfo)) MarkerInfo$Filtered_at_step = character(nrow(MarkerInfo)) else
      MarkerInfo = data.frame(Filtered_at_step = character(MarkerNo(object, Masked=F)), row.names=MarkerNames(object, Masked=F), stringsAsFactors=F)
  }
  f = setdiff(MarkerNames(object, Masked=F), MarkerNames(object, Masked=T))
  i = which(nchar(MarkerInfo[f, "Filtered_at_step"])>0)
  if (length(i)) f = f[-i]
  MarkerInfo[f, "Filtered_at_step"] = step_name
  object@MInfo = MarkerInfo
  if (nrow(getFilterTable(object)) ==0) object = start.filter.table(object)
  UpdateFilterTable(object, step_name=step_name)
})
##---------------------------------------------------------
#' UnFiltering the markers
#'
#' \code{ResetMarkerFilter} reset marker filter, so that all markers are
#'  unmasked
#'
#' @param object \code{Genotype}. The object to be un-filtered.
#' @return The unfiltered \code{Genotype} object
#' @examples
#' data("Mango_raw")
#' MarkerNo(MxS)
#' MxS = SetMarkerFilter(MxS, Filter=c("Contig1987_1443", "Contig2006_1970",
#'     "mango_c15408_840"), step_name="arbitrary filtering")
#' getMarkerInfo(MxS, Masked=F)["Contig1987_1443", ]
#' MxS = ResetMarkerFilter(MxS)
setGeneric("ResetMarkerFilter", function(object){
  standardGeneric("ResetMarkerFilter")})
setMethod("ResetMarkerFilter","Genotype", function(object){
  object@MarkerFilter = 1:MarkerNo(object,Masked=FALSE)
  cat("Number of markers: ",length(object@MarkerFilter),"\n")
  UpdateFilterTable(object, step_name="Reset marker filter")
})

##---------------------------------------------------------
#' Getting the marker Filter
#'
#' @param object \code{Genotype}. The object to get the filter from.
#' @return The indeces of the unmasked markers.
setGeneric("GetMarkerFilter", function(object){
  standardGeneric("GetMarkerFilter")})

setMethod("GetMarkerFilter","Genotype", function(object){
  object@MarkerFilter
})

##---------------------------------------------------------
#' Quality filtering of markers
#'
#' A general function for several standard quality filtrations of markers
#'
#' @param object \code{Genotype}. The object to be filtered.
#' @param quality.criterion chatacter. The name of the quality criterion to
#' filter the markers by. One of the followings:
#' \describe{
#'   \item{no.calls}{filter markers with relative number of no calls > cutoff.}
#'   \item{pic}{filter markers with Polymorphic Information Content < cutoff.}
#'   \item{major.genotype}{filter markers with relative number of samples
#'   having the same call > cutoff. This is important to prevent a situation
#'   where almost all samples are heterozygous.}
#'   \item{minor.genotype}{filter markers where (at least) one of the 3
#'   possible genotype calls has relative number < cutoff. This is important
#'    when we want to make sure all 3 possible genotype calls exist for this
#'     marker.}
#' }
#'
#' @param cutoff \code{numeric}. The cutoff for the filtration
#' @param print.n.markers logical. Should the number of remaining markers be
#' printed on the screen?
#' @return The filtered \code{Genotype} object, with updated filter table and
#' marker info table.
#'
#' @examples
#' data("Mango_raw")
#' MarkerNo(MxS)
#' MxS = FilterMarkerByQuality(MxS, quality.criterion="no.calls", cutoff=0.1)
#' MxS = FilterMarkerByQuality(MxS, quality.criterion="pic", cutoff=0.1)
#' Note: one marker is all heterozygous:
#' genotype_freq(AlleleCounting(MxS))["Contig4311_367", , drop=FALSE]
#' MxS = FilterMarkerByQuality(MxS, quality.criterion="major.genotype",
#' cutoff=0.9)
#' "Contig4311_367" %in% MarkerNames(MxS)
#' getMarkerInfo(MxS, Masked=F)["Contig4311_367", ]
#' genotype_freq(AlleleCounting(MxS))["Contig1309_880", , drop=FALSE]
#' MxS = FilterMarkerByQuality(MxS, quality.criterion="minor.genotype",
#' cutoff=1)
#' getMarkerInfo(MxS, Masked=F)["Contig1309_880", ]

setGeneric("FilterMarkerByQuality", function(object, quality.criterion="",
                                             cutoff=NULL, print.n.markers=T) {
  standardGeneric("FilterMarkerByQuality")})
setMethod("FilterMarkerByQuality", "Genotype",
          function(object, quality.criterion="", cutoff=NULL,
                   print.n.markers=T) {
  return(switch(quality.criterion,
                no.calls = FilterMarkerByNoCall(object, cutoff, print.n.markers=print.n.markers),
                pic = FilterMarkerByPIC(object, cutoff, print.n.markers=print.n.markers),
                major.genotype = FilterMarkerByMajorGenotype(object, cutoff, print.n.markers=print.n.markers),
                minor.genotype = FilterMarkerByMinorGenotype(object, cutoff, print.n.markers=print.n.markers)
  ))
})
##---------------------------------------------------------
#' Filtering Markers by relative number of no calls
#'
#' @param object \code{Genotype}. The object to be filtered.
#' @param cutoff \code{numeric}. Markers with a relative number of missing calls
#'  > cutoff will be masked.
#' @param print.n.markers logical. Should the number of remaining markers be
#' printed on the screen?
#' @return a filtered \code{Genotype} object, with updated filter table and
#' \code{MarkerInfo$Filtered_at_step}
#' @examples
#' data("Mango_raw")
#' MarkerNo(MxS)
#' MxS = FilterMarkerByNoCall(MxS, cutoff=0.1)

setGeneric("FilterMarkerByNoCall",
           function(object, cutoff=NULL, print.n.markers=T) {
             standardGeneric("FilterMarkerByNoCall")
             })
setMethod("FilterMarkerByNoCall", "Genotype",
          function(object, cutoff=0.1, print.n.markers=T) {
  l = SampleNo(object)
  Calls = getCalls(object)
  NoCalls = rowSums(is.na(Calls))
  x = NoCalls/l > cutoff
  cutoff = round(cutoff*100)/100
  step_name = paste("Filtering out SNPs with >", cutoff,  "no calls")
  SetMarkerFilter(object, Filter=x, step_name=step_name, print.n.markers=print.n.markers)
})

##---------------------------------------------------------
#' Filtering Markers by their polymorphic information content
#'
#' @param object \code{Genotype}. The object to be filtered.
#' @param cutoff \code{numeric}. Markers with a PIC < cutoff will be masked.
#' @param print.n.markers logical. Should the number of remaining markers be
#' printed on the screen?
#' @return a filtered \code{Genotype} object, with updated filter table and
#' \code{MarkerInfo$Filtered_at_step}
#' @examples
#' data("Mango_raw")
#' MarkerNo(MxS)
#' MxS = FilterMarkerByPIC(MxS, cutoff=0.1)

setGeneric("FilterMarkerByPIC",
           function(object, cutoff=NULL, print.n.markers=T) {
             standardGeneric("FilterMarkerByPIC")
             })

setMethod("FilterMarkerByPIC", "Genotype",
          function(object, cutoff=0.1, print.n.markers=T) {
  MarkersStat = AlleleCounting(object)
  pic = PIC(MarkersStat,type=1)
  cutoff = round(cutoff*100)/100
  step_name = paste("Filtering out markers with PIC <", cutoff)
  SetMarkerFilter(object, Filter=pic<cutoff, step_name=step_name, print.n.markers=print.n.markers)
})

##---------------------------------------------------------
#' Filtering markers by the relative aboundance of the most aboundant genotype
#'
#' @param object \code{Genotype}. The object to be filtered.
#' @param cutoff \code{numeric}. Markers where relative number of calls having
#' the same genotype call is > \code{cutoff} will be masked.
#' @param print.n.markers logical. Should the number of remaining markers be
#' printed on the screen?
#' @return a filtered \code{Genotype} object, with updated filter table and
#' \code{MarkerInfo$Filtered_at_step}
#' @examples
#' data("Mango_raw")
#' MarkerNo(MxS)
#' genotype_freq(AlleleCounting(MxS))["Contig4311_367", , drop=FALSE]
#' MxS = FilterMarkerByMajorGenotype(MxS, cutoff=0.95)
#' "Contig4311_367" %in% MarkerNames(MxS)
#' getMarkerInfo(MxS, Masked=F)["Contig4311_367", ]

setGeneric("FilterMarkerByMajorGenotype",
           function(object, cutoff=NULL, print.n.markers=T) {
             standardGeneric("FilterMarkerByMajorGenotype")
             })

setMethod("FilterMarkerByMajorGenotype", "Genotype", function(object, cutoff=0.9, print.n.markers=T) {
  genotype_freq = AlleleCounting(object)@genotype_freq
  x = apply(genotype_freq, 1, max)/rowSums(genotype_freq)
  cutoff = round(cutoff*100)/100
  step_name = paste("Filtering out markers with more than", cutoff, "of the samples have the same call")
  SetMarkerFilter(object, Filter=x>cutoff, step_name=step_name, print.n.markers)
})
##---------------------------------------------------------
#' Filtering markers by the aboundance of the least aboundant genotype
#'
#' This function is used when you want to make sure all markers have the 3
#' possible genoytpe calls.
#'
#' @param object \code{Genotype}. The object to be filtered.
#' @param cutoff \code{numeric}. Markers with a less than \code{cutoff}
#' absolute number of calls in some genotype will be masked.
#' @param print.n.markers logical. Should the number of remaining markers be
#' printed on the screen?
#' @return a filtered \code{Genotype} object, with updated filter table and
#' \code{MarkerInfo$Filtered_at_step}
#' @examples
#' data("Mango_raw")
#' MarkerNo(MxS)
#' genotype_freq(AlleleCounting(MxS))["Contig1327_496", , drop=FALSE]
#' MxS = FilterMarkerByMinorGenotype(MxS, cutoff=2)
#' "Contig1327_496" %in% MarkerNames(MxS)
#' getMarkerInfo(MxS, Masked=F)["Contig1327_496", ]

setGeneric("FilterMarkerByMinorGenotype", function(object, cutoff=NULL, print.n.markers=T) {standardGeneric("FilterMarkerByMinorGenotype")})
setMethod("FilterMarkerByMinorGenotype", "Genotype", function(object, cutoff=2, print.n.markers=T) {
  genotype_freq = AlleleCounting(object)@genotype_freq
  x = apply(genotype_freq, 1, min)
  cutoff = round(cutoff*100)/100
  step_name = paste("Filtering out markers with less than", cutoff, "samples in each genotype")
  SetMarkerFilter(object, Filter=x<cutoff, step_name=step_name, print.n.markers=print.n.markers)
})

##---------------------------------------------------------
#' Filtering markers by some \code{MarkerInfo} column
#'
#' This functions allows you to filter markers by any property, such as the
#' Fluidigm team rank, or whether this marker is mapped to the reference genome.
#' To use it, the relevant property should appear in the object's
#' \code{MarkerInfo} table.
#'
#' @param object \code{Genotype}. The object to be filtered.
#' @param property.name \code{character}. The name of the \code{MarkerInfo}
#' column to use for filtering
#' @param property.criterion a single \code{character} or \code{numeric}
#' variable. If \code{MarkerInfo$property.name} is a \code{character}, use a
#' \code{character} \code{property.criterion} to filter all markers with
#' \code{property.name} equal to \code{property.criterion} or filter all
#' markers with \code{property.name} different from \code{property.criterion},
#' depending on \code{remove}. If \code{MarkerInfo$property.name} is a
#' \code{numeric} vector, use a \code{numeric} \code{property.criterion} as a
#' cutoff.
#' @param remove \code{character}. One of the following:
#'     "equal", "different", "greater", "greaterOrEquall", "smaller" or
#' "smallerOrEquall".
#' Which markers to remove with respect to \code{property.criterion}
#' @param print.n.markers logical. Should the number of remaining markers be
#' printed on the screen?
#' @return a filtered \code{Genotype} object, with updated filter table and
#' \code{MarkerInfo$Filtered_at_step}
#'
#' @examples
#' filter markers which Fluidigm team ranked >1
#' data("Mango_raw")
#' MarkerNo(MxS)
#' MxS = FilterMarkerByInfo(MxS, property.name="Category",
#' property.criterion=1, remove="greater")
#'
#' leave only markers from chromosome 1
#' data("citrus_clean")
#' MarkerNo(MxS)
#' MxS = FilterMarkerByInfo(MxS, property.name="Genome_Chromosome",
#' property.criterion="chr1", remove="different")

setGeneric("FilterMarkerByInfo",
           function(object, property.name="", property.criterion=NULL,
                    remove="greaterOrEquall", print.n.markers=T) {
             standardGeneric("FilterMarkerByInfo")
             })

setMethod("FilterMarkerByInfo", "Genotype", function(object, property.name="", property.criterion, remove="greaterOrEquall", print.n.markers=T) {
  MarkerInfo = getMarkerInfo(object)
  if (!(property.name %in% names(MarkerInfo)))
    stop(paste("No property called \'", property.name, "\' in marker info. Please provide a property name which appears in the marker info.", sep=""))
  x = switch(remove,
             equal = (MarkerInfo[[property.name]] == property.criterion),
             different = (MarkerInfo[[property.name]] != property.criterion),
             greater = (MarkerInfo[[property.name]] > property.criterion),
             greaterOrEquall = (MarkerInfo[[property.name]] >= property.criterion),
             smaller = (MarkerInfo[[property.name]] < property.criterion),
             smallerOrEquall = (MarkerInfo[[property.name]] <= property.criterion))
  step_name = "Filtering out SNPs with a "
  step_name = switch(property.name,
                     paste(step_name,  property.name, " ", sep=""),
                     Categories = paste(step_name,  "quality-category ", sep="")
  )
  step_name = switch(remove,
                     "greaterOrEquall" = paste(step_name,  ">=", sep=""),
                     "greater" = paste(step_name,  ">", sep=""),
                     "equall" = paste(step_name,  "=", sep=""),
                     "smallerOrEquall" = paste(step_name,  "<=", sep=""),
                     "smaller" = paste(step_name,  "<", sep=""),
                     "different" = paste(step_name,  "different than", sep="")
  )
  step_name = paste(step_name, property.criterion, sep="")
  SetMarkerFilter(object, Filter=x, step_name=step_name, print.n.markers=print.n.markers)

})

##---------------------------------------------------------
#' Filtering the samples
#'
#' Filter (or 'mask') a specified set of samples from a \code{Genotype} object.
#' Update the filter table and the sample info.
#'
#' @param object \code{Genotype}. The object to be filtered.
#' @param Filter character, integer or logical. The samples to be filtered.
#'  This could be a vector of sample names, sample indeces, or a logical vector
#'  with TRUE values for samples to filter and FALSE for samples to keep.
#' @param print.n.samples logical. Should the number of remaining samples be
#' printed on the screen?
#' @param step_name character. The name of the filtration step. This will
#' in the step name column in the filter table.
#' @return The filtered \code{Genotype} object. The filter table holds a new
#' line with the number of samples left after this step, and the
#' 'Filtered_at_step' column in the SampleInfo of the filtered samples, have
#' the step name.
#'
#' @examples
#' data("Mango_raw")
#' SampleNo(MxS)
#' MxS = SetSampleFilter(MxS, Filter=c("4_139", "77_13", "12"), step_name=
#' "arbitrary filtering")
#' getSampleInfo(MxS, Masked=F)["77_13", ]
setGeneric("SetSampleFilter",
           function(object,Filter=NULL,print.n.samples=T,
                    step_name="Filter samples") {
             standardGeneric("SetSampleFilter")
             })
setMethod("SetSampleFilter","Genotype",
          function(object,Filter=NULL, print.n.samples=T,
                   step_name="Filter samples") {
  cls = class(Filter)
  object@SampleFilter <- switch(cls,
                                NULL= stop("Filter is mandatory"),
                                integer = Filter[Filter%in%(1:SampleNo(object,Masked=FALSE))],
                                logical = object@SampleFilter[!Filter],
                                character = object@SampleFilter[!SampleNames(object)%in%Filter],
                                stop("Filter must be either characater, or logical")
  )
  if (print.n.samples) cat("Number of remaining samples: ",length(object@SampleFilter),"\n")
  SampleInfo = getSampleInfo(object, Masked=FALSE)
  if (!("Filtered_at_step" %in% names(SampleInfo))) {
    if (length(SampleInfo)) SampleInfo$Filtered_at_step = character(nrow(SampleInfo)) else
      SampleInfo = data.frame(Filtered_at_step = character(SampleNo(object)), stringsAsFactors=F)
  }
  f = setdiff(SampleNames(object, Masked=F), SampleNames(object, Masked=T))
  i = which(nchar(SampleInfo[f, "Filtered_at_step"])>0)
  if (length(i)) f = f[-i]
  SampleInfo[f, "Filtered_at_step"] = step_name
  object@SInfo = SampleInfo
  if (nrow(getFilterTable(object)) ==0) object = start.filter.table(object)
  UpdateFilterTable(object, step_name=step_name)
})

##---------------------------------------------------------
#' UnFiltering the samples
#'
#' \code{ResetSampleFilter} reset sample filter, so that all samples are
#'  unmasked
#'
#' @param object \code{Genotype}. The object to be un-filtered.
#' @return The unfiltered (with respect to samples) \code{Genotype} object.
#' @examples
#' data("Mango_raw")
#' SampleNo(MxS)
#' MxS = SetSampleFilter(MxS, Filter=c("Contig1987_1443", "Contig2006_1970",
#'     "mango_c15408_840"), step_name="arbitrary filtering")
#' getSampleInfo(MxS, Masked=F)["Contig1987_1443", ]
#' MxS = ResetSampleFilter(MxS)
setGeneric("ResetSampleFilter", function(object){
  standardGeneric("ResetSampleFilter")
  })
setMethod("ResetSampleFilter","Genotype", function(object){
  object@SampleFilter = 1:SampleNo(object,Masked=FALSE)
  cat("Number of samples: ",length(object@SampleFilter),"\n")
  UpdateFilterTable(object, step_name="Reset sample filter")
})
##---------------------------------------------------------
#' Getting the sample Filter
#'
#' @param object \code{Genotype}. The object to get the filter from.
#' @return The indeces of the unmasked samples.
setGeneric("GetSampleFilter", function(object){
  standardGeneric("GetSampleFilter")
  })
setMethod("GetSampleFilter","Genotype", function(object){
  object@SampleFilter
})
##---------------------------------------------------------
#' Quality filtering of samples
#'
#' A general function for several standard quality filtrations of samples
#'
#' @param object \code{Genotype}. The object to be filtered.
#' @param quality.criterion chatacter. The name of the quality criterion to
#' filter the samples by. Currently there is only one acceptable value:
#'    "no.calls" - filter samples with relative number of no calls > cutoff.
#' @param cutoff \code{numeric}. The cutoff for the filtration
#' @param print.n.samples logical. Should the number of remaining samples be
#' printed on the screen?
#' @return The filtered \code{Genotype} object, with updated filter table and
#' sample info table.
#'
#' @examples
#' data("Mango_raw")
#' SampleNo(MxS)
#' MxS = FilterSampleByQuality(MxS, quality.criterion="no.calls", cutoff=0.3)
setGeneric("FilterSampleByQuality",
           function(object, quality.criterion="", cutoff=NULL,
                    print.n.samples=T) {
             standardGeneric("FilterSampleByQuality")
             })
setMethod("FilterSampleByQuality", "Genotype", function(object, quality.criterion="", cutoff=NULL, print.n.samples=T) {
  return(switch(quality.criterion,
                no.calls = FilterSampleByNoCall(object, cutoff, print.n.samples=print.n.samples)
  ))
})
##---------------------------------------------------------
#' Filtering samples by relative number of no calls
#'
#' @param object \code{Genotype}. The object to be filtered
#' @param cutoff \code{numeric}. Samples with relative number of missing calls
#' higher than this value will be filtered.
#' @param print.n.markers logical. Should the number of remaining markers be
#' printed on the screen?
#' @return The same object with the samples with too many missing calls
#' filtered.
#'
#' @examples
#' data("Mango_raw")
#' MarkerNo(MxS)
#' MxS = FilterSampleByNoCall(MxS, cutoff=0.3)
setGeneric("FilterSampleByNoCall",
           function(object, cutoff=NULL,print.n.samples=T) {
             standardGeneric("FilterSampleByNoCall")
           })
setMethod("FilterSampleByNoCall", "Genotype",
          function(object, cutoff=1/3, print.n.samples=print.n.samples) {
            l = MarkerNo(object)
            Calls = getCalls(object)
            NoCalls = colSums(is.na(Calls))
            x = NoCalls/l > cutoff
            cutoff = round(cutoff*100)/100
            step_name = paste("Filtering out samples with >", cutoff, "no calls")
            SetSampleFilter(object, Filter=x, step_name=step_name, print.n.samples=print.n.samples)
          })

##---------------------------------------------------------
#' Filtering samples by some \code{SampleInfo} column
#'
#' This functions allows you to filter samples by any property. To use it, the
#'  relevant property should appear in the object's \code{SampleInfo} table.
#'
#' @param object \code{Genotype}. The object to be filtered.
#' @param property.name \code{character}. The name of the \code{SampleInfo}
#' column to use for filtering
#' @param property.criterion a single \code{character} or \code{numeric}
#' variable. If \code{SampleInfo$property.name} is a \code{character}, use a
#' \code{character} \code{property.criterion} to filter all samples with
#' \code{property.name} equal to \code{property.criterion} or filter all
#' samples with \code{property.name} different from \code{property.criterion},
#' depending on \code{remove}. If \code{SampleInfo$property.name} is a
#' \code{numeric} vector, use a \code{numeric} \code{property.criterion} as a
#' cutoff.
#' @param remove \code{character}. One of the following:
#'     "equal", "different", "greater", "greaterOrEquall", "smaller" or
#' "smallerOrEquall".
#' Which samples to remove with respect to \code{property.criterion}
#' @param print.n.samples logical. Should the number of remaining samples be
#' printed on the screen?
#' @return a filtered \code{Genotype} object, with updated filter table and
#' \code{SampleInfo$Filtered_at_step}
#'
#' @examples
#' Leaving only samples of graft strains
#' data("Mango_raw")
#' SampleNo(MxS)
#' MxS = FilterSampleByInfo(MxS, property.name="Rootstock/Graft",
#' property.criterion="G", remove="different")
#'

setGeneric("FilterSampleByInfo",
           function(object, property.name="", property.criterion=NULL,
                    remove="greaterOrEquall",print.n.samples=T) {
             standardGeneric("FilterSampleByInfo")
             })

setMethod("FilterSampleByInfo", "Genotype", function(object, property.name="", property.criterion, remove="greaterOrEquall", print.n.samples=T) {
  SampleInfo = getSampleInfo(object)
  if (!(property.name %in% names(SampleInfo)))
    stop(paste("No property called \'", property.name, "\' in sample info. Please provide a property name which appears in the sample info.", sep=""))
  x = switch(remove,
             equal = (SampleInfo[[property.name]] == property.criterion),
             different = (SampleInfo[[property.name]] != property.criterion),
             greater = (SampleInfo[[property.name]] > property.criterion),
             greaterOrEquall = (SampleInfo[[property.name]] >= property.criterion),
             smaller = (SampleInfo[[property.name]] < property.criterion),
             smallerOrEquall = (SampleInfo[[property.name]] <= property.criterion))

  step_name = paste("Filtering out Samples with a",  property.name)
  step_name = switch(remove,
                     "greaterOrEquall" = paste(step_name,  ">=", sep=""),
                     "greater" = paste(step_name,  ">", sep=""),
                     "equall" = paste(step_name,  "=", sep=""),
                     "smallerOrEquall" = paste(step_name,  "<=", sep=""),
                     "smaller" = paste(step_name,  "<", sep=""),
                     "different" = paste(step_name,  " different than", sep="")
  )
  step_name = paste(step_name, property.criterion)
  SetSampleFilter(object, Filter=x, step_name=step_name, print.n.samples=print.n.samples)
})
##---------------------------------------------------------
#' Eliminating linkage-disequilibrium redundancy among markers
#'
#' Leaving only one marker from each pair of linked markers. Using
#' \code{genetics} package function \code{LD} to calculate linkage
#' disequilibrium among marker.
#'
#' @param object \code{Genotype}. The object to be filtered
#' @param by \code{character}. The measure to determine if two markers are
#' linked. One of the followings:
#' \describe{
#'   \item{"R^2"}{Filter by r-square between markres}
#'   \item{"p" or "P"}{Filter by Chi-square p-value for marker independence}
#' }
#' @param cutoff \code{numberic} - The cutoff for filtering
#' @param score \code{numeric}. An optional score value for each marker, so
#' when choosing which marker to filter from a pair of linked markers, the one
#' with the higher score will be filtered, and the lower score will be kept.
#'
#' @return The same object with the filtered markers, so that all markers are
#' unlinked.
#'
#' @examples
#' data("Genotype_43X80")
#' MarkerNo(MxS)
#' MxS = FilterLDMarkers(MxS)
setGeneric("FilterLDMarkers",
           function(object, by="R^2", cutoff=NULL, score=NULL,
                    print.n.markers=T) {
             standardGeneric("FilterLDMarkers")
             })
setMethod("FilterLDMarkers","Genotype",function(object, by="R^2", cutoff=NULL, score=NULL, print.n.markers=T) {
  genotype_call = apply(getCalls(object), 1, genetics::genotype, sep=1)
  genotype_call = as.data.frame(genotype_call)
  for (i in 1:length(genotype_call))
    genotype_call[[i]] = genetics::as.genotype(genotype_call[[i]])
  LD_genotype = genetics::LD(genotype_call)
  if (by %in% c("p", "P")) by="P-value"
  if (is.null(cutoff)) {
    cutoff = switch(by,
                    "R^2" =  0.7,
                    "P-value" = 1e-14)
  }
  pairs = switch(by,
                 "R^2" = which(LD_genotype[[by]]>cutoff, arr.ind=T),
                 "P-value" = which(LD_genotype[[by]]<cutoff, arr.ind=T)
  )
  p = as.numeric(LD_genotype[["P-value"]])
  p = p[!is.na(p)]
  m = min(p[p>0])
  p = apply(pairs, 1, function(x) LD_genotype[["P-value"]][x[1], x[2]])
  P = -log10(p+m/2)
  ind = unique(c(pairs[, 1], pairs[, 2]))
  rP = tapply(P, factor(pairs[, 1], levels=ind), sum)
  rP[is.na(rP)] = 0
  cP = tapply(P, factor(pairs[, 2], levels=ind), sum)
  cP[is.na(cP)] = 0
  indP = rP+cP
  if (is.null(score)) score = rep(1, MarkerNo(object))
  o = order(score[ind], indP, decreasing=T)
  remove = integer()
  while (nrow(pairs)) {
    # i = as.numeric(names(which.max(indP)))
    i = ind[o[1]]
    pairsind = pairs[, 1]==i | pairs[, 2]==i
    x = c(pairs[, 1][pairsind], pairs[, 2][pairsind])
    P = P[!(pairs[, 1] %in% x | pairs[, 2] %in% x)]
    pairs = pairs[!(pairs[, 1] %in% x | pairs[, 2] %in% x), , drop=F]
    remove = c(remove, setdiff(x,i))
    ind = unique(c(pairs[, 1], pairs[, 2]))
    rP = tapply(P, factor(pairs[, 1], levels=ind), sum)
    rP[is.na(rP)] = 0
    cP = tapply(P, factor(pairs[, 2], levels=ind), sum)
    cP[is.na(cP)] = 0
    indP = rP+cP
    o = order(score[ind], indP, decreasing=T)
  }
  step_name = "Leaving only one marker from each pair of linked markers"
  step_name = switch(by,
                     "R^2" = paste(step_name,  " R^2 >", cutoff, ")", sep=""),
                     "P-value" = paste(step_name,  " P-value <", cutoff, ")", sep="")
  )
  SetMarkerFilter(object, Filter=MarkerNames(object)[remove], step_name=step_name, print.n.markers=print.n.markers)
})

##-----------------------------------------------------------------------------
#' Remove a specific filter from a a \code{Genotype} object.
#'
#' All samples or markers that were removed by this filtration step will be
#' un-masked. Use this carefully. If other filtration steps were taking place
#' after the step that is being removed, the un-masked markers or samples will
#' not be affected by them.
#'
#' @param object \code{Genotype}. The object to remove the filter from.
#' @param step_name character. The name of the filtration step to be cancelled.
#' @param print.n.markers logical. Should the new number of remaining markers
#' be printed on the screen?
#' @return The \code{Genotype} object with the newly unfiltered markers or
#' samples.
#'
#' @examples
#' data("citrus_clean")
#' MarkerNo(MxS)
#' MxS = RemoveFilter(MxS, "Filtering out SNPs with a Category >=3")
#' getFilterTable(MxS)
setGeneric("RemoveFilter",function(object,step_name="no_step",print.n=T) {
  standardGeneric("RemoveFilter")
  })
setMethod("RemoveFilter","Genotype",
          function(object, step_name="no_step", print.n=T) {
  filter.table = getFilterTable(object)
  if (!step_name %in% filter.table$Step) stop("step_name not in filter table.")
  MarkerInfo = getMarkerInfo(object, Masked=F)
  if (step_name %in% MarkerInfo$Filtered_at_step) {
    x = which(MarkerInfo$Filtered_at_step == step_name)
    if (sum(x %in% object@MarkerFilter)) stop ("Some filtered markers were not actually filtered.")
    object@MarkerFilter = union(object@MarkerFilter, x)
    if (!("Removed_Filter" %in% names(MarkerInfo))) MarkerInfo$Removed_Filter = rep("", nrow(MarkerInfo))
    MarkerInfo$Removed_Filter[x] = MarkerInfo$Filtered_at_step[x]
    MarkerInfo$Filtered_at_step[x] = rep("", length(x))
    object@MInfo = MarkerInfo
    if (print.n) cat("Number of remaining markers: ",length(object@MarkerFilter),"\n")
  }
  SampleInfo = getSampleInfo(object, Masked=F)
  if (step_name %in% SampleInfo$Filtered_at_step) {
    x = which(SampleInfo$Filtered_at_step == step_name)
    if (sum(x %in% object@SampleFilter)) stop ("Some filtered samples were not actually filtered.")
    object@SampleFilter = union(object@SampleFilter, x)
    if (!("Removed_Filter" %in% names(SampleInfo))) SampleInfo$Removed_Filter = rep("", nrow(SampleInfo))
    SampleInfo$Removed_Filter[x] = SampleInfo$Filtered_at_step[x]
    SampleInfo$Filtered_at_step[x] = rep("", length(x))
    object@SInfo = SampleInfo
    if (print.n) cat("Number of remaining samples: ",length(object@SampleFilter),"\n")
  }
  return(UpdateFilterTable(object, step_name=paste("Undo:", step_name)))
})

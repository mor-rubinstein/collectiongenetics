#' Creating a filter table for a \code{Genotype} object
#'
#' This is a function which usually runs automatically with the first call to
#' SetMarkerFilter or to SetSampleFilter
setGeneric("start.filter.table",function(object, anyway=F) {
  standardGeneric("start.filter.table")
  })
setMethod("start.filter.table","Genotype",function(object, anyway=F){
  object <- callNextMethod()
  object
})
setMethod("start.filter.table","filterTable",function(object, anyway=F){
  if ("filter_table" %in% slotNames(object))
    if (nrow(object@filter_table)>1)
      if (!anyway) {
        cat("object has already a non-empty filter table.\nTo start a filter table any way, call \'start.filter.table\' again with \'anyway=TRUE\'")
        return(object)
      }
  filter_table = data.frame(Step="begining", No_Markers=MarkerNo(object, Masked=F),
                            Markers_From_begining=1, Markers_From_last_step=1,
                            No_Samples=SampleNo(object, Masked=F),
                            Samples_From_begining=1, Samples_From_last_step=1)
  object@filter_table = filter_table
  object
})

##-----------------------------------------------------------------------------
#' Keeping the filter table updated
#'
#' the filter table keeps track of all markers and samples filtration steps.
#' each such step is followed by a call to \code{UpdateFilterTable}.
#'
#' @param object \code{Genotype}. The object which has the table you need to
#' update.
#' @param step_name character. The name of the filtration step.
#' @return The \code{Genotype} object with the new line in its filter table
#' specifying the step name and the number of markers/samples left after it.
#'
setGeneric("UpdateFilterTable",function(object, step_name) {
  standardGeneric("UpdateFilterTable")
  })
setMethod("UpdateFilterTable", "Genotype",function(object, step_name) {
  filter_table = getFilterTable(object)
  i = nrow(filter_table)
  if (step_name=="Reset sample filter" & i<2) return(object)
  if (step_name=="Reset marker filter" & i<2) return(object)
  filter_table =
    rbind(filter_table, data.frame(Step=step_name, No_Markers=MarkerNo(object),
                                   Markers_From_begining=MarkerNo(object)/filter_table$No_Markers[1],
                                   Markers_From_last_step=MarkerNo(object)/filter_table$No_Markers[i],
                                   No_Samples=SampleNo(object),
                                   Samples_From_begining=SampleNo(object)/filter_table$No_Samples[1],
                                   Samples_From_last_step=SampleNo(object)/filter_table$No_Samples[i]))
  object@filter_table = filter_table
  object
})

##-----------------------------------------------------------------------------
#' Getting the filter table
#'
#' @param object \code{Genotype} or \code{filterTable}. The object to get the
#' filter table from
#' @return The filter table

setGeneric("getFilterTable",function(object) {
  standardGeneric("getFilterTable")
  })
setMethod("getFilterTable","filterTable",function(object) {
  object@filter_table
})
setMethod("getFilterTable","Genotype",function(object) {
  callNextMethod()
})
##-----------------------------------------------------------------------------
#' Printing the filter table to a csv file.
#'
#' @param object \code{Genotype}. The object which has the table.
#' @param path \code{character}. The path to the folder to write the table in.
#' @param filenam \code{character}. the filename
#' @return \code{NULL}
#'
setGeneric("write.filter.table",
           function(object, path="./", filename="filter_table.csv") {
  standardGeneric("write.filter.table")
             })
setMethod("write.filter.table", "filterTable",
          function(object, path="./", filename="filter_table.csv") {
  filter_table <- getFilterTable(object)
  for (i in grep("From", names(filter_table)))
    filter_table[[i]] = round(100 * filter_table[[i]])
  names(filter_table)[grep("From", names(filter_table))] = paste("%", names(filter_table)[grep("From", names(filter_table))])
  write.table(filter_table, file=paste(path, filename, sep="/"), quote=F, sep="\t", row.names=F)
})
setMethod("write.filter.table", "Genotype",
          function(object, path="./", filename="filter_table.csv") {
            callNextMethod()
          })

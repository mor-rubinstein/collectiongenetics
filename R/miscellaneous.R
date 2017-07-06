##-----------------------------------------------------------------------------
#' Find a minimal set of markers to distinguish between all samples
#'
#' @param object \code{Genotype}. The object conatining the genotype calls to
#' use.
#' @param marker.score \code{numeric}. Optional. A vector of marker score. In
#' case there is a choice between several markers, a marker with a higher score
#' will be preferred.
#' @param filter.nocall \code{character}. Should markers with missing calls be
#' filtered? options are:
#' \describe{
#' \item{quantile}{filter out the 25\% of markers with highest rate of no call}
#' \item{all}{filter out the markers with any missing call}
#' \item{none}{no filtering.}
#' }
#'
#' @return \code{character}. The names of the markers in the minimal set.
#'
#' @examples
#' data("citrus_clean")
#' myset = MinimalSet(MxS)
#'
setGeneric("MinimalSet",function(object, marker.score=NULL, filter.nocall=T) {
  standardGeneric("MinimalSet")
  })
setMethod("MinimalSet","Genotype",
          function(object, marker.score=NULL, filter.nocall="quantile") {
  Calls = getCalls(object)
  if (is.null(marker.score))
    marker.score = rep(0, nrow(Calls))
  names(marker.score) = rownames(Calls)
  genotype_call = apply(Calls, 1, genetics::genotype, sep=1)
  genotype_call = as.data.frame(genotype_call)
  cat("Calculating LD between markers.\n")
  for (i in 1:length(genotype_call))
    genotype_call[[i]] = genetics::as.genotype(genotype_call[[i]])
  LD_genotype = genetics::LD(genotype_call)[["R^2"]]
  cat("Filtering markers with No Calls.\n")
  allele = c("XX", "YY", "XY")
  n_nocall = apply(Calls, 1, function(x) sum(!(x %in% allele)))
  q = quantile(n_nocall, 0.75)
  x = switch(filter.nocall,
             quantile = n_nocall<=q,
             all = n_nocall==0,
             none = rep(T, length(n_nocall))
  )
  if (filter.nocall=="all" && sum(x)==0) {
    filter.nocall = "quantile"
    warning("All markers had No Calls enteries.")
  }
  Calls = Calls[x, ]
  LD_genotype = LD_genotype[x, x]
  marker.score = marker.score[x]
  nocall.message = switch(filter.nocall,
                          quantile = paste("Removed", sum(!x), "Markers with more than", q, "No Calls each.\n"),
                          all = paste("Removed", sum(!x), "Markers with more than 0 No Calls each.\n"),
                          none = paste("Removed No Markers for having No Calls enteries.\n")
  )
  cat(nocall.message)
  ## Creating a matrix of where each row is a pair of samples and each column
  ## is a marker If the marker differentiates between the pair of samples, the
  ## corresponding entry will have '1'.
  A = character(SampleNo(object)*(SampleNo(object)-1)/2)
  i=0
  for(S1 in SampleNames(object)[-SampleNo(object)]){
    for (S2 in SampleNames(object)[(which(SampleNames(object)==S1)+1):SampleNo(object)]){
      i = i+1
      A[i] = paste(S1,"/",S2,sep="")
    }
  }
  M = matrix(0,nrow=length(A),ncol=nrow(Calls),
             dimnames=list(A, rownames(Calls)))
  for(S1 in SampleNames(object)[-SampleNo(object)]) {
    for (S2 in SampleNames(object)[(which(SampleNames(object)==S1)+1):SampleNo(object)]) {
      M[paste(S1,"/",S2,sep=""), ] = as.numeric(apply(Calls[,c(S1,S2)],1,function(x) x[1]!=x[2]))
    }
  }
  DiffNo = apply(M ,2 ,sum, na.rm = TRUE)
  Marker = character()
  cat("Finding a minimal set of markers:\n")
  while (max(DiffNo, na.rm = T)) {
    x = DiffNo > 0
    M = M[, x, drop=F]
    DiffNo = DiffNo[x]
    marker.score = marker.score[x]

    x = marker.score == max(marker.score)

    if (length(Marker)) {
      x = colnames(M)[x]
      LD_max = suppressWarnings(apply(LD_genotype[x, Marker, drop=F], 1, max, na.rm=T))
      LD_max = pmax(LD_max, suppressWarnings(apply(LD_genotype[Marker, x, drop=F], 2, max, na.rm=T)))
      marker =  x[which.min(LD_max)]
    } else marker = names(which.max(DiffNo[x]))
    Marker = c(Marker, marker)
    cat("number of additional differing pairs for marker No.", length(Marker), ": ", DiffNo[marker], "\n", sep="")
    x = M[, marker] == 0
    x[is.na(x)] = T
    M = M[x, ,drop=FALSE]
    DiffNo = apply(M ,2 ,sum, na.rm = TRUE)
  }
  cat("number of pairs without any differentiating markers: ", nrow(M), ".\nFinished.\n")
  return(Marker)
})

##-----------------------------------------------------------------------------
#' Markers with shared alleles between a pair of samples
#'
#' @param v1 \code{character}. a vector of genotype calls for one sample.
#' @param v2 \code{character}. a vector of genotype calls for the other sample.
#'
#' @return \code{logical}. A vector of the same length as the number of markers
#' with no missing calls in both samples, specifying for each marker whether
#' the two samples have a common allele.
#'
pair_common_alleles <- function(v1, v2) {
  na = is.na(v1) | is.na(v2)
  v1 = v1[!na]
  v2 = v2[!na]
  s = (v1=="XX" & v2 =="YY") | (v1=="YY" & v2 =="XX")
  return(!s)
}

##-----------------------------------------------------------------------------
#' Find potential selfings
#'
#' Finding potential pairs of a parent and an offspring, which is the result of
#' the parent selfing, within a \code{Genotype} object.
#'
#' @param object \code{Genotype}. The object to find selfings in.
#' @param error.freq \code{numeric}. A value between 0 and 1, with the
#' estimated genotype call error rate.
#'
#' @return \code{data.frame}. A table with potential selfings. Each row
#' represents a pair of selfing parent and its offpring. The columns are:
#' \describe{
#'   \item{selfing_parent}{The potential selfing parent}
#'   \item{offspring}{The potential offspring}
#'   \item{agreement}{The relative number of markers in agreement with this
#'    possibility}
#'   \item{NMarkers}{The total number of markers compared}
#'   \item{identity}{The relative number of markers identical in both samples}
#'   \item{parent_heterozygousity}{The relative number of heterozygous loci in
#'    the potential parent}
#'   \item{offspring_heterozygosity}{The relative number of heterozygous loci
#'   in the potential offspring}
#'   \item{generation}{The total number selfing generations, as calculated from
#'    the reduction in heterozygousity.}
#' }
#' @examples
#' data("citrus_clean")
#' self = FindSelfing(MxS, error.freq=0.05)

setGeneric("FindSelfing",function(object,error.freq=0) {
  standardGeneric("FindSelfing")
  })
setMethod("FindSelfing","Genotype", function(object, error.freq=0) {
  compare_pair<-function(v1, v2) {
    na = is.na(v1) | is.na(v2)
    v1 = v1[!na]
    v2 = v2[!na]
    v1 = v1[v2 != "XY"]
    v2 = v2[v2 != "XY"]
    v1 = strsplit(v1, split="")
    v2 = strsplit(v2, split="")
    s = sapply(1:length(v1), function(x) all(v1[[x]] %in% v2[[x]]))
    missmatch = sum(!s)/length(s)
    return(missmatch)
  }
  M = getCalls(object)
  colnames(M) = paste("sample", colnames(M), sep="_")
  Names = paste("sample", SampleNames(object), sep="_")
  candidate_pairs_table = data.frame(selfing_parent=character(), offspring=character(),
                                     agreement=numeric(), NMarkers=numeric(), identity=numeric(),
                                     parent_heterozygousity=numeric(), offspring_heterozygosity=numeric(),
                                     generation=numeric(), stringsAsFactors=F)
  for (S1 in Names) {
    V1 = M[, S1]
    for (S2 in setdiff(Names, S1)) {
      V2 = M[, S2]
      x = compare_pair(V1, V2)
      if (x <= error.freq)  {
        na = is.na(V1) | is.na(V2)
        v1 = V1[!na]
        v2 = V2[!na]
        id = sum(v1 == v2)/length(v2)
        Hv2 = (v2 == "XY")
        Hremain = sum(v1[Hv2] == "XY")/sum(Hv2)
        generation = log(Hremain, base=0.5)

        candidate_pairs_table = rbind(candidate_pairs_table, data.frame(selfing_parent=S2,
                                                                        offspring=S1, agreement=1-x, NMarkers=length(v1), identity=id,
                                                                        parent_heterozygousity=sum(V2=="XY", na.rm=T)/sum(!is.na(V2)),
                                                                        offspring_heterozygosity=sum(V1=="XY", na.rm=T)/sum(!is.na(V1)),
                                                                        generation=generation, stringsAsFactors=F))
      }
    }
  }
  candidate_pairs_table$selfing_parent = gsub("^sample_", "", candidate_pairs_table$selfing_parent)
  candidate_pairs_table$offspring = gsub("^sample_", "", candidate_pairs_table$offspring)
  candidate_pairs_table
})

##-----------------------------------------------------------------------------
#' Find potential parent-offspring pairs
#'
#' Finding potential pairs of a parent and an offspring, which is the result of
#' the parent selfing, within a \code{Genotype} object.
#'
#' @param object \code{Genotype}. The object to find selfings in.
#' @param error.freq \code{numeric}. A value between 0 and 1, with the
#' estimated genotype call error rate.
#'
#' @return \code{data.frame}. A table with each row representing a potential
#' pair of parent and offpring. The columns are:
#' \describe{
#'   \item{parent}{The potential parent}
#'   \item{offspring}{The potential offspring}
#'   \item{agreement}{The relative number of markers in agreement with this
#'    possibility}
#'   \item{NMarkers}{The total number of markers compared}
#'   \item{identity}{The relative number of markers identical in both samples}
#'   \item{expected_identity}{The expected relative number of markers identical
#'   between the parent and a potential offspring.}
#'}
#' @examples
#' data("citrus_clean")
#' PO = FindParentOffspring(MxS, error.freq=0.05)

setGeneric("FindParentOffspring", function(object, error.freq=0, MaskedSample=T) {
  standardGeneric("FindParentOffspring")
  })
setMethod("FindParentOffspring","Genotype",function(object, error.freq=0, MaskedSample=TRUE){
  M = getCalls(object, MaskedSample=MaskedSample)
  colnames(M) = paste("sample", colnames(M), sep="_")
  Names = paste("sample", SampleNames(object, Masked=MaskedSample), sep="_")
  candidate_pairs_table = data.frame(parent=character(), offspring=character(),
                                     agreement=numeric(), NMarkers=numeric(),
                                     identity=numeric(), expected_identity=numeric(),
                                     stringsAsFactors=F)
  for (S1 in 1:(SampleNo(object, Masked=MaskedSample)-1)) {
    V1 = M[, S1]
    for (S2 in (S1+1):SampleNo(object, Masked=MaskedSample)) {
      V2 = M[, S2]
      s = pair_common_alleles(V1, V2)
      if (length(s)) {
        x = sum(!s)/length(s)
        if (x <= error.freq) {
          na = is.na(V1) | is.na(V2)
          v1 = V1[!na]
          v2 = V2[!na]
          id = sum(v1 == v2)/length(v2)
          expected_id1 = 1-(sum(v1=="XY")/(2*length(v1)))
          expected_id2 = 1-(sum(v2=="XY")/(2*length(v2)))
          candidate_pairs_table = rbind(candidate_pairs_table,
                                        data.frame(parent=Names[c(S1, S2)], offspring=Names[c(S2, S1)],
                                                   agreement=1-x, NMarkers=length(v1), identity=id,
                                                   expected_identity=c(expected_id1, expected_id2), stringsAsFactors=F))
        }
      }
    }
  }
  candidate_pairs_table$parent = gsub("^sample_", "", candidate_pairs_table$parent)
  candidate_pairs_table$offspring = gsub("^sample_", "", candidate_pairs_table$offspring)
  candidate_pairs_table
})

##-----------------------------------------------------------------------------
#' Find the match of a potential parents-offspring trio
#'
#' @param voffspring \code{character}. The genotype calls of the offspring, in
#' a "XX", "XY", "YY" format.
#' @param vparent1 \code{character}. The genotype calls of the first parent, in
#' a "XX", "XY", "YY" format.
#' @param vparent1 \code{character}. The genotype calls of the second parent,
#' in a "XX", "XY", "YY" format.
#'
#' @return \code{logical}. A vector specifying for each marker, if it supports
#' the possibility of the offspring being the offspring of the two parents.
#'
trio_genotype_match <- function(voffspring, vparent1, vparent2) {
  na = is.na(voffspring) | is.na(vparent1) | is.na(vparent2)
  v = voffspring[!na]
  v1 = vparent1[!na]
  v2 = vparent2[!na]
  match_trio_table = build_match_trio_table()
  match_trio_table = paste(match_trio_table$parent1, match_trio_table$parent2,
                           match_trio_table$offspring, sep="_")
  trio_match = paste(v1, v2, v, sep="_") %in% match_trio_table
  return(trio_match)
}

# An auxiliary function used by 'trio_genotype_match'
build_match_trio_table =  function() {
  match_trio_table = data.frame(parent1=rep(c("XX", "XY", "YY"), each=12),
                                parent2=rep(rep(c("XX", "XY", "YY"), each=4), 3),
                                stringsAsFactors=F)
  parent1 = strsplit(match_trio_table$parent1, split="")
  parent2 = strsplit(match_trio_table$parent2, split="")
  s1 = rep(1:2, 18)
  s2 = rep(c(1,1,2,2), 9)
  match_trio_table$offspring = sapply(1:nrow(match_trio_table), function(x) {
    a1 = parent1[[x]][s1[x]]
    a2 = parent2[[x]][s2[x]]
    paste(a1, a2, sep="")
  })
  match_trio_table$offspring[match_trio_table$offspring=="YX"] = "XY"
  x = paste(match_trio_table$parent1, match_trio_table$parent2, match_trio_table$offspring, sep="_")
  match_trio_table = match_trio_table[!duplicated(x), ]
  rm(parent1, parent2, s1, s2, x)
  rownames(match_trio_table) = NULL
  match_trio_table
}

#' A \code{data.frame} for a list of pairs of genetically identical samples
#'
#' A data frame with the following variables:
#'    sample1: the first sample in the identical pair
#'    sample2: the second sample in the identical pair
#'    identity: the relative number of markers in which the pair is identical.
#'
#' @seealso \code{\link{FindIdentical}} for finding pairs of identical samples,
#' \code{\link{IdenticalGroups}} for identifying groups of identical samples,
#'  \code{\link{AddGroup2SampleInfo}} for adding the identity group to the
#'  sample information, and \code{\link{Unique}} for filtering all samples but
#'  one from each identity group.

setClass(Class="IdenticalSamples",
         representation("data.frame"),
         prototype(data.frame(sample1=character(0), sample2=character(0),
                              identity=numeric(), stringsAsFactors=FALSE))
)

##-----------------------------------------------------------------------------

#' Find pairs of identical samples
#'
#' finding pairs of individuals which are genetically identical accross the
#' panel of marker set
#'
#' @param object \code{Genotype}. The name of the object to find the identical
#' samples in.
#' @param error.freq \code{numeric}. A single numeric value: the estimated
#' genotype call error frequency.
#' @return \code{IdenticalSamples}. A table of identical pairs.
#'
#' @details \code{FindIdentical} considers any pair of samples, that have
#' identical genotype call in \code{1-error.freq} or more markers, as
#' genetically identical.
#'
#' @examples
#' data("Mango_row")
#' F0 = FindIdentical(MxS)
#' F0.05 = FindIdentical(MxS, error.freq=0.05)
#' @seealso \code{\link{IdenticalGroups}} for creating groups of identical
#' samples, \code{\link{AddGroup2SampleInfo}} for adding the identity group to
#' the sample information, and \code{\link{Unique}} for filtering all samples
#' but one from each identity group.

setGeneric("FindIdentical",function(object,error.freq=0) {
  standardGeneric("FindIdentical")
  })
setMethod("FindIdentical","Genotype",function(object, error.freq=0) {
  compare_pair<-function(v1, v2) {
    na = is.na(v1) | is.na(v2)
    v1 = v1[!na]
    v2 = v2[!na]
    return(sum(v1 != v2)/length(v1))
  }
  M = getCalls(object)
  colnames(M) = paste("sample", colnames(M), sep="_")
  Names = paste("sample", SampleNames(object), sep="_")
  sample1 <- sample2 <- character()
  identity <- numeric()
  for (S1 in 1:(SampleNo(object)-1)) {
    V1 = M[, S1]
    for (S2 in (S1+1):SampleNo(object)) {
      V2 = M[, S2]
      cp = compare_pair(V1, V2)
      if (!is.na(cp) && cp <= error.freq) {
        sample1 = c(sample1, Names[S1])
        sample2 = c(sample2, Names[S2])
        identity = c(identity, 1-cp)
      }
    }
  }
  a = data.frame(sample1 = sample1, sample2 = sample2, identity = identity,
                 stringsAsFactors = F)
  identical_pairs_table = new("IdenticalSamples", a)
  identical_pairs_table$sample1 = gsub("^sample_", "", identical_pairs_table$sample1)
  identical_pairs_table$sample2 = gsub("^sample_", "", identical_pairs_table$sample2)
  new("IdenticalSamples", identical_pairs_table)
})
##-------------------------------------------------------------------------------

#' Grouping together identical samples
#'
#' @param id_table \code{IdenticalSamples}. A special kind of data.frame which
#' specifies pairs of identical samples
#' @return \code{list}. A list of \code{character} vectors. Each one is a set
#' of sample names which are all genetically identical
#'
#' @details \code{IdenticalGroups} uses the \code{MCL} package \code{mcl}
#' function to partition the samples into groups.
#'
#' @seealso \code{\link{FindIdentical}} for finding pairs of identical samples,
#' , \code{\link{AddGroup2SampleInfo}} for adding the identity group to the
#' sample information and \code{\link{Unique}} for filtering all samples but
#' one from each identity group.
#'
#' @examples
#' data("Mango_row")
#' F0.05 = FindIdentical(MxS, error.freq=0.05)
#' group = IdenticalGroups(F0.05)

setGeneric("IdenticalGroups",function(id_table) {
  standardGeneric("IdenticalGroups")
  })
setMethod("IdenticalGroups", "IdenticalSamples", function(id_table) {
  s = unique(c(id_table$sample1, id_table$sample2))
  M = matrix(0, nrow=length(s), ncol=length(s), dimnames=list(s, s))
  diag(M) = 1
  for (i in 1:nrow(id_table)) {
    M[id_table$sample1[i], id_table$sample2[i]] = id_table$identity[i]
    M[id_table$sample2[i], id_table$sample1[i]] = id_table$identity[i]
  }
  x = try(MCL::mcl(M, addLoops=T, allow1=T), silent=T)
  if (class(x) == "try-error" |
      length(grep("An Error occurred at iteration", x))>0) {
    x = colSums(M)
    y = rep(1, ncol(M)) %*% t(x)
    M1 <- M2 <- M3 <-M/y
    Ma = M1 %*% M1
    Mb = Ma ^2
    x = colSums(Mb)
    y = rep(1, ncol(Mb)) %*% t(x)
    Mb <-Mb/y
    i=1
    while (max(abs(Mb-M1))>0 & max(abs(Mb-M2))>0 & max(abs(Mb-M3))>0 & i<=100) {
      # cat(i, "")
      M3 = M2
      M2 = M1
      M1 = Mb
      Ma = M1 %*% M1
      Mb = Ma ^2
      x = colSums(Mb)
      y = rep(1, ncol(Mb)) %*% t(x)
      Mb <-Mb/y
      i = i+1
    }
    id_list =  apply(Mb, 1, function(x) which(x>0))
    if (!all(1:nrow(M) %in% unlist(id_list))) stop("Grouping failed")
    i=1; while(i<length(id_list)) {
      j=i+1; while(j<=length(id_list)) {
        if (setequal(id_list[[i]], id_list[[j]])) id_list = id_list[-j] else j=j+1
      }
      i = i+1
    }
    if (sum(duplicated(unlist(id_list)))) stop("Grouping failed")
    id_list = id_list[sapply(id_list, length)>0]
    id_list = lapply(id_list, function(x) s[x])

    } else id_list = tapply(s, x$Cluster, function(x) x)
  return(id_list)
})

##-------------------------------------------------------------------------------
#' Adding information about group of identical samples to a Genotype object
#'
#' After applying the functions \code{FindIdentical} and \code{IdenticalGroups},
#' this function allows you to add an "Identity_group" column to the object's
#' \code{SampleInfo}, which specifies for samples that are a part of genetically
#' identical group, the group they belong to.
#'
#' @param id_table \code{Genotype}. The object to add the information to.
#' @param \code{list}. A list of \code{character} vectors. Each one is a set
#' of sample names which are all genetically identical
#' @return \code{Genotype}. The object with the sample information containing
#' the identity group.
#'
#' @seealso \code{\link{FindIdentical}} for finding pairs of identical samples,
#' \code{\link{IdenticalGroups}} for creating groups of identical samples, and
#' \code{\link{Unique}} for filtering all samples but one from each identity
#' group.
#'
#' @examples
#' data("Mango_row")
#' F0.05 = FindIdentical(MxS, error.freq=0.05)
#' group = IdenticalGroups(F0.05)
#' names(group) = paste("Group", names(group))
#' MxS = AddGroup2SampleInfo(MxS, group)
#' head(getSampleInfo(MxS))
#' getSampleInfo(MxS)[group[[7]], ]

setGeneric("AddGroup2SampleInfo",function(object, group=NULL) {
  standardGeneric("AddGroup2SampleInfo")
  })
setMethod("AddGroup2SampleInfo", "Genotype", function(object, group=NULL) {
  if (mode(group) != "list") stop("group must be a list of groups of samples with identical genotype")
  if (!all(unlist(group) %in% SampleNames(object, Masked=F))) stop("Not all the sample names in 'group' appear in the object sample names.")
  SampleInfo = getSampleInfo(object, Masked=F)
  SampleInfo$Identity_group = character(nrow(SampleInfo))
  for (i in 1:length(group)) SampleInfo[group[[i]], "Identity_group"] = rep(names(group)[i], length(group[[i]]))
  object@SInfo = SampleInfo
  return(object)
})

##-------------------------------------------------------------------------------
#' Remove redundency due to genetically identical samples
#'
#' This function filters the samples in a \code{Genotype} object, such that
#' from each group of genetically identical samples, only one representative
#' remains, and all the others are masked.
#'
#' @param object \code{Genotype}. The name of the object to remove the
#' redundency from.
#' @param error.freq \code{numeric}. A single numeric value: the estimated
#' genotype call error frequency.
#' @return \code{Genotype}. The filtered object.
#'
#' @details Any pair of samples, that have identical genotype call in
#' \code{1-error.freq} or more markers, are considered genetically identical.
#' \code{Unique} calls \code{FindIdentical} and \code{IdenticalGroups} to find
#' groups of identical samples. It chooses a sample with minimal number of
#' No-Calls as the of representative to keep from each group.
#'     After the filtration, the object slots 'SampleFilter' and 'filter_table'
#' are updated appropriately.
#'
#' @seealso \code{\link{FindIdentical}} for finding pairs of identical samples,
#' \code{\link{IdenticalGroups}} for creating groups of identical samples and
#' \code{\link{AddGroup2SampleInfo}} for adding the identity group to the
#' sample information.
#'
#' @examples
#' data("Mango_row")
#' SampleNo(MxS)
#' MxS = Unique(MxS, error.freq=0.05)
setGeneric("Unique",function(object,error.freq=0) {standardGeneric("Unique")})
setMethod("Unique","Genotype",function(object, error.freq=0) {
  remove = character()
  Identical_table = FindIdentical(object, error.freq)
  if (nrow(Identical_table) ==0) return(object)
  group = IdenticalGroups(Identical_table)
  M = getCalls(MxS)
  alleles = c("XX","XY","YY")
  for (i in 1:length(group)) {
    groupi = group[[i]]
    m = M[, groupi]
    na = colSums(is.na(m))
    groupi = groupi[na==min(na)]
    if (length(groupi) > 1) {
      m =m[, na==min(na)]
      y = t(apply(m, 2, function(x) table(x)[alleles]))
      colnames(y)= alleles
      y[is.na(y)]=0
      groupi = groupi[which.max(apply(y, 1, min))]
    }
    remove = c(remove, setdiff(group[[i]], groupi))
  }
  step_name = paste("Leaving only one sample from group of identical samples. (identity >=", 1-error.freq, ")", sep="")
  SetSampleFilter(object,Filter=remove, step_name=step_name)
})

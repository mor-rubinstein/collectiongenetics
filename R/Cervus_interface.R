
#' Write an input file with genotype calls for 'Cervus' software
#'
#' Cervus is a computer program for assignment of parents to their offspring
#' using genetic markers. It first calculates allele frequencies from a
#' genotype call file. Than it runs a simulation to determine the level of log
#' likelihood and delta required for parental assignment. Finally it assigns
#' uses a list of offsprings, and a list of candidate parents to each one of
#' them, to assign parents to offsprings, both when one of the parents is
#' known, and without a known parent.
#'
#' @param object \code{Genotype}. The object to use as input for Cervus. Must
#' be in "ACGT" format rather than in a "XY" format.
#' Masked samples and markers will not be included in this input.
#' @param fn \code{character}. The file name for the Cervus input file.
#' @param type \code{character}. The format (nuclotides or numbers) to use for
#' the Cervus input file.
#'
#' @return \code{NULL}.
#'
#' @seealso \code{\link{xy2atgc}}, \code{\link{write.candidates}},
#' \code{\link{read.cervus}}
#'
#' @examples
#' data("citrus_clean")
#' Alleles = getMarkerInfo(MxS, Masked=FALSE)$Alleles
#' names(Alleles) <- MarkerNames(MxS, Masked=F)
#' nucMxS=xy2atgc(MxS, alleles=Alleles)
#' write.cervus(nucMxS, fn="citrus_genotype.txt", type="numbers")

setGeneric("write.cervus",
           function(object, fn, type=c("nuc","numbers")) {
             standardGeneric("write.cervus")
             })
setMethod("write.cervus","Genotype",
          function(object, fn, type=c("nuc","numbers")) {
  if (!is.character(fn)) stop("file name is required!\n")
  MxS = getCalls(object)
  missing = unique(as.vector(MxS)[!grepl("[ATCG][ATCG]",as.vector(MxS))])
  MxS[MxS%in%missing] = "0"
  MxS1 = apply(MxS,2,function(x) stringr::str_replace(x,"^[ATGC]",""))
  MxS2 = apply(MxS,2,function(x) stringr::str_replace(x,"[ATGC]$",""))
  rownames(MxS1) <- rownames(MxS2) <-rownames(MxS)
  SxM_cervus = cbind(t(MxS1),t(MxS2))
  SxM_cervus = SxM_cervus[,order(colnames(SxM_cervus))]
  colnames(SxM_cervus) = paste(colnames(SxM_cervus),rep(c("a","b"),nrow(MxS)),sep="")
  M = apply(SxM_cervus,2,function(x) chartr("ACGT","1234",x))
  write.table(switch(type,nuc=SxM_cervus,numbers=M),
              file=fn,
              sep="\t",
              row.names=TRUE,
              col.names=NA, quote=F
  )
})

##-----------------------------------------------------------------------------

#' Write an input file with candidate parents for each offspring, for 'Cervus'
#' software
#'
#' Cervus is a computer program for assignment of parents to their offspring
#' using genetic markers. It first calculates allele frequencies from a
#' genotype call file. Than it runs a simulation to determine the level of log
#' likelihood and delta required for parental assignment. Finally it assigns
#' uses a list of offsprings, and a list of candidate parents to each one of
#' them, to assign parents to offsprings, both when one of the parents is
#' known, and without a known parent.
#'
#' @param object \code{Genotype}. The object to take the samples (offsprings
#' and candidate parents) from. Masked samples and markers will not be
#' included.
#' @param filename \code{character}. The file name for the candidate parents
#' file.
#' @param path \code{character}. The path for the candidate parents file.
#' @param identical_group \code{list}. A list of \code{character} vectors. Each
#' one is a set of sample names which are all genetically identical. Can be
#' created using the function \code{\link{IdenticalGroups}}.
#' @param ID_identical \code{IdenticalSamples}. A table of identical sample
#' pairs. Each row corresponds to one identical pair. The table must contain
#' columns called sample1 and sample2 with pair sample names. It can be created
#' using the function \code{\link{FindIdentical}}. No need to use this if
#' \code{identical_group} is present.
#' @param filter_by_generation \code{logical}. Should an information about the
#' generations of samples be used to eliminate candidate parents? If TRUE,
#' samples of younger generations cannot be candidate parents for sample from
#' older generations. To use this option, there must be a column of generation
#' information in the object's \code{SampleInfo} slot.
#' @param mm_cutoff \code{numeric}. The mis-match cutoff. A value between 0 and
#' 1. The relative number of markers with no common allele that is allowed for
#' a candidate parent-offspring pair.
#'
#' @return \code{NULL}.
#'
#' @examples
#' data("citrus_clean")
#' data("cit_id_groups")
#' write.candidates(MxS, identical_group=group, filter_by_generation=1,
#' mm_cutoff=0.95)
#'
#' @seealso \code{\link{write.cervus}}, \code{\link{read.cervus}}
#'
setGeneric("write.candidates", function(object, ...) {
  standardGeneric("write.candidates")
  })
setMethod("write.candidates","Genotype",
          function(object, ID_identical=NULL, identical_group=NULL,
                   filter_by_generation=F, path="Cervus",
                   filename="CandidateParent.txt", mm_cutoff=1) {
  ### Filtering out candidate parents that are identical to other candidate parents
  SampleInfo = getSampleInfo(object)
  if (filter_by_generation) {
    generation_col = grep("generation", names(SampleInfo), ignore.case=T)
    if (length(generation_col) >1) stop("Sample Info has too many columns with \'Generation\' in their names.")
    if (length(generation_col) < 1) {
      warning("Sample Info has no \'Generation\' column name.")
      filter_by_generation = F
    }
  }
  if (length(identical_group)) {
    Names_identical = unique(unlist(identical_group))
  } else {
    if (length(ID_identical)) {
      Names_identical = unique(c(ID_identical$sample1, ID_identical$sample2))
      Names_identical = Names_identical[Names_identical %in% SampleNames(object)]
    } else Names_identical = character()
  }
  if (filter_by_generation) {
    CandidateParrentID = lapply(rownames(SampleInfo), function(os) {
      candidate = rownames(SampleInfo)
      candidate = setdiff(candidate, os)
      ### filtering out samples of newer generations:
      osgeneration = SampleInfo[os, generation_col]
      candidate_generation = SampleInfo[candidate, generation_col]
      candidate_generation[is.na(candidate_generation)] = min(candidate_generation, na.rm=T)
      if (is.na(osgeneration)) osgeneration = max(candidate_generation,na.rm=T)
      candidate = candidate[candidate_generation <= osgeneration]
      ### filtering out samples that are identical to the offspring:
      if (length(identical_group)) {
        ind = which(sapply(identical_group, function(x) os %in% x))
        if (length(ind) >1) stop(paste(os, "appears in more than 1 identity groups."))
        if (length(ind)) {
          s = setdiff(identical_group[[ind]], os)
          candidate = setdiff(candidate, s)
        }
      } else {
        if (length(ID_identical)) {
          s1 = ID_identical$sample1 %in% os
          candidate = setdiff(candidate, ID_identical$sample2[s1])
          s2 = ID_identical$sample2 %in% os
          candidate = setdiff(candidate, ID_identical$sample1[s2])
        }
      }
      Tokeep = rep(T, length(candidate))
      names(Tokeep) = candidate
      Tokeep[names(Tokeep) %in% Names_identical] <- F
      Identical = candidate[candidate %in% Names_identical]
      if (length(identical_group)) {
        while (length(Identical)) {
          Name = Identical[1]
          Tokeep[Name] = T
          ind = which(sapply(identical_group, function(x) Name %in% x))
          if (length(ind) >1) stop(paste(Name, "appears in more than 1 identity groups."))
          if (length(ind)) Identical = setdiff(Identical, identical_group[[ind]])
        }
      } else {
        if (length(ID_identical)) {
          Calls = getCalls(object)[, Identical]
          Identical = Identical[order(colSums(is.na(Calls)))]
          while (length(Identical)) {
            Name = Identical[1]
            Tokeep[Name] = T
            x = ID_identical$sample2[ID_identical$sample1 == Name]
            x = c(x, ID_identical$sample1[ID_identical$sample2 == Name])
            Identical = setdiff(Identical, c(x, Name))
          }
        }
      }
      candidate = candidate[Tokeep]
      return(candidate)
    })
  } else {
    Tokeep = rep(T, nrow(SampleInfo))
    names(Tokeep) = rownames(SampleInfo)
    Tokeep[names(Tokeep) %in% Names_identical] <- F
    if (length(identical_group)) {
      while(length(Names_identical)) {
        Name = Names_identical[1]
        Tokeep[Name] = T
        ind = which(sapply(identical_group, function(x) Name %in% x))
        if (length(ind) >1) stop(paste(Name, "appears in more than 1 identity groups."))
        if (length(ind)) Names_identical = setdiff(Names_identical, identical_group[[ind]])
      }
    } else {
      if (length(ID_identical)) {
        Calls = getCalls(object)[, Names_identical]
        Names_identical = Names_identical[order(colSums(is.na(Calls)))]
        while(length(Names_identical)) {
          Name = Names_identical[1]
          Tokeep[Name] = T
          x = ID_identical$sample2[ID_identical$sample1 == Name]
          x = c(x, ID_identical$sample1[ID_identical$sample2 == Name])
          Names_identical = setdiff(Names_identical, c(x, Name))
        }
      }
    }
    CandidateParrentID = lapply(rownames(SampleInfo), function(x) {
      candidate = names(Tokeep)[Tokeep]
      candidate = setdiff(candidate, x)
      return(candidate)
    })
  }
  if (mm_cutoff < 1) {
    p_o = FindParentOffspring(object, error.freq=mm_cutoff, MaskedSample=F)
    for (i in 1:length(CandidateParrentID)) {
      CandidateParrentID[[i]] = intersect(CandidateParrentID[[i]], p_o$parent[p_o$offspring==rownames(SampleInfo)[i]])
    }
  }
  if (!dir.exists(path)) dir.create(path)
  con = file(paste(path, filename, sep="/"), open="w")
  for (i in 1:length(CandidateParrentID)) {
    x = paste(CandidateParrentID[[i]], collapse="\t")
    x = paste(rownames(SampleInfo)[i], x, sep="\t")
    writeLines(x, con=con)
  }
  close(con)
})

##-----------------------------------------------------------------------------
#' Read 'Cervus' results
#'
#' Cervus is a computer program for assignment of parents to their offspring
#' using genetic markers. It first calculates allele frequencies from a
#' genotype call file. Than it runs a simulation to determine the level of log
#' likelihood and delta required for parental assignment. Finally it assigns
#' uses a list of offsprings, and a list of candidate parents to each one of
#' them, to assign parents to offsprings, both when one of the parents is
#' known, and without a known parent.
#' Cervus outputs a .csv file where every row specifies one offspring-mother-
#' father trio, and includes the number of mismatching markers, log likelihood
#' and delta
#'
#' @param file_path \code{character}. The path and file name of the citrus .csv
#' file.
#'
#' @return \code{data.frame}. A table with the cervus results.
#'
#' @seealso \code{\link{write.cervus}}, \code{\link{write.candidates}}
#'
setGeneric("read.cervus", function(file_path) {standardGeneric("read.cervus")})
setMethod("read.cervus", "character", function(file_path) {
  ### Cervus added a ',' at the end of each line, which makes it unsuitable to read with 'read.delim' sort of functions, so I'm using 'scan' instead.
  x = scan(file_path, what="", sep="\n")
  A = strsplit(x, split=",")
  l = max(sapply(A, length))
  A[sapply(A, length)<l] = lapply(A[sapply(A, length)<l], function(a) {
    a = c(a, rep(NA, l-length(a)))
    a
  })
  if (substr(x[1], nchar(x[1]), nchar(x[1])) ==",")
    A = lapply(A, function(a) a[-l])
  A = t(as.data.frame(A, stringsAsFactors=F))
  cervus = as.data.frame(A[-1, ], stringsAsFactors=F)
  names(cervus) = gsub(" +", "_", A[1, ])
  rownames(cervus)=NULL
  for (i in grep("(ID|confidence)", names(cervus), invert=T)) cervus[[i]] = as.numeric(cervus[[i]])
  while(sum(duplicated(names(cervus)))) {
    x = which(duplicated(names(cervus)))
    names(cervus)[x] = paste(names(cervus)[x], ".1", sep="")
  }
  return(cervus)
})

##-----------------------------------------------------------------------------
#' Filtering 'Cervus' results
#'
#' Cervus is a computer program for assignment of parents to their offspring
#' using genetic markers. This function allows you to filter its results and
#' leave only most probable candidate parents.
#'
#' @param mm_cutoff \code{numeric}. The mis-match cutoff. The number of non-
#' mathcing markers that is allowed for a candidate parents-offspring trio.
#' @param lod_cutoff \code{numeric}. The lod score cutoff. A single numeric
#' value. The minimum lod score accepted for a candidate parents-offspring
#' trio. This value can be obtained from a simulation run of Cervus.
#'
#' @return \code{NULL}.
#'
#' @examples
#' data("cervus")
#' cervus = filter.cervus.results(cervus, mm_cutoff=9, lod_cutoff=9.84)
#'
setGeneric("filter.cervus.results",
           function(cervus, mm_cutoff=Inf, lod_cutoff=0) {
             standardGeneric("filter.cervus.results")
             })
setMethod("filter.cervus.results", "data.frame",
          function(cervus, mm_cutoff=Inf, lod_cutoff=0) {
  mm_column = grep("mismatching", names(cervus))
  mm_column = mm_column[length(mm_column)]
  lod_column = grep("LOD", names(cervus))
  lod_column = lod_column[length(lod_column)]
  check_parents = function(id) {
    ind_parents = which(cervus$Offspring_ID ==id)
    mm = cervus[[mm_column]][ind_parents]
    if (!is.na(mm[1])) {
      ind_parents = ind_parents[mm <= mm_cutoff]
        lod = cervus[[lod_column]][ind_parents]
      ind_parents = ind_parents[lod >= lod_cutoff]
    }
    if (length(ind_parents)) return(ind_parents) else return(integer())
  }
  parents_test = lapply(unique(cervus$Offspring_ID), check_parents)

  x = sapply(parents_test, length)==0
  best_candidate = cervus[unlist(parents_test), ]
  A = cervus[1:sum(x), ]
  A$Offspring_ID = unique(cervus$Offspring_ID)[x]
  off_col = names(A)[grep("Offspring", names(A))]
  for (i in off_col[-1]) {
    A[[i]] = sapply(A$Offspring_ID, function(x) {
      ind = which(cervus$Offspring_ID==x)
      ind = ind[1]
      cervus[[i]][ind]
    })
  }
  for (i in setdiff(names(A), off_col)) A[[i]] = rep(NA, nrow(A))
  return(best_candidate)
})

##-----------------------------------------------------------------------------
#' Using known parent-offspring pairs to filter non-matching markers
#'
#' If we have previous knowledge of parent-offspring pairs in our data (usually
#' mother-offspring) we can use this knowledge to filter out markers that don't
#' support it. A parent offspring pair must have a common allele. If a marker
#' does not match in a significant number of known and well established parent-
#' offspring pair, this marker may not be trusted. Provided that you have a
#' large enough number of previously known parent-offspring pairs, you can run
#' 'Cervus' software with known parent on your genotype data, and then use this
#'  function to get the match of each marker to these pairs.
#'
#' @param object \code{data.frame}. A table with Cervus results. Use
#' \code{\link{read.cervus}} to obtain it.
#' @param genotype \code{Genotype} or \code{matrix}. The object to take the
#' genotype calls from.
#' @param mm_cutoff \code{numeric}. The maximum number of mismatches between an
#' offspring and a parent, to count them as a true pair.
#' @param lod_cutoff \code{numeric}. The minumum lod score of an offspring-
#' parent pair, to count them as a true pair.
#'
#' @return \code{numeric}. A vector the same length as the number of markers,
#' with the relative number of true parent-offspring pair that match according
#' to each marker.
#'
#' @seealso \code{\link{write.cervus}}, \code{\link{write.candidates}},
#' \code{\link{read.cervus}}
#'
#' @examples
#' data("known_mother")
#' data("citrus_clean")
#' marker_match = Marker_Match_from_Cervus(known_mother, genotype=MxS, mm_cutoff=5,
#' lod_cutoff=10)

setGeneric("Marker_Match_from_Cervus",
           function(object, genotype, mm_cutoff=Inf, lod_cutoff=-Inf) {
             standardGeneric("Marker_Match_from_Cervus")
             })
setMethod("Marker_Match_from_Cervus", "data.frame",
          function(object, genotype, mm_cutoff=Inf, lod_cutoff=-Inf) {
  if (class(genotype) == "Genotype") genotype = getCalls(genotype)
  if (class(genotype) != "matrix") stop("genotype must be an object of class \'Genotype\' or a matrix of calls")
  if (class(object) != "data.frame") stop("object must be an output table of Cervus analysis with known mother.")
  if (sum(names(object) %in% c("Mother_ID", "Mother.ID")) == 0)
    stop("object must be an output table of Cervus analysis with known mother.\nIt must have a column called \'Mother_ID\' or \'Mother.ID\'")
  names(object) = gsub("\\.1", "*1", names(object))
  names(object) = gsub("\\.", "_", names(object))
  names(object) = gsub("\\*1", ".1", names(object))

  if (!all(object$Offspring_ID %in% colnames(genotype)))
    stop("Not all Offspring_IDs appear in the genotype sample names.")
  if (!all(object$Mother_ID %in% colnames(genotype)))
    stop("Not all Offspring_IDs appear in the genotype sample names.")

  mm_column = intersect(grep("mismatching", names(object)),
                        grep("Pair", names(object)))[1]
  object = object[object[[mm_column]]<=mm_cutoff, ]
  lod_column = intersect(grep("LOD", names(object)), grep("Pair", names(object)))[1]
  object = object[object[[lod_column]]>=lod_cutoff, ]
  object = object[!duplicated(object$Offspring_ID), , drop=F]

  pair_common_alleles<-function(v1, v2) {
    s = (v1=="XX" & v2 =="YY") | (v1=="YY" & v2 =="XX")
    return(!s)
  }

  match = sapply(1:nrow(object), function(i) {
    offspring = object$Offspring_ID[i]
    mother = object$Mother_ID[i]
    V1 = genotype[, offspring]
    V2 = genotype[, mother]
    pair_common_alleles(V1, V2)
  })
  return(rowSums(match, na.rm=T)/rowSums(!is.na(match)))

})

##-----------------------------------------------------------------------------
#' Adding sample information to the cervus result table
#'
#' @param cervus \code{data.frame}. The cervus results table, obtained with
#' \code{\link{read.cervus}}.
#' @param SampleInfo \code{SampleInfo}. Sample information table
#' @param add.columns \code{character}. Names of columns in \code{SampleInfo}
#' to add to the Cervus results table
#'
#' @return \code{data.frame}. The cervus results table including the new
#' columns.

setGeneric("addSampleInfo2cervus", function(cervus, SampleInfo, add.columns) {
  standardGeneric("addSampleInfo2cervus")
  })
setMethod("addSampleInfo2cervus", "data.frame",
          function(cervus, SampleInfo, add.columns) {
  names(cervus) = gsub("\\.", "_", names(cervus))
  names(cervus) = gsub(" +", "_", names(cervus))
  ID_columns = names(cervus)[grep("ID", names(cervus))]
  for (i in ID_columns) {
    IDind = which(names(cervus) ==i)
    ind = which((!is.na(cervus[[IDind]])) & nchar(cervus[[IDind]])>0)
    for (add in add.columns) {
      a = character(nrow(cervus))
      a[ind] = sapply(cervus[[i]][ind], function(x) SampleInfo[[add]][rownames(SampleInfo) ==x])
      indadd = which(add.columns==add)
      if ((IDind+indadd-1) < length(cervus))
        cervus = data.frame(cervus[1:(IDind+indadd-1)], a, cervus[(IDind+indadd):length(cervus)]) else
          cervus = data.frame(cervus, a)
      names(cervus)[IDind+indadd] = gsub("ID", add, names(cervus)[IDind])
    }
  }
  return(cervus)
})

#' Change genotype call allele format from call from X and Y to A, C, T and G
#'
#' In order to use some downstream tools, it is sometimes required to transorm
#' the coding of the genotype calls from "XX", "XY" and "YY" to nucleotide
#' based, like "CC", "GT" and so on.
#'
#' @param object \code{Genotype}. The object to change the genotype call
#' encoding in.
#' @param alleles \code{character}. A named character vector of length
#' \code{MakrerNo(object, Masked=F)} with the 2 alleles of each marker. For a
#' marker with the two alleles, "A" and "G" the value will be "AG". The names
#' of the values are the marker names, which should be the same as the
#' \code{Genotype} object mareker names.
#' @return \code{Genotype}. The same object but with genotype calls specified
#' in nucleotides.
#'
#' @examples
#' data("citrus_clean")
#' Alleles = getMarkerInfo(MxS, Masked=FALSE)$Alleles
#' names(Alleles) <- MarkerNames(MxS, Masked=F)
#' nucMxS=xy2atgc(MxS, alleles=Alleles)

setGeneric("xy2atgc", function(object,alleles=character()) {
  standardGeneric("xy2atgc")
  })
setMethod("xy2atgc","Genotype",function(object,alleles=character()){
  if (MarkerNo(object,Masked=FALSE)!=length(alleles))
    stop("alleles length should be as long as object rows")
  if (!all(names(alleles) == MarkerNames(object,Masked=FALSE))) stop("not all snps exists")
  nMxS=NULL
  for(SNP in MarkerNames(object,Masked=FALSE)) {
    nMxS = rbind(nMxS, chartr("XY",alleles[SNP],getCalls(object,MaskedMarker=FALSE, MaskedSample=F)[SNP,]))
  }
  object@Calls<-nMxS
  object
})

##-----------------------------------------------------------------------------

#' Write an input file for 'STRUCTURE' software
#'
#' STRUCTURE is a free software package for using multi-locus genotype data to
#' investigate population structure. It infers the presence of distinct
#' populations, assigns individuals to populations, identifies migrants and
#' admixed individuals, and estimates population allele frequencies. It is
#' therefore desireable for your project to write an input for it, and read its
#' output.
#' \code{write.structure} writes an input file for STRUCTURE.

#' @param object \code{Genotype}. The object to use as input for STRUCTURE.
#' Masked samples and markers will not be included in this input.
#' @param fn \code{character}. The file name for the STRUCTURE input file.
#' @param type \code{character}. The format (nuclotides or numbers) to use for
#' the STRUCTURE input file.
#' @param popdata \code{numeric}. An optional argument. A named numeric vector.
#' The population each sample belongs to. Use 0 for samples of unknown
#' population. The names are the names of samples in \code{object}.
#'
#' @return \code{NULL}.
#' @examples
#' data("citrus_clean")
#' Alleles = getMarkerInfo(MxS, Masked=FALSE)$Alleles
#' names(Alleles) <- MarkerNames(MxS, Masked=F)
#' nucMxS=xy2atgc(MxS, alleles=Alleles)
#' write.structure(nucMxS, fn="Structure_input.txt", type="numbers")

setGeneric("write.structure",
           function(object,fn,type=c("nuc","numbers"), popdata=NULL) {
             standardGeneric("write.structure")
             })
setMethod("write.structure","Genotype",
          function(object, fn, type=c("nuc","numbers"), popdata=NULL) {
  if (!is.character(fn)) stop("file name is required!\n")
  MxS = getCalls(object)
  missing = unique(as.vector(MxS)[!grepl("[ATCG][ATCG]",as.vector(MxS))])
  offset = ncol(MxS)
  SxM_struct = rbind(t(MxS),t(MxS))
  for(R in 1:offset){
    SxM_struct[R,][SxM_struct[R,]%in%missing] = "-9"
    SxM_struct[R+offset,][SxM_struct[R+offset,]%in%missing] = "-9"
    SxM_struct[R,] =
      stringr::str_replace(SxM_struct[R,],"^[ATGC]","")
    SxM_struct[R+offset,] =
      stringr::str_replace(SxM_struct[R+offset,],"[ATGC]$","")
  }
  o = lapply(SampleNames(object), function(x) which(rownames(SxM_struct) ==x))
  if (!all(sapply(o, length)==2)) stop("A probelm with sample names")
  if (!all(sapply(o, diff)==SampleNo(object))) stop("A probelm with sample names")
  o = unlist(o)
  # SxM_struct = SxM_struct[order(rownames(SxM_struct)),]
  SxM_struct = SxM_struct[o, ]
  M = apply(SxM_struct,2,function(x) chartr("ACGT","1234",x))
  if (length(popdata)) {
    popdata = popdata[rownames(M)]
    popflag = rep(0, length(popdata))
    popflag[popdata>0] <- 1
    M = cbind(popdata, popflag, M)
    colnames(M)[1:2] <- ""
  }
  write.table(switch(type,nuc=SxM_struct,numbers=M), file=fn, sep="\t", row.names=TRUE, col.names=NA, quote = F)
})

##-----------------------------------------------------------------------------

#' writing an executable file to run the 'structure' commands
#'
#' Writing an 'bash' file with the command to run a series of STRUCTURE
#' simulations, and files with parameters for STRUCTURE. This allows you to run
#' a number of simulations for each K (the number of ancestry populations),
#' which is recommended for the analysis of STRUCTURE results.
#'
#' @param object \code{Genotype}. The object containing the samples and markers
#' to run STRUCTURE on.
#' @param fn \code{character}. The file name for the file with the STRUCTURE commands.
#' @param k.range \code{numeric}. A vector with 2 values: the minimum and
#' maximum 'K's to run. The program will run all the values of K between and
#' including the minimum and maximum.
#' @param n.runs \code{numeric}. The number of runs to run each K. In any case,
#' there will not be more than one run with K=1, because for this K the results
#' are degenerated and deterministic.
#' @param path \code{character}. The directory to write the STRUCTURE files in.
#' @param parameter_list \code{list}. An optional list of additional parameters
#' for STRUCTURE.
#'
#' @details After running \code{write.structure.bin}, a new file named
#' \code{fn} is formed in your working directory, with the calls to STRUCTURE.
#' to run it, you first need to turn it to an executable file, for example by
#' typing \code{chmod a=rwx <fn>} in your Unix command line terminal.  You can
#' then run it by typing in your Unix command line terminal
#'     \code{./<fn>}.
#'
#' @examples
#' data("citrus_clean")
#' Alleles = getMarkerInfo(MxS, Masked=FALSE)$Alleles
#' names(Alleles) <- MarkerNames(MxS, Masked=F)
#' nucMxS=xy2atgc(MxS, alleles=Alleles)
#' write.structure.bin(MxS, fn="structure_sim", k.range = c(1, 10), n.runs=20)

setGeneric("write.structure.bin",
           function(object, fn="structure_sim", k.range = c(1, 10), n.runs=20,
                    path="Structure", parameter_list=list()) {
             standardGeneric("write.structure.bin")
             })

setMethod("write.structure.bin","Genotype",
          function(object, fn="structure_sim", k.range = c(1, 10), n.runs=20,
                   path="Structure", parameter_list=list()) {
  myfile = file(fn, open="w")
  writeLines(c("#!/bin/bash", ""), con=myfile)
  if (!is.null(path)) {
    if (class(path) == "character") {
      writeLines(c(paste("mkdir -p", path), ""), con=myfile)
      writeLines(c(paste("mkdir -p", paste(path, "Results", sep="/")), ""),
                 con=myfile)
    } else stop("\'path\' should be an object of class \'character\'")
  } else writeLines(c("mkdir -p Results", "", ""), con=myfile)
  mainparams = "mainparams"
  if (!is.null(path)) mainparams = paste(path, mainparams, sep="/")
  writeLines(paste("mainparams=\"", mainparams, "\"", sep=""), con=myfile)
  writeLines("cat <<EOF > $mainparams", con=myfile)
  if ("outfile" %in% names(parameter_list)) outfile = parameter_list$outfile else
    outfile = ""
  outfile = paste(getwd(), path, outfile, sep="/")
  writeLines(paste("#define OUTFILE", outfile), con=myfile)
  if ("infile" %in% names(parameter_list)) infile = parameter_list$infile else {
    infile = list.files(pattern="Structure_input.txt", ignore.case=T)
    if (length(infile) == 0) stop("Found no genotype call file with suitable name.")
    if (length(infile)>1) stop("Found more than one genotype call file with suitable name.")
  }
  infile = paste(getwd(), infile, sep="/")
  writeLines(paste("#define INFILE", infile), con=myfile)
  writeLines(c(paste("#define NUMINDS", SampleNo(object)),
               paste("#define NUMLOCI", MarkerNo(object))), con=myfile)
  parameter = c(label=1, popdata=0, popflag=0, locdata=0, phenotype=0,
              markernames=1, mapdistances=0, onerowperind=0, phaseinfo=0,
              phased=0, recessivealleles=0, extracols=0, missing=-9, ploidy=2,
              maxpops=1, burnin=10000, numreps=100000, noadmix=0, linkage=0,
              usepopinfo=0, locprior=0, inferalpha=1, alpha=1.0, popalphas=0,
              unifprioralpha=1, alphamax=10.0, alphapropsd=0.025, freqscorr=1,
              onefst=0, fpriormean=0.01, fpriorsd=0.05, inferlambda=0,
              lambda=1.0, computeprob=1, pfrompopflagonly=0, ancestdist=0,
              startatpopinfo=0, metrofreq=10, updatefreq = 1)
  if (length(parameter_list)) {
    for (a in names(parameter_list))
      parameter[tolower(a)] <- parameter_list[[a]]
  }
  for (a in names(parameter))
    writeLines(paste("#define ", toupper(a), parameter[a]), con=myfile)
  writeLines("EOF", con=myfile)
  extraparams = "extraparams"
  if (!is.null(path)) extraparams = paste(path, extraparams, sep="/")
  writeLines(paste("extraparams=\"", extraparams, "\"", sep=""), con=myfile)
  writeLines(c("touch $extraparams", ""), con=myfile)

  if (k.range[1] ==1) {
    k.range[1] = 2
    writeLines("  echo \"now running K: 1 and iteration: 1\"", con=myfile)
    structure_command = "	structure -K 1 -m $mainparams -e $extraparams -o 10000_100000/Results/K_1_run_1"
    if (!is.null(path)) structure_command = paste("	structure -K 1 -m $mainparams -e $extraparams -o ", path, "/Results/K_1_run_1", sep="")
    writeLines(structure_command, con=myfile)
    writeLines("", con=myfile)
  }
  if (length(k.range) == 1)
    writeLines(paste("k=", k.range, sep=""), con=myfile) else {
      writeLines(paste("for ((k=",  k.range[1], ";k<=", k.range[2], ";k++))", sep=""), con=myfile)
      writeLines("do", con=myfile)
    }
  writeLines(paste("  for ((i=1;i<=", n.runs, ";i++))", sep=""), con=myfile)
  writeLines("  do", con=myfile)
  writeLines("  echo \"now running K:  $k and iteration: $i\"", con=myfile)
  structure_command = "	structure -K $k -m $mainparams -e $extraparams -o 10000_100000/Results/K_${k}_run_${i}"
  if (!is.null(path)) structure_command = paste("	structure -K $k -m $mainparams -e $extraparams -o ", path, "/Results/K_${k}_run_${i}", sep="")
  writeLines(c(structure_command, "	done"), con=myfile)
  if (length(k.range)>1) writeLines("done", con=myfile)

  close(con=myfile)
})

##-----------------------------------------------------------------------------

#' Getting the matrix of sub-population ancestry from STRUCTURE output
#'
#' STRUCTURE calculates for each individual sample, how much of its genome
#' belongs to each one of K sub-populations, or so-called clusters. This
#' function extracts this information from the STRUCTURE output file, and
#' returns a matrix with K columns and n rows, n being the number of samples.
#'    This function is used by \code{\link{Sum_Structure_Results}}, which sums
#' over multiple simulations with the same data and K.
#'
#' @param f \code{character}. The file name
#' @param path \code{character}. The path to the directory where the file is.
#' @param k \code{numeric}. The number of sub-population or 'clusters'. This
#' number is a pre-defined parameter required by STRUCTURE.
#' @param print_file \code{logical}. Should the name of the file be printed?
#'
#' @return the matrix of relative ancestry.

get_matrix = function(f, path, k, print_file) {
  s = scan(paste(path, f, sep="/"), what="", sep="\n", quiet=T)
  nind = grep("individuals", s)[1]
  nind = gsub("[A-Za-z ]+", "", s[nind])
  if (is.na(nind) || length(nind)==0)
    stop(paste("The line with the number of individuals is missing in file", f))
  nind = as.integer(nind)
  g = grep("Label", s)
  g = g[grep("(%Miss)", s[g])]
  g = g[grep("Pop", s[g])]
  if (is.na(g) || length(g)==0)
    stop(paste("The headline 'Label (%Miss) Pop' is missing in file", f))
  s = s[g + 1:nind]
  s = gsub("^ +", "", s)
  s = gsub("^[0-9]+ +", "", s)
  s = strsplit(s, split=" +")
  A = data.frame(Sample_name=sapply(s, "[[", 1))
  A$pre_defined = sapply(s, "[[", 3)
  M = matrix(0, nrow=nrow(A), ncol=k)
  x = s[!(A$pre_defined %in% 1:k)]
  if (length(x)) {
    x = lapply(x, function(v) v[4+1:k])
    x = sapply(x, as.numeric)
    x = t(x)
  }
  M[!(A$pre_defined %in% 1:k), ] <- x
  for (i in 1:k)
    M[A$pre_defined==i, i] <- 1
  if (sum(is.na(M)))
    stop(paste("There is a problem with the ancestry matrix in file", f))
  if (print_file) cat(f, "\n")
  return(data.frame(A, M))
}

##-----------------------------------------------------------------------------

#' Summing up multiple STRUCTUE simulations with predifined populations
#'
#' STRUCTURE calculates for each individual sample, how much of its genome
#' belongs to each one of K sub-populations, or so-called clusters.
#' \code{Sum_Structure_Results} sums up the results of a number of STRUCTURE
#' simulations, all ran with on the same data and with the same K. It can only
#' be used if STRCTURE used some samples with pre-defined populations.
#' Otherwise the 'phase' of the STRUCTURE output will be different in each run,
#' that is, population number 1 from one simulation, may have a different
#' number in another population. In this case, use tools like CLUMPP or CLUMPAK
#' which solve the phase problem first.
#'
#' @param path \code{character}. The path to the directory where the STRUCTURE
#' output files are.
#' @param postfix \code{character}. The postfix of the file names. If STRUCTURE
#' run used \code{"write.structure.bin"}, then the file name is
#' \code{"K_<K>_run_<itteration number><postfix>"}.
#' @param k \code{numeric}. The number of sub-population or 'clusters'. This
#' number is a pre-defined parameter required by STRUCTURE.
#' @param print_file \code{logical}. Should the name of the each file read be
#' printed?
#'
#' @return the matrix of mean relative ancestry, summing over all simulations.
#'
#' @seealso \code{\link{get_matrix}}
#'
Sum_Structure_Results = function(path, postfix="_f", k, print_file=T) {
  if (missing(k)) stop("Please provide k, the number of ancestral populations")
  cat("This function is for summing up several runs of STRUCTURE, with some samples having PRE-DEFINED POPULATIONS")
  files = list.files(path, pattern=paste(postfix, "$", sep=""))
  if (length(files) ==0) stop("No appropriate files found")
  results = lapply(files, get_matrix, path=path, k=k, print_file=print_file)
  meanM = results[[1]][, 2+1:k]
  for (i in 2:length(results))
    meanM =meanM + results[[i]][, 2+1:k]
  meanM = meanM/length(results)
  x = results[[1]]
  x[, 2+1:k] <- meanM
  names(x)[2+1:k] <- paste("cluster", 1:k, sep="_")
  names(x)[names(x) == "pre_defined"] <- "user-assigned_population"
  return(x)
}

##-----------------------------------------------------------------------------

setClass(Class="STRUCTURE",
         representation(name = "character", path = "character", param= "list",
                        data= "matrix", loci = "character",
                        individuals = "character")
)

setMethod("initialize","STRUCTURE",function(.Object,name, path) {
  .Object@name =  name
  .Object@path = path
  .Object = read.project(.Object,name,path)
  .Object
})
setGeneric("lociNo",function (object){standardGeneric("lociNo")})
setMethod("lociNo","STRUCTURE",function (object) {
  length(object@loci)
})
setGeneric("lociNames",function (object){standardGeneric("lociNames")})
setMethod("lociNames","STRUCTURE",function (object){
  object@loci
})
setGeneric("IndividualsNo",function (object){standardGeneric("IndividualsNo")})
setMethod("IndividualsNo","STRUCTURE",function (object){
  length(object@individuals)/as.numeric(object@param[["PLOIDY"]])
})
setGeneric("IndividualsNames",function (object){
  standardGeneric("IndividualsNames")
  })
setMethod("IndividualsNames","STRUCTURE",function (object) {
  unique(object@individuals)
})
setGeneric("read.structure.results",
           function(object, path=NULL, name = NULL, itr = 20,
                    k.range = c(1,10), show.itr=FALSE) {
             standardGeneric("read.structure.results")
           })
setMethod("read.structure.results","STRUCTURE",function(object, path=NULL, name=NULL, itr=20, k.range=c(1,10), show.itr=FALSE){
  fp = paste(object@path, path, paste(object@name,name, sep=""),sep="/")
  if (!file.exists(fp)) stop("can not find directory\n",fp)

  ####  Q-values #####
  Q.values = list()
  for(k in k.range[1]:k.range[2]){
    Q.values[[k]]= array(NA,dim=c(IndividualsNo(object),k,itr))
    if (show.itr) cat("K:",k,"\n")
    for(j in 1:itr){
      if (show.itr) cat("Itr:",j,"\n")
      fn = paste("K", k, "run", j, "f",sep="_")
      topFile = readLines(con=paste(fp, fn, sep=""), n=400)
      s.table = grep("Inferred ancestry of individuals",topFile) + 1
      e.table = grep("Estimated Allele Frequencies in each cluster",topFile)-3
      # widths=c(3,11,-13,rep(6,k))
      Q.table =read.table(file=paste(fp, fn, sep=""),header=FALSE,skip=s.table, nrows=(e.table-s.table),row.names=1)
      ID= Q.table[,1]
      pop = paste("P",1:k,sep="")
      iterations = paste("Itr",1:itr,sep="")
      dimnames(Q.values[[k]]) = list(ID,pop,iterations)
      Q.values[[k]][,,j]= as.matrix(Q.table[,-(1:3),drop=FALSE])
    }
  }
  Q.mean.values = list()

  for(k in k.range[1]:k.range[2]){
    if (k==1){
      Q.mean.values[[k]] <- apply(Q.values[[k]],c(1,2),mean)
    } else {
      Q.mean.values[[k]] <- t(apply(apply(Q.values[[k]],c(1,3),sort),c(1,2),mean))
      colnames(Q.mean.values[[k]]) = paste("P",1:k,sep="")
      #reorder the Q means according to subpopulations
      for (r in 1:nrow(Q.mean.values[[k]])){
        Q.mean.values[[k]][r,] =
          Q.mean.values[[k]][r,][rank(Q.values[[k]][r,,1])]
      }
    }
  }
  Q.mean.values
})
setGeneric("read.project",function(object, name = NULL, path = NULL) {
  standardGeneric("read.project")
  })
setMethod("read.project","STRUCTURE",
          function(object, name = NULL, path = NULL) {
  fn = paste(path,name,sep="/")
  if (!file.exists(fn)) stop("can not find directory\n",fn)
  tab=read.table(paste(fn,paste(name,"spj",sep="."),sep="/"),header=F,stringsAsFactor=F)
  for (r in 1:nrow(tab)){
    object@param[[tab[r,1]]]=tab[r,2]
  }
  tab=read.delim(paste(fn,"project_data",sep="/"),stringsAsFactor=F)
  object@data = as.matrix(tab[,-1])
  object@individuals = as.character(tab[,1])
  object@loci = colnames(tab)
  object
})

setGeneric("read.structure.input", function(object, name = NULL) {
  standardGeneric("read.structure.input")
  })
setMethod("read.structure.input","STRUCTURE",
          function(object, name="mainparams") {
  path = object@path
  fn = paste(path,name,sep="/")
  if (!file.exists(fn)) stop("can not find directory\n",fn)
  params = scan(fn, what="", sep="\n", quiet=T)
  params = params[grep("^#define", params)]
  params = gsub("^#define +", "", params)
  params = strsplit(params, split=" +")
  Params = lapply(params, function(x) x[-1])
  names(Params) = sapply(params, "[[", 1)
  object@param = Params
  tab = read.delim(file=Params$INFILE, stringsAsFactor=F)
  object@individuals = as.character(tab[,1])
  object@loci = colnames(tab)[-1]
  Data =  as.matrix(tab[,-1])
  Data[Data==1] = "A"
  Data[Data==2] = "C"
  Data[Data==3] = "G"
  Data[Data==4] = "T"
  Data[Data == Params$MISSING] = NA
  object@data = Data
  object
})
setGeneric("read.simulation",
           function(object, name = NULL, itr = 20, sim.range = c(1,10),
                    show.itr=FALSE) {
  standardGeneric("read.simulation")
  })
setMethod("read.simulation","STRUCTURE",
          function(object, name = NULL, itr = 20, sim.range = c(1,10),
                   show.itr=FALSE){

  fp = paste(object@path,object@name,name,sep="/")
  if (!file.exists(fp)) stop("can not find directory\n",fp)

  fn.prefix = paste(fp,"Results","sim_5000_50000_run",sep="/")

  ####  Q-values #####
  Q.values = list()
  for(k in sim.range[1]:sim.range[2]){
    Q.values[[k]]= array(NA,dim=c(IndividualsNo(object),k,itr))
    if (show.itr) cat("K:",k,"\n")
    for(j in 1:itr){
      if (show.itr) cat("Itr:",j,"\n")
      fn = paste(fn.prefix,(k-1)*itr+j,"f",sep="_")
      topFile = readLines(con=fn,n=400)
      s.table = grep("Inferred ancestry of individuals",topFile)+1
      e.table = grep("Estimated Allele Frequencies in each cluster",
                     topFile)-3
      # widths=c(3,11,-13,rep(6,k))
      Q.table =read.table(file=fn,header=FALSE,skip=s.table, nrows=(e.table-s.table),row.names=1)
      if (show.itr) cat("run:",(k-1)*itr+j,"\n")
      ID= Q.table[,1]
      pop = paste("P",1:k,sep="")
      iterations = paste("Itr",1:itr,sep="")
      dimnames(Q.values[[k]]) = list(ID,pop,iterations)
      Q.values[[k]][,,j]= as.matrix(Q.table[,-(1:3),drop=FALSE])
    }
  }

  Q.mean.values = list()

  for(k in sim.range){
    if (k==1){
      Q.mean.values[[k]] <- apply(Q.values[[k]],c(1,2),mean)
    } else {
      Q.mean.values[[k]] <- t(apply(apply(Q.values[[k]],c(1,3),sort),c(1,2),mean))
      colnames(Q.mean.values[[k]]) = paste("P",1:k,sep="")
      #reorder the Q means according to subpopulations
      for (r in 1:nrow(Q.mean.values[[k]])){
        Q.mean.values[[k]][r,] =
          Q.mean.values[[k]][r,][rank(Q.values[[k]][r,,1])]
      }
    }
  }
  Q.mean.values
})
setGeneric("read.simulation.summary",
           function(object, name = NULL) {
             standardGeneric("read.simulation.summary")
             })
setMethod("read.simulation.summary","STRUCTURE",function(object, name = NULL){
  tmp=readLines(paste(object@path,object@name,"simulation_summary.txt",sep="/"))
  headers = unlist(strsplit(tmp[1],split="  +"))
  records = strsplit(tmp[2:length(tmp)],split=" +")
  summary.table=as.data.frame(
    do.call(rbind,records),stringsAsFactors=FALSE)

  headers =  gsub("\\?","alpha",headers)
  headers =  gsub(" +","",headers)
  headers =  gsub("\\[|\\(","_",headers)
  headers =  gsub("\\]|\\)","",headers)

  colnames(summary.table) = headers

  for (cl in 3:ncol(summary.table)){
    summary.table[[cl]] = as.numeric(summary.table[[cl]])
  }
  summary.table
})

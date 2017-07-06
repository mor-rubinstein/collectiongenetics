#' A S4 class for a matrix genotype frequencies.
#'
#' @slot genotype_freq \code{matrix}: A matrix of genotype frequencies. Each
#' row corresponds to a marker, and the columns list the number of individual
#' samples with "XX", "XY" and "YY" genotypes.
#'
setClass(Class="AlleleStat",
         representation(genotype_freq= "matrix")
)

##-----------------------------------------------------------------------------
setMethod("show","AlleleStat",function(object){
  cat("Genotype frequency: \n")
  print(head(object@genotype_freq))
  cat("\n")
  cat("Observed genotype counts summary: \n")
  avg=apply(object@genotype_freq,2,mean)
  print(rbind(mean_counts=avg,propotion_of_mean_counts=avg/sum(avg)))
  cat("\n")
  cat("average of allele statistics: \n")
  cat("pX: ",mean(pX(object)),"\n")
  cat("pY: ",mean(pY(object)),"\n")
  cat("Hexp: ",mean(Hexp(object)),"\n")
})
##-----------------------------------------------------------------------------

#' Genotype frequencies
#'
#' @param object \code{AlleleStat}. The object to get the genotype frequencies
#' from.
#' @return \code{matrix}. The matrix of genotype frequencies. Each row
#' corresponds to a marker, and the columns list the number of individual
#' samples with "XX", "XY" and "YY" genotypes.
#' @examples
#' data("Mango_raw")
#' MarkersStat = AlleleCounting(MxS)
#' GF = genotype_freq(MarkersStat)
#' head(GF)

setGeneric("genotype_freq", function(object) {
  standardGeneric("genotype_freq")
  })
setMethod("genotype_freq", "AlleleStat", function(object) {
  object@genotype_freq
})

##-----------------------------------------------------------------------------

#' Counting genotype frequencies
#'
#' Counting the frequency of each genotype, per each marker
#' @param object. \code{Genotype}. The object to count the genotypes in.
#' @return \code{AlleleStat}. An object contianing the matrix of genotype
#'  frequencies. Each row corresponds to a marker, and the columns list the
#'  number of individual samples with "XX", "XY" and "YY" genotypes.
#' @examples
#' data("Mango_raw")
#' MarkersStat = AlleleCounting(MxS)
#' GF = genotype_freq(MarkersStat)
#' head(GF)

setGeneric("AlleleCounting",function(object) {
  standardGeneric("AlleleCounting")
  })
setMethod("AlleleCounting","Genotype",function(object){
  alleles = c("XX","XY","YY")
  y=t(apply(getCalls(object),1,function(x) table(x)[alleles]))
  colnames(y)= alleles
  y[is.na(y)]=0
  return(new("AlleleStat",genotype_freq= y))
})

##-----------------------------------------------------------------------------
#' Number of individuals with a genotype call
#'
#' Counting the number of samples with valid genotype calls, per each marker.
#'
#' @param object. \code{AlleleStat}. The object to count the samples in.
#' @return \code{numeric}. A numeric vector, with the number of genotype calls
#' available for each marker.
#' @examples
#' data("citrus_clean")
#' MarkersStat = AlleleCounting(MxS)
#' pz = PopSize(MarkersStat)
#' SampleNo(MxS)
#' summary(pz)

setGeneric("PopSize",function(object) {
  standardGeneric("PopSize")
  })
setMethod("PopSize","AlleleStat",function(object){
  apply(object@genotype_freq,1, sum)
})

##-----------------------------------------------------------------------------
#' Allele 1 frequency
#'
#' Calculates the frequency of the 1st, per each marker.
#'
#' @param object. \code{AlleleStat}. The object to calculate allele frequency
#' from.
#' @return \code{numeric}. A numeric vector, with the frequency of the 1st
#' allele for each marker.
#' @seealso \code{\link{pY}}
#' @examples
#' data("citrus_clean")
#' MarkersStat = AlleleCounting(MxS)
#' px = pX(MarkersStat)

setGeneric("pX",function(object) {standardGeneric("pX")})
setMethod("pX","AlleleStat",function(object){
  apply(object@genotype_freq,1,
        function(x) (x["XX"]*2+x["XY"])/(sum(x)*2))
})

##-----------------------------------------------------------------------------
#' Allele 2 frequency
#'
#' Calculates the frequency of the 2nd, per each marker.
#'
#' @param object. \code{AlleleStat}. The object to calculate allele frequency
#' from.
#' @return \code{numeric}. A numeric vector, with the frequency of the 2nd
#' allele for each marker.
#' @seealso \code{\link{pX}}
#' @examples
#' data("citrus_clean")
#' MarkersStat = AlleleCounting(MxS)
#' py = pY(MarkersStat)

setGeneric("pY",function(object) {standardGeneric("pY")})
setMethod("pY","AlleleStat",function(object){
  apply(object@genotype_freq,1,
        function(x) (x["YY"]*2+x["XY"])/(sum(x)*2))
})

##-----------------------------------------------------------------------------
#' The observed frequency of heterozygous
#'
#' Calculates the relative frequency of heterozygous, per each marker.
#'
#' @param object. \code{AlleleStat}. The object to calculate heterozygousity
#' frequency from.
#' @return \code{numeric}. A numeric vector, with the frequency of the
#' heterozygous for each marker.
#' @seealso \code{\link{pX}} and \code{\link{pY}} for allele frequencies,
#' \code{\link{Hexp}} for expected heterozygousity.
#' @examples
#' data("citrus_clean")
#' MarkersStat = AlleleCounting(MxS)
#' h = Hobs(MarkersStat)

setGeneric("Hobs",function(object) {standardGeneric("Hobs")})
setMethod("Hobs","AlleleStat",function(object){
  apply(object@genotype_freq,1,function(a) a["XY"]/sum(a))
})

##-----------------------------------------------------------------------------
#' The expected frequency of heterozygous
#'
#' Calculates the expected relative frequency of heterozygous, per each marker.
#'
#' @param object. \code{AlleleStat}. The object to calculate heterozygousity
#' frequency from.
#' @return \code{numeric}. A numeric vector, with the expected frequency of
#' heterozygous for each marker.
#' @seealso \code{\link{pX}} and \code{\link{pY}} for allele frequencies,
#' \code{\link{Hobs}} for observed heterozygousity.
#' @examples
#' data("citrus_clean")
#' MarkersStat = AlleleCounting(MxS)
#' h = Hexp(MarkersStat)

setGeneric("Hexp",function(object) {standardGeneric("Hexp")})
setMethod("Hexp","AlleleStat",function(object){
  1-(pX(object)^2+pY(object)^2)
})

##-----------------------------------------------------------------------------
#' Fixation index
#'
#' Fixation index, as the reduction in heterozygousity, compared to expectd
#' heterozygousity in a random mating population
#'
#' @param object. \code{AlleleStat}. The object to calculate the fixation index
#' from.
#' @return \code{numeric}. A numeric vector, with the fixation index for each
#' marker.
#' @seealso \code{\link{pX}} and \code{\link{pY}} for allele frequencies,
#' \code{\link{Hobs}} for observed heterozygousity and  \code{\link{Hexp}} for
#' observed heterozygousity.
#' @examples
#' data("citrus_clean")
#' MarkersStat = AlleleCounting(MxS)
#' fs = Fs(MarkersStat)

setGeneric("Fs",function(object) {standardGeneric("Fs")})
setMethod("Fs","AlleleStat",function(object){
  (Hexp(object)-Hobs(object))/Hexp(object)
})

##-----------------------------------------------------------------------------
#' Polymorphic information content
#'
#' @param object. \code{AlleleStat}. The object to calculate the polymorphic
#' information content from.
#' @return \code{numeric}. A numeric vector, with the polymorphic information
#' content for each marker.
#' @examples
#' data("citrus_clean")
#' MarkersStat = AlleleCounting(MxS)
#' pic = pic(MarkersStat)

setGeneric("PIC",function(object,type=1) {standardGeneric("PIC")})
setMethod("PIC","AlleleStat",function(object,type=1){
  if(!type[1]%in%c(1,2)) stop("type should be either 1 or 2")
  if(type==1)
    pic = 1-(pX(object)^2+pY(object)^2)
  if(type==2)
    pic =1-(pX(object)^2+pY(object)^2)-2*pX(object)^2*pY(object)^2
  return(pic)
})

##-----------------------------------------------------------------------------
#' Intercross pollination
#'
#' @param object \code{AlleleStat}. The object to calculate the polymorphic
#' information content from.
#' @param test \code{character}. a vector of genotype calls for the same
#' markers that exist in \code{object}.
#' @param N \code{numeric}. A single numeric value. The maximum number of
#' markers to test
#' @return \code{data.frame}. A table with maximum N rows corresponding to
#' markers, in which test is homozygous for the less frequent allele in the
#' object.
#' @examples
#' data("citrus_clean")
#' MarkersStat = AlleleCounting(MxS)
#' test = getCalls(MxS)[, 1]
#' ocpol = ocPol(MarkersStat, test=test, N=10)

setGeneric("ocPol",function(object,test=character(),N=numeric()) {
  standardGeneric("ocPol")
  })
setMethod("ocPol","AlleleStat",function(object,test=character(),N=numeric()){
  if (!is.character(test)|!any(names(table(test))%in%c("XX","YY","XY")))
    stop("test must be a character vector of XX, YY and XY")
  prop = cbind(pX(object),pY(object))
  colnames(prop) = c("pX","pY")
  if(nrow(prop)!=length(test)) stop("prop and test must be of the same length")
  MAf=data.frame(Prop=apply(prop[,c("pX","pY")],1,min),Allele=apply(prop[,c("pX","pY")],1,function(x) ifelse(x[1]<x[2],"X","Y")),stringsAsFactors=F)
  homo_loci = sapply(test,function(x) switch(x,XX="X",YY="Y",XY=NA,NA))
  loci_of_interest = apply(cbind(MAf$Allele,homo_loci),1,function(x) x[1]==x[2])
  loci_of_interest[is.na(loci_of_interest)] = FALSE
  Markers = MAf[loci_of_interest,]
  Markers = Markers[order(Markers$Prop),]
  if (N>nrow(Markers)) N=nrows(Markers)
  cat("The probability for intercross pollination is",prod(Markers[1:N,"Prop"]),"\n")
  return(Markers[1:N,])
})


#' A S4 class for a sub-populations.
#'
#' @slot Pop \code{list}: A named list of \code{AlleleStat} objects, each
#' specifying allele genotype frequencies of the same markers in a different
#' population.
#'
setClass(Class="Populations",
         representation(Pop= "list")
)

setMethod("initialize","Populations",function(.Object,Pop){
  .Object@Pop<-Pop
  if (length(.Object@Pop)==0) stop("A list of AlleleStat objects is mandatory")
  if(!all(sapply(.Object@Pop,function(x) is(x)=="AlleleStat")))  stop("List of 'AlleleStat' object is required ")
  if(is.null(names(.Object@Pop)))
    names(.Object@Pop)<- paste("P",1:length(.Object@Pop),sep="")
  .Object
})

##-----------------------------------------------------------------------------
#' Number of populatoins
#'
#' The number or names of populations in a \code{Populations} object
#'
#' @param object. \code{Populations}. The object with the populations.
setGeneric("PopNo",function(object) {standardGeneric("PopNo")})
setMethod("PopNo","Populations",function(object){
  length(object@Pop)
})
#' @rdname PopNo
setGeneric("PopNames",function(object) {standardGeneric("PopNames")})
setMethod("PopNames","Populations",function(object){
  object@PopNames
})

##-----------------------------------------------------------------------------
#' Heterozygosity in a global population containing several sub-populations
#'
#' Global observed (\code{Hi}) heterozygousity, global expected (\code{Ht})
#' heterozygousity and weighted mean expected heterozygousity among
#' sub-populations (\code{Hs})
#'
#' @param object. \code{Populations}. The object to calculate the
#' heterozygosity in.
#' @return \code{numeric}. A numeric vector, with the heterozygosity in each
#' marker.
#' @examples
#' data("citrus_clean")
#' group.size= round(SampleNo(MxS)/3)
#' Total = 1:SampleNo(MxS)
#' group.samples=list()
#' for (i in 1:3) {
#'   if (i<3) {
#'     group.samples[[i]]=sample(Total,group.size)
#'     Total = Total[!Total%in%group.samples[[i]]]
#'   } else group.samples[[i]]= Total
#' }
#' group.samples = lapply(group.samples, function(x) SampleNames(MxS)[x])
#' A1 = AlleleCounting(SetSampleFilter(MxS, unlist(group.samples[-1])))
#' A2 = AlleleCounting(SetSampleFilter(MxS, unlist(group.samples[-2])))
#' A3 = AlleleCounting(SetSampleFilter(MxS, unlist(group.samples[-3])))
#' RandomPop = new("Populations",Pop=list(A1, A2, A3))
#' hi = Hi(RandomPop)
#' ht = Ht(RandomPop)
#' hs = Hs(RandomPop)
#' summary(abs(hi - Hobs(AlleleCounting(MxS))))
#' summary(abs(ht - Hexp(AlleleCounting(MxS))))
#' summary(abs(hs - Hexp(AlleleCounting(MxS))))

setGeneric("Hi",function(object) {standardGeneric("Hi")})
setMethod("Hi","Populations",function(object){
  N = sapply(object@Pop,PopSize)
  H.obs = sapply(object@Pop,Hobs)
  H  = rowSums(H.obs*N)/rowSums(N)
  names(H)=rownames(H.obs)
  return(H)
})

#' @rdname Hi
setGeneric("Hs",function(object) {standardGeneric("Hs")})
setMethod("Hs","Populations",function(object){
  N = sapply(object@Pop,PopSize)
  H.exp = sapply(object@Pop,Hexp)
  H = rowSums(H.exp*N)/rowSums(N)
  names(H)=rownames(H.exp)
  return(H)
})
#' @rdname Hi
setGeneric("Ht",function(object) {standardGeneric("Ht")})
setMethod("Ht","Populations",function(object){
  N = sapply(object@Pop,PopSize)
  pX = sapply(object@Pop,pX)
  pX.bar = rowSums(pX*N)/rowSums(N)
  pY = sapply(object@Pop,pY)
  pY.bar = rowSums(pY*N)/rowSums(N)
  H = 1-(pX.bar^2+pY.bar^2)
  names(H)=rownames(pX)
  return(H)
})
##-----------------------------------------------------------------------------

#' Fixation index: F-statistics
#'
#' The inbreeding coefficient, as the reduction in heterozgous frequency
#' compared to the expected frequency, due to inbreeding within each sub-
#' population (Fis), due to the effect of subpopulations (Fst), and globaly
#' among the whole population (Fit)
#'
#' @param object. \code{Populations}. The object to calculate the
#' inbreeding coefficient in.
#' @return \code{numeric}. A numeric vector, with the inbreeding coefficient
#' in each marker.
#' @examples
#' in random sub-populations:
#' data("citrus_clean")
#' group.size= round(SampleNo(MxS)/3)
#' Total = 1:SampleNo(MxS)
#' group.samples=list()
#' for (i in 1:3) {
#'   if (i<3) {
#'     group.samples[[i]]=sample(Total,group.size)
#'     Total = Total[!Total%in%group.samples[[i]]]
#'   } else group.samples[[i]]= Total
#' }
#' group.samples = lapply(group.samples, function(x) SampleNames(MxS)[x])
#' A1 = AlleleCounting(SetSampleFilter(MxS, unlist(group.samples[-1])))
#' A2 = AlleleCounting(SetSampleFilter(MxS, unlist(group.samples[-2])))
#' A3 = AlleleCounting(SetSampleFilter(MxS, unlist(group.samples[-3])))
#' RandomPop = new("Populations",Pop=list(A1, A2, A3))
#' fis = Fis(RandomPop)
#' fst = Fst(RandomPop)
#' fit = Fit(RandomPop)
#' summary(fis)
#' summary(fst)
#' summary(fit)
#'
#' in real sub-populations:
#' data("Avocado_populations")
#' fis = Fis( )
#' fst = Fst(AvPop)
#' fit = Fit(AvPop)
#' summary(fis)
#' summary(fst)
#' summary(fit)

setGeneric("Fis",function(object) {standardGeneric("Fis")})
setMethod("Fis","Populations",function(object){
  (Hs(object)-Hi(object))/Hs(object)
})

#' @rdname Fis
setGeneric("Fst",function(object) {standardGeneric("Fst")})
setMethod("Fst","Populations",function(object){
  (Ht(object)-Hs(object))/Ht(object)
})

#' @rdname Fis
setGeneric("Fit",function(object) {standardGeneric("Fit")})
setMethod("Fit","Populations",function(object){
  (Ht(object)-Hi(object))/Ht(object)
})

##-----------------------------------------------------------------------------

#' Write a \code{gdsfmt} package gds file from \code{Genotype}
#'
#' Changing from a \code{Genotype} object to a gds file allows us to use
#' \code{gdsfmt} package functions
#'
#' @param object \code{Genotype}. The object to write to a gds file.
#' @param gds.fn \code{character}. The name of the gds file.
#' @param Chrom \code{character}. Chromosome names markers.
#' @param Pos \code{numeric}. Chromosomal positions markers.
#' @param Allele \code{character}. Marker alleles.
#'
#' @return \code{logical}
#'
setGeneric("xy2gds",
           function(object, gds.fn="SNP.gds", Chrom=NULL, Pos=NULL, Allele) {
             standardGeneric("xy2gds")
             })
setMethod("xy2gds", "Genotype",
          function(object, gds.fn="SNP.gds", Chrom=NULL, Pos=NULL, Allele=NULL) {
  gdsGenotype = getCalls(object)
  gdsGenotype[gdsGenotype=="XX"] = 0
  gdsGenotype[gdsGenotype=="XY"] = 1
  gdsGenotype[gdsGenotype=="YX"] = 1
  gdsGenotype[gdsGenotype=="YY"] = 2
  gdsGenotype[gdsGenotype=="No Call"] = 3
  gdsGenotype = apply(gdsGenotype, 2, as.numeric)
  snp.pos <- if(is.null(Pos)) 1:MarkerNo(object) else getMarkerInfo(object)[,Pos]
  chrom <- if(is.null(Chrom))
    c(rep(1:10,each=MarkerNo(object)%/%10),rep(11,MarkerNo(object)%%10)) else
      getMarkerInfo(object)[,Chrom]
  gfile <- gdsfmt::createfn.gds(gds.fn)
  gdsfmt::add.gdsn(gfile, "sample.id", SampleNames(object),
                   compress = "ZIP.max", closezip = TRUE)
  gdsfmt::add.gdsn(gfile, "snp.id",MarkerNames(object), compress = "ZIP.max",
                   closezip = TRUE)
  gdsfmt::add.gdsn(gfile, "snp.position", snp.pos, compress = "ZIP.max",
                   closezip = TRUE)
  gdsfmt::add.gdsn(gfile, "snp.chromosome", chrom, compress = "ZIP.max",
                   closezip = TRUE)
  if(!is.null(Allele))
    gdsfmt::add.gdsn(gfile, "snp.allele", Allele, compress = "ZIP.max",
                     closezip = TRUE)
  node.geno <- gdsfmt::add.gdsn(gfile, "genotype", gdsGenotype,
                                storage = "bit2", compress = "ZIP.max")
  gdsfmt::closefn.gds(gfile)
  return(invisible(TRUE))
})

##-----------------------------------------------------------------------------

#' Write a \code{SNPRelate} package \code{snpgds} file from \code{Genotype}
#'
#' Changing from a \code{Genotype} object to a gds file allows us to use
#' \code{gdsfmt} package functions
#'
#' @param object \code{Genotype}. The object to write to a gds file.
#' @param filename \code{character}. The name of the snpgds file.
#' @return \code{character}. The name of the file created.
#'
#' @examples
#' data("Genotype_43X80")
#' write_snpgds_file(Genotype, filename="Genotype.gds")

write_snpgds_file = function(Genotype, filename="Genotype.gds") {
  g = getCalls(Genotype)
  A = matrix(data=rep(3, length(g)), nrow=nrow(g), ncol=ncol(g), dimnames=dimnames(g))
  A[g=="XX"] <- 0
  A[g=="XY"] <- 1
  A[g=="YY"] <- 2
  SNPInfo = getMarkerInfo(Genotype)
  chr = grep("chr", names(SNPInfo), ignore.case=T)
  if (length(chr)) {
    if (length(chr) > 1)
      stop("Cannot distinguish which is the chromosome column in the Marker Info.")
    chr = SNPInfo[[chr]]
    pos = grep("pos", names(SNPInfo), ignore.case=T)
    if (length(pos) != 1)
      stop("Cannot distinguish which is the position column in the Marker Info.")
    pos = SNPInfo[[pos]]
  } else {
    chr = rep(1, nrow(g))
    pos = as.integer(rep(0, nrow(g)))
  }
  # Create a gds file
  SNPRelate::snpgdsCreateGeno(filename, genmat=A, sample.id=colnames(g),
                              snp.id=1:nrow(g), snp.rs.id=rownames(g),
                              snp.chromosome=chr, snp.position=pos,
                              snpfirstdim=TRUE)
  return(filename)
}

##-----------------------------------------------------------------------------

#' Calculate kinship using \code{SNPRelate} package
#'
#' @param object \code{Genotype}. The object to calculate kinship from
#' @param method \code{character}. The method to use for identity-by-descent
#' calculation. Either "King" for using the "KING-robust" method, or "MoM" for
#' PLINK method of moment.
#' @param maf \code{numeric}. Cutoff for minor allele frequency.
#' @param missing.rate \code{numeric}. Cutoff for missing call rate.
#' @param rm.negatives \code{logical}. Should negative values be changed to 0?
#' @return \code{numeric}. the estimated kinship coefficients.
#'
setGeneric("kinship",
           function(object, method=c("King", "MoM"), maf=0.05,
                    missing.rate=0.05, rm.negatives=T) {
             standardGeneric("kinship")
             })
setMethod("kinship","Genotype",
          function(object, method=c("King", "MoM"), maf=0.05,
                   missing.rate=0.05, rm.negatives=T) {
  g = getCalls(object)
  A = matrix(data=rep(3, length(g)), nrow=nrow(g), ncol=ncol(g),
             dimnames=dimnames(g))
  A[g=="XX"] <- 0
  A[g=="XY"] <- 1
  A[g=="YY"] <- 2
  # Create a gds file
  SNPRelate::snpgdsCreateGeno("Genotype.gds", genmat=A, sample.id=colnames(g),
                   snp.id=1:nrow(g), snp.rs.id=rownames(g),
                   snp.chromosome=rep(1, nrow(g)),
                   snp.position=as.integer(rep(0, nrow(g))), snpfirstdim=TRUE)
  King_ibd = function(G, maf, missing.rate, rm.negatives) {
    ibd <- SNPRelate::snpgdsIBDKING(G, maf=maf, missing.rate=missing.rate,
                                    type="KING-robust")
    M = ibd$kinship
    if (rm.negatives) M[M<0] = 0
    return(M)
  }
  MoM_ibd = function(G, maf, missing.rate) {
    ibd <- SNPRelate::snpgdsIBDMoM(G, maf=maf, missing.rate=missing.rate)
    ibd.coeff <- SNPRelate::snpgdsIBDSelection(ibd)
    M = matrix(0, nrow=ncol(g), ncol=ncol(g), dimnames=list(colnames(g),
                                                            colnames(g)))
    for (s in unique(ibd.coeff$ID1)) {
      M[s, ibd.coeff$ID2[ibd.coeff$ID1==s]] <- 1 - ibd.coeff$kinship[ibd.coeff$ID1==s]
      M[ibd.coeff$ID2[ibd.coeff$ID1==s], s] <- 1 - ibd.coeff$kinship[ibd.coeff$ID1==s]
    }
    return(M)
  }
  G <- SNPRelate::snpgdsOpen("Genotype.gds")
  K = switch(method,
             King = King_ibd(G, maf=maf, missing.rate=missing.rate,
                             rm.negatives=rm.negatives),
             MoM = MoM_ibd(G, maf=maf, missing.rate=missing.rate)
  )
  SNPRelate::snpgdsClose(G)
  file.remove("Genotype.gds")
  colnames(K) <- rownames(K) <- SampleNames(MxS)
  return(K)
})

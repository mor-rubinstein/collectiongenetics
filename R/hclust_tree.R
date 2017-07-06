setClass("phyloS4",
         representation(edge = "matrix",
                        Nnode = "integer",
                        tip.label = "character",
                        edge.length = "numeric"))
setOldClass("phylo", S4Class="phyloS4")
selectMethod("show", "phylo")
removeClass("phyloS4")



#' Genetic distance between samples
#'
#' @param object \code{Genotype}. The object to calculate the distance in.
#' @param method \code{character}. The method of distance. One of the
#' following:
#' #' \describe{
#'   \item{"euclidean"}{euclidean distance between 2 vectors of genotype calls,
#'    where XX is translated to 1, XY to 2, and YY to 3.}
#'   \item{"PSA"}{proportion of shared alleles, and}
#'   \item{"IBD"}{kinship, where kinship is calculated using \code{SNPRelate}
#'   package \code{snpgdsIBDKING} function.}
#' }
#' @param maf \code{numeric}. Minor allele frequency threshold for
#' \code{snpgdsIBDKING}.
#' @param missing.rate \code{numeric}. Use SNPs with "<= missing.rate" only in
#' \code{snpgdsIBDKING}.
#'
#' @return \code{dist}
#'
#' @details This function is used by \code{\link{bootstrap}}, which uses
#' genetic distance to build an hierarchical clustered tree.
#'
#' @examples
#' data("citrus_clean")
#' De = GeneticDistance(MxS, method="euclidean")
#' Dp = GeneticDistance(MxS, method="PSA")
#' Di = GeneticDistance(MxS, method="IBD")

setGeneric("GeneticDistance",
           function(object, method=c("euclidean","PSA", "IBD"), maf=0.05,
                    missing.rate=0.05) {
             standardGeneric("GeneticDistance")
             })
setMethod("GeneticDistance","Genotype",
          function(object, method=c("euclidean","PSA", "IBD"), maf=0.05,
                   missing.rate=0.05) {
  PSA.dist<-function(g){
    M = matrix(numeric(), ncol(g), ncol(g),
               dimnames=list(colnames(g), colnames(g)))
    g[g=="XX"]=0
    g[g=="XY"]=0.5
    g[g=="YX"]=-0.5
    g[g=="YY"]=1
    mode(g)="numeric"
    for(i in 1:(ncol(g)-1)){
      for (j in (i+1):ncol(g)){
        PSA = 1-abs(g[,i]-g[,j])
        M[j,i] = 1-mean(ifelse(PSA==1.5,0.5,PSA),na.rm=T)
      }
    }
    return(as.dist(M))
  }
  IBD.dist = function(g, maf=0.05, missing.rate=0.05) {
    A = matrix(data=rep(3, length(g)), nrow=nrow(g), ncol=ncol(g), dimnames=dimnames(g))
    A[g=="XX"] <- 0
    A[g=="XY"] <- 1
    A[g=="YY"] <- 2
    # Create a gds file
    SNPRelate::snpgdsCreateGeno("Genotype.gds", genmat=A, sample.id=colnames(g),
                     snp.id=1:nrow(g), snp.rs.id=rownames(g), snp.chromosome=rep(1, nrow(g)),
                     snp.position=as.integer(rep(0, nrow(g))), snpfirstdim=TRUE)
    G <- SNPRelate::snpgdsOpen("Genotype.gds")
    ibd <- SNPRelate::snpgdsIBDKING(G, maf=maf, missing.rate=missing.rate, type="KING-robust")
    SNPRelate::snpgdsClose(G)
    M = ibd$kinship
    M[M<0] = 0
    M = 1-M
    return(as.dist(M))
  }
  D = switch(method,
    euclidean = dist(t(apply(getCalls(object),2,
                             function(x){
                               x[x=="XX"]=1
                               x[x=="XY"]=2
                               x[x=="YY"]=3;
                               as.numeric(x)}
                             ))),
    PSA = PSA.dist(getCalls(object)),
    IBD = IBD.dist(getCalls(object), maf=maf, missing.rate=missing.rate)
    )
  return(D)
})
##---------------------------------------------------------
#' creating bootstraped tree, based on hierarchical clustering
#'
#' \code{bootstrap} creats a 'tree' by clustering the samples according to
#' their similarity. It than calculates the reproducibility of each branch in
#' the tree by bootstrapping over the markers and building a tree with the
#' sampled set of markers.
#'
#' @param object \code{Genotype}. The object to build the tree from.
#' @param B \code{numeric}. The number of bootstraps.
#' @param treeM \code{character}. The method for calculating the tree.
#' \emph{"nj"} for using the \emph{ape} package neighbour joining function
#' \emph{nj}, or \emph{"upgma"} for using the \emph{stats} package
#' \emph{hclust} function, with 'average' agglomeration method.
#' @param distM \code{character}. The method for calculating the distance
#' matrix between samples. \emph{"PSA"} for proportion of shared alleles,
#' \emph{"euclidean"} for euclidean distance, where XX genotype is evaluated as
#' 1, XY as 2 and YY as 3, or \emph{"IBD"} for 1-kinship, calculated using the
#' \emph{"SNPRelate"} pacakge function \emph{"snpgdsIBDKING"}.
#' @param show.itr \code{logical}. Wheter to print each itteration.
#'
#' @return \code{phylo}. An \code{ape} package \code{phylo} object.
#'
#' @seealso \code{\link{invert_tree}} to change the leaf order of the tree.
#' \code{\link{GeneticDistance}}, \code{\link{getTreeLeafOrder}},
#' \code{\link{DrowPopulationTree}}
#'
#' @examples
#' data("citrus_clean")
#' b.tree = bootstrap(MxS, treeM="upgma")
#'
setGeneric("bootstrap",
           function(object, B=100, treeM="nj", distM="PSA", show.itr=FALSE) {
             standardGeneric("bootstrap")
             })
setMethod("bootstrap", "Genotype",
          function(object, B=100, treeM="nj", distM="PSA", show.itr=FALSE) {
  B.trees = list()
  length(B.trees) = B
  distMatrix = GeneticDistance(object, method=distM)
  if (treeM=="nj") {
    Tree = ape::nj(distMatrix)
    for (i in 1:B){
      f = sample(GetMarkerFilter(object),replace=TRUE)
      itt.object = SetMarkerFilter(object, Filter=f, print.n.markers=F)
      B.trees[[i]] = ape::nj(GeneticDistance(itt.object,method=distM))
      if (show.itr) cat("Tree No:",i,"\n")
    }
  }
  if (treeM=="upgma") {
    Tree = try(hclust(distMatrix, method="ward.D"), silent=T)
    if (class(Tree) == "try-error")
      Tree = hclust(distMatrix, method="ward")
    for (i in 1:B) {
      f = sample(GetMarkerFilter(object),replace=TRUE)
      itt.object = SetMarkerFilter(object,Filter=f, print.n.markers=F)
      x = try(hclust(GeneticDistance(itt.object, method=distM),method="ward.D"), silent=T)
      if (class(x) == "try-error")
        x = hclust(GeneticDistance(itt.object, method=distM),method="ward")
      B.trees[[i]] = x
      if (show.itr) cat("Tree No:",i,"\n")
    }

    ## reordering tree, so that nearby object pairs are adjacent
    oTree = gclus::reorder.hclust(Tree, distMatrix)
    o.B.trees = lapply(B.trees, gclus::reorder.hclust, dis=distMatrix)

    # changing from class 'hclust' to class 'phylo'
    oTree = ape::as.phylo(oTree)
    o.B.trees = lapply(o.B.trees, ape::as.phylo)
  }
  if (length(B.trees)==0) stop("No tree method was chosen")

  cat("Calculating boot values")
  oTree$node.label=ape::prop.clades(oTree, o.B.trees)
  Tree = oTree
  return(Tree)
})

##---------------------------------------------------------
#' Getting the order of leaves in a 'phylo' tree
#'
#' This allows you to re-order a Genotype object according to the similarity
#' between samples
#'
#' @param object \code{phylo}. The tree of interest
#' @return \code{integer}. The order of the leaves in the tree.
#'
#' @examples
#' data("Avocado_clean123")
#' library(ape)
#' data("Avocado_phylo")
#' o = getTreeLeafOrder(b.tree)
#' MxS = SetSampleOrder(MxS, Order=o)

setGeneric("getTreeLeafOrder", function(object) {
  standardGeneric("getTreeLeafOrder")
  })
setMethod("getTreeLeafOrder","ape::phylo",function(object) {
  x = object$edge
  x = x[x[, 2] %in% 1:length(object$tip.label), 2]
  x
})
##---------------------------------------------------------
#' Inverting the order of a 'phylo' tree or a branch of it
#'
#' This function allows you to invert the order of leaves in a 'phylo' object.
#' You can invert the whole tree, or just all the leaves that descend from a
#' certain node. It is sometimes desireable when plotting the tree.
#'
#' @param tree1 \code{phylo}. The tree to invert
#' @param node \code{numeric}. The node number around which to invert the tree.
#' If this argument is missing or NULL, the entire tree will be inverted.
#' @return \code{phylo}. The re-ordered tree.
#'
#' @seealso \code{\link{bootstrap}}, \code{\link{GeneticDistance}},
#' \code{\link{getTreeLeafOrder}}, \code{\link{DrowPopulationTree}}
#'
#' @examples
#' library("ape")
#' data("Avocado_phylo")
#' c.tree = invert_tree(b.tree)
#' d.tree = invert_tree(b.tree, 128)
#' plot(b.tree)
#' plot(c.tree)
#' plot(d.tree)

invert_tree = function(tree1, node=NULL) {
  rotateNodes = function(node, edge, new_edge) {
    x = edge[, 1]==node
    if (sum(x)) {
      l = min(which(new_edge[, 1]==0))
      ch = edge[x, 2]
      new_edge[l, ] = c(node, ch[2])
      new_edge = rotateNodes(ch[2], edge, new_edge)
      l = min(which(new_edge[, 1]==0))
      new_edge[l, ] = c(node, ch[1])
      new_edge = rotateNodes(ch[1], edge, new_edge)
    }
    return(new_edge)
  }
  edge <- new_edge <- tree1$edge
  if (is.null(node)) {
    new_edge[edge!=0] <-0
    node = edge[1,1]
  } else {
    E = which(edge[, 1] == node)
    if (length(E)==0) {
      cat("No such node in the tree")
      return(tree1)
    } else {
      E1 = which(edge[, 1] %in% edge[E, 2])
      while (!all(E1 %in% E)) {
        E = sort(union(E, E1))
        E1 = which(edge[, 1] %in% edge[E, 2])
      }
    }
    new_edge[E,] <-0
  }
  new_edge = rotateNodes(node, edge, new_edge)
  tree2 = tree1
  x = max(edge[, 2])
  x1 = edge[, 1] + (edge[, 2])/(x+1)
  x2 = new_edge[, 1] + (new_edge[, 2])/(x+1)
  o = order(x1)[rank(x2)]
  tree2$edge.length = tree1$edge.length[o]
  tree2$edge = new_edge
  return(tree2)
}

##---------------------------------------------------------
#' Drows a population \code{phylo} 'tree' with sample info specifications
#'
#' This function allows you to set the leave names and to colour the leaves
#' according to any category in the \code{SampleInfo}.
#'
#' @param object \code{phylo}. The tree to plot
#' @param SampleInfo \code{SampleInfo}. A table of sample information. The row
#' names of the table should be identical to the \code{object$tip.label}.
#' @param tip.label.by \code{character}. Which column in SampleInfo to colour
#' the tree tip labels by.
#' @param colourby \code{character}. Colour tip lables according to this column
#' in \code{SampleInfo}.
#' @param ColourPalette \code{character}. ColourPalette for tip lables.
#' @param legend \code{logical}. Should a legend be added?
#' @param legend.x \code{numeric}. legend x position
#' @param legend.y \code{numeric}. legend y position
#' @param legend.cex \code{numeric}. legend font size
#' @param mark.node.labels \code{logical}. Should the nodes be labeled?
#'
#' @return \code{NULL}.
#' @seealso \code{\link{bootstrap}} to create the tree, \code{\link{invert_tree}} to
#'  change the leaf order of the tree and  \code{\link{getTreeLeafOrder}} to get
#'  the tree leaf order.
#' @examples
#' data("Avocado_phylo")
#' data("Avocado_clean123")
#' SampleInfo = getSampleInfo(MxS)
#' DrowPopulationTree(b.tree, SampleInfo, colourby="Category", legend=TRUE)

setGeneric("DrowPopulationTree",
           function(object, SampleInfo=NULL, tip.label.by="Cultivar",
                    colourby=NULL, ColourPalette=NULL, legend=TRUE, legend.x=NULL,
                    legend.y=NULL, legend.cex=1, mark.node.labels=FALSE) {
             standardGeneric("DrowPopulationTree")
             })

setMethod("DrowPopulationTree", "ape::phylo",
          function(object, SampleInfo=NULL, tip.label.by="Cultivar",
                   colourby=NULL, ColourPalette=NULL, legend=TRUE, legend.x=NULL,
                   legend.y=NULL, legend.cex=1, mark.node.labels=FALSE) {
  if (!is.null(SampleInfo) && !(is.null(colourby))) {
    if (!(colourby %in% names(SampleInfo)))
      stop(paste(colourby, "is not in the names of SampleInfo."))
    category = SampleInfo[[colourby]]
    category.name = unique(category)
    category.name = sort(category.name)
    grey = F
    if ("" %in% category.name) {
      category.name = category.name[order(nchar(category.name), decreasing=T)]
      grey = T
    }
    x = grep("Miscellaneous", category.name, ignore.case=T)
    if (length(x)) {
      category.name = c(category.name[-x], category.name[x])
      grey=T
    }
    x = grep("unknown", category.name, ignore.case=T)
    if (length(x)) {
      category.name = c(category.name[-x], category.name[x])
      grey=T
    }
    x = grep("not ", category.name, ignore.case=T)
    if (length(x)) {
      category.name = c(category.name[-x], category.name[x])
      grey=T
    }
    categoryColour = character(length(category))
    if (is.null(ColourPalette)) {
      Set = ""
      if (length(category.name)<=8) Set = "Dark2"
      if (length(category.name)<=10) Set = "Set1" else
        if (length(category.name)<=13) Set = "Set3"
        if (nchar(Set)) {
          if (grey) ColourPalette = c(RColorBrewer::brewer.pal(length(category.name)-1, Set), "grey20") else
            ColourPalette = RColorBrewer::brewer.pal(length(category.name), Set)
        } else
          stop("Too many categories to find a suitable colour palette.  Please provide a colour palette.")
    }
    for (i in 1:length(category.name)) {
      categoryColour[category==category.name[i]] = ColourPalette[i]
    }
  } else categoryColour="black"
  if (!is.null(SampleInfo) && !(is.null(tip.label.by))) {
    if (!(tip.label.by %in% names(SampleInfo)))
      stop(paste(tip.label.by, "is not in the names of SampleInfo."))
    object$tip.label = SampleInfo[object$tip.label, tip.label.by]
  }
  op = par(mar=c(2, 2, 2, 3)+0.1)

  population.tree = ape::plot.phylo(object,direction="downwards", show.tip.label = TRUE,cex=0.5, edge.color = "black",tip.color = categoryColour)
  if (legend & !is.null(colourby)) {
    if (is.null(legend.x)) legend.x = (population.tree$x.lim[2]-population.tree$x.lim[1])*4/5
    if (is.null(legend.y)) legend.y = (population.tree$y.lim[2]-population.tree$y.lim[1])
    legend(legend.x, legend.y, legend=category.name, text.col=ColourPalette, cex=legend.cex)
  }
  if (mark.node.labels) {
    node.labels = object$node.label
    node.labels[node.labels<75] = NA
    ape::nodelabels(node.labels,adj=c(+0.5,-0.5),frame="n",bg="white", cex=0.6)
  }
  par(op)
  x = NULL
})

##---------------------------------------------------------
#' Drows a population \code{phylo} 'tree' with STRUCTRE results at the bottom
#'
#' This function creates a tiff image file, with a dendogram of the samples,
#' and STRUCTURE results below it.
#'
#' @param tree \code{phylo}. The tree to plot
#' @param anc \code{data.frame}. A table with the summary of STRUCTURE results.
#' For STRUCURE with K populations, this table should contain K columns, named
#' "cluster_1", "cluster_2", ... "cluster_K", specifying the relative
#' ancestries of each individual. The number of rows should be identical to the
#' number of individuals, and it has to be ordered according to the order of
#' tips in the tree.
#' @param addlegend \code{logical}. Should a legend be added?
#' @param parameter_list \code{list}. A named list of additional arguments to
#' pass to the function. possible arguments are:
#' \describe{
#'   \item{"plotname"}{\code{character}. The tiff image file name. default is
#' "cluster_STRUCTURE".}
#'   \item{"plotname_prefix"}{\code{character}. A prefix for the tiff image
#'   file name. default is NULL.}
#'   \item{"BarplotColour"}{The colours to use for the STRUCTURE results bar
#'   plot.}
#'   \item{"SampleInfo"}{\code{SampleInfo}. A table of sample information. The
#'    row names of the table should be identical to the \code{object$tip.label}
#'    .}
#'   \item{"colourby"}{\code{character}. Colour the tree tip lables according
#'   to this column in \code{SampleInfo}.}
#'   \item{"ColourPalette"}{\code{character}. ColourPalette for the tree tip
#' lables.}
#'   \item{"tip.label.by"}{\code{character}. Which column in SampleInfo to
#'   colour the tree tip labels by.}
#'   \item{"legend.x"}{\code{numeric}. legend x position}
#'   \item{"legend.y"}{\code{numeric}. legend y position}
#'   \item{"legend.cex"}{\code{numeric}. legend font size}
#'   \item{"mark.node.labels"}{\code{logical}. Should the nodes be labeled?}
#'   \item{"node.label.cutoff"}{\code{numeric}. What is the minimum value of
#'   node label to include in the tree graph? The node labels numeric values
#'   specifying the reproducibility of each node.}
#' }
#'
#' @return \code{NULL}.
#' @seealso \code{\link{bootstrap}} to create the tree, \code{\link{invert_tree}}
#' to change the leaf order of the tree, \code{\link{DrowPopulationTree}} to
#' plot only the tree, and \code{\link{getTreeLeafOrder}} to get the tree leaf
#' order.
#' @examples
#' data("Avocado_phylo")
#' data("Avocado_clean123")
#' data("Avocado_ancestry")
#' SampleInfo = getSampleInfo(MxS)
#' cluster_STRUCTURE_plot(b.tree, anc=ancestry, parameter_list=list(SampleInfo=
#' getSampleInfo(MxS), colourby="Category", mark.node.labels=T))
cluster_STRUCTURE_plot = function(tree, anc, addlegend=T,
                                  parameter_list=list()) {
  if ("plotname" %in% names(parameter_list)) plotname = parameter_list$plotname else
    plotname="cluster_STRUCTURE"
  plotname = paste(plotname, "tiff", sep=".")
  if ("plotname_prefix" %in% names(parameter_list))
    plotname_prefix = parameter_list$plotname_prefix else plotname_prefix=NULL
  if (length(plotname_prefix)) plotname = paste(plotname_prefix, plotname, sep="_")
  x = anc[grep("cluster", names(anc))]
  x = t(x)
  if ("Cultivar" %in% names(anc)) colnames(x) = anc$Cultivar else
    if ("cultivar" %in% names(anc)) colnames(x) = anc$Cultivar else {
      name = grep("name", names(anc), ignore.case=T)
      if (length(name)==1) colnames(x) = anc[[name]] else {
        sam = grep("sample", names(anc), ignore.case=T)
        if (length(intersect(sam, name))==1)
          colnames(x) = anc[[intersect(sam, name)]] else {
            if (length(sam)==1) colnames(x) = anc[[sam]] else
              colnames(x) <- NULL
          }
      }
    }
  tiff(plotname, width=13, height=8, units="in", res=600)
  split.screen(c(2,1))
  screen(2)
  op = par(mar=c(7,2,0,2), oma=c(0,1,0,1))
  if ("BarplotColour" %in% names(parameter_list)) colour=parameter_list$BarplotColour else colour=2:(nrow(x)+1)
  barplot(x, col=colour,las=2, yaxt="n", border=NA, space=0, cex.names=0.7)
  if ("SampleInfo" %in% names(parameter_list)) {
    SampleInfo = parameter_list$SampleInfo
    if ("colourby" %in% names(parameter_list)) {
      colourby = parameter_list$colourby
      if (!(colourby %in% names(SampleInfo)))
        stop(paste(colourby, "is not in the names of SampleInfo."))
      category = SampleInfo[tree$tip.label, colourby]
      category.name = unique(category)
      category.name = sort(category.name)
      grey = F
      if ("" %in% category.name) {
        category.name = category.name[order(nchar(category.name), decreasing=T)]
        grey = T
      }
      x = grep("Miscellaneous", category.name, ignore.case=T)
      if (length(x)) {
        category.name = c(category.name[-x], category.name[x])
        grey=T
      }
      x = grep("unknown", category.name, ignore.case=T)
      if (length(x)) {
        category.name = c(category.name[-x], category.name[x])
        grey=T
      }
      x = grep("not ", category.name, ignore.case=T)
      if (length(x)) {
        category.name = c(category.name[-x], category.name[x])
        grey=T
      }
      categoryColour = character(length(category))
      if ("ColourPalette" %in% names(parameter_list)) {
        ColourPalette = parameter_list$ColourPalette
        if (all(names(ColourPalette) %in% category.name))
          category.name = c(names(ColourPalette), setdiff(category.name, names(ColourPalette))) else
            if (all(category.name %in% names(ColourPalette)))
              category.name = names(ColourPalette)[ColourPalette %in% category.name]
            for (i in category.name)
              categoryColour[category == i] = ColourPalette[i]
            categoryColour[!(category) %in% category.name] <- "black"
      } else {
        Set = ""
        if (length(category.name)<=8) Set = "Dark2"
        if (length(category.name)<=10) Set = "Set1" else
          if (length(category.name)<=13) Set = "Set3"
          if (nchar(Set)) {
            if (grey) ColourPalette = c(RColorBrewer::brewer.pal(length(category.name)-1, Set), "grey20") else
              ColourPalette = RColorBrewer::brewer.pal(length(category.name), Set)
          } else
            stop("Too many categories to find a suitable colour palette.  Please provide a colour palette.")
          for (i in 1:length(category.name)) {
            categoryColour[category==category.name[i]] = ColourPalette[i]
          }
      }
    } else categoryColour="black"
    if ("tip.label.by" %in% names(parameter_list)) {
      tip.label.by = parameter_list$tip.label.by
      if (!(tip.label.by %in% names(SampleInfo)))
        stop(paste(tip.label.by, "is not in the names of SampleInfo."))
      tree$tip.label = SampleInfo[tree$tip.label, tip.label.by]
    }
  } else categoryColour="black"
  screen(1)
  par(mar=c(0,2,2,2), oma=c(0,1,0,1))
  population.tree = ape::plot.phylo(tree, direction = "downwards",
                                    show.tip.label = TRUE, cex = 0.5,
                                    edge.color = "black",
                                    tip.color = categoryColour)
  if (addlegend & "colourby" %in% names(parameter_list)) {
    if ("legend.x" %in% names(parameter_list)) legend.x=parameter_list$legend.x else
      legend.x = (population.tree$x.lim[2]-population.tree$x.lim[1])*4/5
    if ("legend.y" %in% names(parameter_list)) legend.y=parameter_list$legend.y else
      legend.y = (population.tree$y.lim[2]-population.tree$y.lim[1])*4/5
    if ("legend.cex" %in% names(parameter_list)) legend.cex=parameter_list$legend.cex else
      legend.cex = 1
    category.name = category.name[nchar(category.name)>0]
    if (all(names(ColourPalette) %in% category.name))
      category.name = c(names(ColourPalette), setdiff(category.name, names(ColourPalette))) else
        if (all(category.name %in% names(ColourPalette)))
          category.name = names(ColourPalette)[ColourPalette %in% category.name]
    legend(legend.x, legend.y, legend=category.name, text.col=ColourPalette, cex=legend.cex)
  }
  if ("mark.node.labels" %in% names(parameter_list)) {
    if (parameter_list$mark.node.labels) {
      node.labels = tree$node.label
      if ("node.label.cutoff" %in% names(parameter_list))
        node.labels[node.labels<parameter_list$node.label.cutoff] = NA
      ape::nodelabels(node.labels,adj=c(+0.5,-0.5),frame="n",bg="white", cex=0.6)

    }
  }
  par(op)
  close.screen(all = TRUE)
  dev.off()
  cat("created a new plot file \'", plotname, "\'.\n", sep="")
}

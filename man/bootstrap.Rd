% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hclust_tree.R
\name{bootstrap}
\alias{bootstrap}
\title{creating bootstraped tree, based on hierarchical clustering}
\usage{
bootstrap(object, B = 100, treeM = "nj", distM = "PSA",
  show.itr = FALSE)
}
\arguments{
\item{object}{\code{Genotype}. The object to build the tree from.}

\item{B}{\code{numeric}. The number of bootstraps.}

\item{treeM}{\code{character}. The method for calculating the tree.
\emph{"nj"} for using the \emph{ape} package neighbour joining function
\emph{nj}, or \emph{"upgma"} for using the \emph{stats} package
\emph{hclust} function, with 'average' agglomeration method.}

\item{distM}{\code{character}. The method for calculating the distance
matrix between samples. \emph{"PSA"} for proportion of shared alleles,
\emph{"euclidean"} for euclidean distance, where XX genotype is evaluated as
1, XY as 2 and YY as 3, or \emph{"IBD"} for 1-kinship, calculated using the
\emph{"SNPRelate"} pacakge function \emph{"snpgdsIBDKING"}.}

\item{show.itr}{\code{logical}. Wheter to print each itteration.}
}
\value{
\code{phylo}. An \code{ape} package \code{phylo} object.
}
\description{
\code{bootstrap} creats a 'tree' by clustering the samples according to
their similarity. It than calculates the reproducibility of each branch in
the tree by bootstrapping over the markers and building a tree with the
sampled set of markers.
}
\examples{
data("citrus_clean")
b.tree = bootstrap(MxS, treeM="upgma")

}
\seealso{
\code{\link{invert_tree}} to change the leaf order of the tree.
\code{\link{GeneticDistance}}, \code{\link{getTreeLeafOrder}},
\code{\link{DrowPopulationTree}}
}


#' Get orthologous cluster ID
#' 
#' Find orthologous cluster given gene/protein IDs
#' 
#' @param query.id query IDs
#' @param id.type type.id which type of IDs are used: Entrez Gene ID, Ensembl gene ID, Ensembl peptide ID
#' @param taxid query species
#' @return a vector of orthologous cluster IDs (NAs for failed query IDs)
#' @keywords ortholog gortholog
#' @export
#' @examples
#' GetClusterID(query.id=c("YIL106W", "YOR358W"), id.type="ensembl_pep_id")
#' GetClusterID(query.id=c(854700, 854540), id.type="entrez_gene_id")
#' @author Mark J Chen \email{chenj99@@gene.com}
#' 
GetClusterID <- function(query.id=query.id, id.type=c("entrez_gene_id", "ensembl_pep_id", "ensembl_gene_id"), 
                         taxid=taxid, clusterid.cidx=1, taxid.cidx=2, id.cidx=3,
                         no.id=NA, diagnose=FALSE) {
  # no input?
  empty <- ( is.na(query.id) | query.id =='' )
  if ( sum(empty) > 0 ) warning('One or more empty gene ID/cluster in input')  
  # find cluster of query 
  gortholog.data <- LoadGorthologData(id.type)
  query.ridx <- match(query.id, gortholog.data[ , id.cidx])
  cluster.id <- gortholog.data[query.ridx, clusterid.cidx]
  # fail to map to cluster?
  noentry <- (is.na(query.ridx) & !empty)
  if ( sum(noentry) > 0 ) warning('One or more gene input gene ID/cluster not found in gortholog')
  return( as.vector(cluster.id) )
}


#' Get alignment
#' 
#' Given a query, fetch the alignment of its orthologous cluster from gOrtholog website
#' 
#' @param query.id query IDs
#' @param id.type type.id which type of IDs are used: Entrez Gene ID, Ensembl gene ID, Ensembl peptide ID
#' @param taxid query species
#' @return an object of class AAMultipleAlignment (see \code{\link{MultipleAlignment-class}})
#' @keywords ortholog gortholog
#' @import Biostrings
#' @export
#' @seealso \code{\link{MultipleAlignment-class}}
#' @examples
#' ptpn4a1.aln <- GetAlignment(query.id="ENSP00000359685", id.type="ensembl_pep_id")
#' @author Mark J Chen \email{chenj99@@gene.com}
#' 
GetAlignment <- function(query.id=query.id, id.type=c("entrez_gene_id", "ensembl_pep_id", "ensembl_gene_id"), 
                         taxid=taxid, clusterid.cidx=1, taxid.cidx=2, id.cidx=3,
                         no.id=NA, diagnose=FALSE, 
                         url.base="http://resdev.gene.com/gOrtholog/") {
  # no input?
  empty <- ( is.na(query.id) | query.id =='' )
  if ( sum(empty) > 0 ) warning('No query gene')
  # more input?
  if ( sum(empty) > 1 ) warning('More than one query genes; Fetch the first.')
  # find cluster of query 
  gortholog.data <- LoadGorthologData(id.type)
  query.ridx <- match(query.id, gortholog.data[ , id.cidx])
  cluster.id <- gortholog.data[query.ridx, clusterid.cidx]
  # fail to map to cluster?
  noentry <- (is.na(query.ridx) & !empty)
  if ( sum(noentry) > 0 ) warning('Query is not found in gortholog')
  # fetch cluster
  url.fasta <- paste(url.base, 'view/cluster/', cluster.id[1], '/sequence/alignment/muscle', sep="")
  library(Biostrings)
  alignment <- readAAMultipleAlignment(url.fasta, "fasta")
  return(alignment)
}


#' Get tree
#' 
#' Given a query, fetch the tree of it and its orthologous cluster from 
#' gOrtholog website
#' 
#' @param query.id query IDs
#' @param id.type type.id which type of IDs are used: Entrez Gene ID, Ensembl gene ID, Ensembl peptide ID
#' @param taxid query species
#' @return An object of class "phylo" (see \code{\link{read.tree}}).
#' @keywords ortholog gortholog
#' @import ape
#' @export
#' @seealso \code{\link{read.tree}}
#' @examples
#' ptpn4a1.aln <- GetAlignment(query.id="ENSP00000359685", id.type="ensembl_pep_id")
#' @author Mark J Chen \email{chenj99@@gene.com}
#' 
GetTree <- function(query.id=query.id, id.type=c("entrez_gene_id", "ensembl_pep_id", "ensembl_gene_id"), 
                         taxid=taxid, clusterid.cidx=1, taxid.cidx=2, id.cidx=3,
                         no.id=NA, diagnose=FALSE, 
                         url.base="http://resdev.gene.com/gOrtholog/") {
  # no input?
  empty <- ( is.na(query.id) | query.id =='' )
  if ( sum(empty) > 0 ) warning('No query gene')
  # more input?
  if ( sum(empty) > 1 ) warning('More than one query genes; Fetch the first.')
  # find cluster of query 
  gortholog.data <- LoadGorthologData(id.type)
  query.ridx <- match(query.id, gortholog.data[ , id.cidx])
  cluster.id <- gortholog.data[query.ridx, clusterid.cidx]
  # fail to map to cluster?
  noentry <- (is.na(query.ridx) & !empty)
  if ( sum(noentry) > 0 ) warning('Query is not found in gortholog')
  # fetch cluster
  url.ninja.newick <- paste(url.base, 'view/cluster/', cluster.id[1], '/sequence/tree/ninja/newick', sep="")
  library(ape)
  tre <- read.tree(url.ninja.newick)
  return(tre)
}
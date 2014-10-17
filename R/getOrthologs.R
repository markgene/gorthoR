#' Get all orthologs in given species
#' 
#' Get all orthologs in given species by taxonomy IDs.
#' 
#' @param taxid NCBI taxonomy ID of given species
#' @param type.id which type of IDs are used: Entrez Gene ID, Ensembl gene ID, Ensembl peptide ID
#' @return data frame of three columns: cluster ID, taxonomy ID, Ensembl peptide ID.
#' @keywords ortholog gortholog
#' @export
#' @examples
#' human.mouse <- GetAllOrthologs(taxid=c(9606, 10090), id.type="ensembl_pep_id")
#' scer.pombe <- GetAllOrthologs(taxid=c(4932, 4896), id.type="entrez_gene_id")
#' @author Mark J Chen \email{chenj99@@gene.com}
#' 
GetAllOrthologs <- function(taxid=taxid, id.type=c("entrez_gene_id", "ensembl_pep_id", "ensembl_gene_id"), taxid.ridx=2) {
  gortholog.data <- LoadGorthologData(id.type)
  #row.index <- match(taxid, gortholog.data[ , taxid.index])
  row.index <- unlist(
    sapply(taxid, function(x, gd=gd, idx=idx) { return( which(gd[, idx] == x )) }, 
           gd=gortholog.data, idx=taxid.ridx) )
  return(gortholog.data[row.index,])
}


#' Find orthologs in other species by gene/protein IDs
#' 
#' Find orthologs in other species given IDs and their source species taxonomy ID
#'
#' @param query.id query IDs
#' @param id.type which type of IDs are used: Entrez Gene ID, Ensembl gene ID, Ensembl peptide ID
#' @param source.taxid NCBI taxonomy ID of source species
#' @param target.taxid NCBI taxonomy ID of target species
#' @return list of elements each of which is a vector of Ensembl peptide IDs in target species
#' @keywords ortholog gortholog
#' @export
#' @examples
#' GetOrthologs(source.taxid=4932, target.taxid=4896, query.id=c("YIL106W", "YOR358W"), id.type="ensembl_pep_id")
#' GetOrthologs(source.taxid=4932, target.taxid=4896, query.id=c(854700, 854540), id.type="entrez_gene_id")
#' @author Mark J Chen \email{chenj99@@gene.com}
#' 
GetOrthologs <- function(source.taxid=source.taxid, target.taxid=target.taxid, 
                         query.id=query.id, id.type=c("entrez_gene_id", "ensembl_pep_id", "ensembl_gene_id"), 
                         clusterid.cidx=1, taxid.cidx=2, id.cidx=3,
                         no.id=NA, diagnose=FALSE) {
  # no input?
  empty <- ( is.na(query.id) | query.id =='' )
  if ( sum(empty) > 0 ) warning('One or more empty gene ID in input')	
  # find cluster of query 
  gortholog.data <- LoadGorthologData(id.type)
  query.ridx <- match(query.id, gortholog.data[ , id.cidx])
  cluster.id <- gortholog.data[query.ridx, clusterid.cidx]
  # fail to map to cluster?
  noentry <- (is.na(query.ridx) & !empty)
  if ( sum(noentry) > 0 ) warning('One or more gene input gene ID not found in gortholog')
  # create list since one-to-many
  target.id <- as.list( rep(no.id, length(query.id)) )
  # target species row index
  target.taxid.ridx <- gortholog.data[ , taxid.cidx] == target.taxid
  # which row to process
  idx <- which(!empty & !noentry)
  if (sum(!empty & !noentry) > 0) {
    for (i in 1:sum(!empty & !noentry)) {
      # For each row of query and its cluster, select target species AND have the cluster ID; return ensembl pep id column
      target.id[[idx[i]]] <- gortholog.data[ which( target.taxid.ridx & gortholog.data[ , clusterid.cidx] == cluster.id[idx[i]]), id.cidx]
    }
  }
  # map to no gene (i.e. though query map to a cluster, but the cluster is absent in target species)
  target.id.nb<-sapply(target.id, function(x) { length(x) })
  no.target.id <- target.id.nb == 0
  if (sum(no.target.id) > 0) warning('One or more gene ID with no target provided in homologue table')
  # assign to NA
  target.id[noentry | no.target.id] <- no.id
  # 
  if (diagnose) { 
    return(list(target.id,empty,noentry,no.target.id))
  } else {
    return(target.id)
  }
}

#' Load gOrtholog data
#' 
#' Load gOrtholog data
#' 
#' @param type.id which type of IDs are used: Entrez Gene ID, Ensembl gene ID, Ensembl peptide ID
#' @return loaded data
#' @export
#' @examples
#' LoadGorthologData()
#' @author Mark J Chen \email{chenj99@@gene.com}
#'
LoadGorthologData <- function(id.type=id.type) {
  data(gorthologData)
  if (id.type == "entrez_gene_id") {
    data(gorthologDataEntrezGeneid)  # where to load this in package-wide scope?
    return(gortholog.entrez.geneid)
  } else if (id.type == "ensembl_gene_id") {
    data(gorthologDataEnsemblGeneid)  # where to load this in package-wide scope?
    return(gortholog.ensembl.geneid)
  } else if (id.type == "ensembl_pep_id") {
    data(gorthologDataEnsemblPepid)  # where to load this in package-wide scope?
    return(gortholog.ensembl.pepid)
  }
}


# # ID converter
# IdConverter <- function(id.type=id.type) {
# 
# }

# # Generate R data file
# GenerateRdataFile <- function() {
#   # gortholog.data
#   gortholog.data <- read.table("extdata//cid_taxid_pgi.tab", header=FALSE, stringsAsFactors=FALSE, sep="\t")
#   colnames(gortholog.data) <- c("cluster_id", "taxonomy_id", "primary_id")
#   save("gortholog.data", file="data/gorthologData.rda")
#   # gortholog.entrez.geneid
#   gortholog.entrez.geneid <- read.table("extdata//cid_taxid_gene_id.tab", header=FALSE, stringsAsFactors=FALSE, sep="\t")
#   colnames(gortholog.entrez.geneid) <- c("cluster_id", "taxonomy_id", "entrez_geneid")
#   save("gortholog.entrez.geneid", file="data/gorthologDataEntrezGeneid.rda")
#   # gortholog.ensembl.geneid
#   gortholog.ensembl.geneid <- read.table("extdata//cid_taxid_ens_gene_id.tab", header=FALSE, stringsAsFactors=FALSE, sep="\t")
#   colnames(gortholog.ensembl.geneid) <- c("cluster_id", "taxonomy_id", "ensembl_geneid")
#   save("gortholog.ensembl.geneid", file="data/gorthologDataEnsemblGeneid.rda")
#   # gortholog.ensembl.peptide.id
#   gortholog.ensembl.pepid <- read.table("extdata//cid_taxid_ens_pep_id.tab", header=FALSE, stringsAsFactors=FALSE, sep="\t")
#   colnames(gortholog.ensembl.pepid) <- c("cluster_id", "taxonomy_id", "ensembl_pepid")
#   save("gortholog.ensembl.pepid", file="data/gorthologDataEnsemblPepid.rda")
#   # ensembl peptide id
#   ensembl.pepid <- read.table("extdata//pgi_ens_pep_id.tab", header=FALSE, stringsAsFactors=FALSE, sep="\t")
#   colnames(ensembl.pepid) <- c("primary_id", "ensembl_pepid")
#   save("ensembl.pepid", file="data/ensemblPepid.rda")
#   # ensembl gene id
#   ensembl.geneid <- read.table("extdata//pgi_ens_gene_id.tab", header=FALSE, stringsAsFactors=FALSE, sep="\t")
#   colnames(ensembl.geneid) <- c("primary_id", "ensembl_geneid")
#   save("ensembl.geneid", file="data/ensemblGeneid.rda")
#   # entrez gene
#   entrez.geneid <- read.table("extdata//pgi_gene_id.tab", header=FALSE, stringsAsFactors=FALSE, sep="\t")
#   colnames(entrez.geneid) <- c("primary_id", "entrez_geneid")
#   save("entrez.geneid", file="data/entrezGeneid.rda")
# }
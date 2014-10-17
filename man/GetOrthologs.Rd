\name{GetOrthologs}
\alias{GetOrthologs}
\title{Find orthologs in other species by gene/protein IDs}
\usage{
GetOrthologs(source.taxid = source.taxid, target.taxid = target.taxid,
  query.id = query.id, id.type = c("entrez_gene_id", "ensembl_pep_id",
  "ensembl_gene_id"), clusterid.cidx = 1, taxid.cidx = 2, id.cidx = 3,
  no.id = NA, diagnose = FALSE)
}
\arguments{
  \item{query.id}{query IDs}

  \item{id.type}{which type of IDs are used: Entrez Gene
  ID, Ensembl gene ID, Ensembl peptide ID}

  \item{source.taxid}{NCBI taxonomy ID of source species}

  \item{target.taxid}{NCBI taxonomy ID of target species}
}
\value{
list of elements each of which is a vector of Ensembl
peptide IDs in target species
}
\description{
Find orthologs in other species given IDs and their source
species taxonomy ID
}
\examples{
GetOrthologs(source.taxid=4932, target.taxid=4896, query.id=c("YIL106W", "YOR358W"), id.type="ensembl_pep_id")
GetOrthologs(source.taxid=4932, target.taxid=4896, query.id=c(854700, 854540), id.type="entrez_gene_id")
}
\author{
Mark J Chen \email{chenj99@gene.com}
}
\keyword{gortholog}
\keyword{ortholog}


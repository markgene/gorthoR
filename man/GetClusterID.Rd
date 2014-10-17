\name{GetClusterID}
\alias{GetClusterID}
\title{Get orthologous cluster ID}
\usage{
GetClusterID(query.id = query.id, id.type = c("entrez_gene_id",
  "ensembl_pep_id", "ensembl_gene_id"), taxid = taxid, clusterid.cidx = 1,
  taxid.cidx = 2, id.cidx = 3, no.id = NA, diagnose = FALSE)
}
\arguments{
  \item{query.id}{query IDs}

  \item{id.type}{type.id which type of IDs are used: Entrez
  Gene ID, Ensembl gene ID, Ensembl peptide ID}

  \item{taxid}{query species}
}
\value{
a vector of orthologous cluster IDs (NAs for failed query
IDs)
}
\description{
Find orthologous cluster given gene/protein IDs
}
\examples{
GetClusterID(query.id=c("YIL106W", "YOR358W"), id.type="ensembl_pep_id")
GetClusterID(query.id=c(854700, 854540), id.type="entrez_gene_id")
}
\author{
Mark J Chen \email{chenj99@gene.com}
}
\keyword{gortholog}
\keyword{ortholog}


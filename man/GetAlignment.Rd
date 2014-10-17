\name{GetAlignment}
\alias{GetAlignment}
\title{Get alignment}
\usage{
GetAlignment(query.id = query.id, id.type = c("entrez_gene_id",
  "ensembl_pep_id", "ensembl_gene_id"), taxid = taxid, clusterid.cidx = 1,
  taxid.cidx = 2, id.cidx = 3, no.id = NA, diagnose = FALSE,
  url.base = "http://resdev.gene.com/gOrtholog/")
}
\arguments{
  \item{query.id}{query IDs}

  \item{id.type}{type.id which type of IDs are used: Entrez
  Gene ID, Ensembl gene ID, Ensembl peptide ID}

  \item{taxid}{query species}
}
\value{
an object of class AAMultipleAlignment (see
\code{\link{MultipleAlignment-class}})
}
\description{
Given a query, fetch the alignment of its orthologous
cluster from gOrtholog website
}
\examples{
ptpn4a1.aln <- GetAlignment(query.id="ENSP00000359685", id.type="ensembl_pep_id")
}
\author{
Mark J Chen \email{chenj99@gene.com}
}
\seealso{
\code{\link{MultipleAlignment-class}}
}
\keyword{gortholog}
\keyword{ortholog}


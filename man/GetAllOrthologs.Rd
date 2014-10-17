\name{GetAllOrthologs}
\alias{GetAllOrthologs}
\title{Get all orthologs in given species}
\usage{
GetAllOrthologs(taxid = taxid, id.type = c("entrez_gene_id",
  "ensembl_pep_id", "ensembl_gene_id"), taxid.ridx = 2)
}
\arguments{
  \item{taxid}{NCBI taxonomy ID of given species}

  \item{type.id}{which type of IDs are used: Entrez Gene
  ID, Ensembl gene ID, Ensembl peptide ID}
}
\value{
data frame of three columns: cluster ID, taxonomy ID,
Ensembl peptide ID.
}
\description{
Get all orthologs in given species by taxonomy IDs.
}
\examples{
human.mouse <- GetAllOrthologs(taxid=c(9606, 10090), id.type="ensembl_pep_id")
scer.pombe <- GetAllOrthologs(taxid=c(4932, 4896), id.type="entrez_gene_id")
}
\author{
Mark J Chen \email{chenj99@gene.com}
}
\keyword{gortholog}
\keyword{ortholog}


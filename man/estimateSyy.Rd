\name{estimateSyy}
\alias{estimateSyy}

\title{Function to estimate correlations between phenotypic variables from summary statistics}

\description{
This function computes phenotypic correlation
matrix \code{S_YY} based on univariate summary statistics \code{S_XY}.
}

\usage{
estimateSyy( S_XY )
}



\arguments{
\item{S_XY}{
  Univariate summary statistics.

  Data frame with row names corresponding to SNP IDs (e.g., position or rs_id)
  and the following columns:

  - \code{allele_0} - string composed of "A", "C", "G" or "T",

  - \code{allele_1} - string composed of "A", "C", "G" or "T",

  - then, two columns for each trait (phenotypic variable) to be included in the
  analysis; in turn:

  1) \code{traitID_b} with linear regression coefficients,

  2) \code{traitID_se} with corresponding standard errors

  ("traitID" in the column name must be an ID of a trait specified by a user;
  do not use underscores "_" in trait IDs outside "_b"/"_se" in order for the
  IDs to be processed correctly).}
}


\value{
\item{S_YY}{ Matrix containing correlations between traits given as input.
             Row and column names correspond to trait IDs.}


}
\references{
Cichonska et al. (2016) metaCCA: Summary statistics-based multivariate
meta-analysis of genome-wide association studies using canonical
correlation analysis. Bioinformatics, btw052 (in press, to be updated).
}

\author{
Anna Cichonska
}

\note{
In practice, summary statistics of at least one chromosome should be used
in order to ensure good quality of the estimate of phenotypic correlation
structure.
}



\examples{
# Estimating correlations between 10 traits given their
# univariate summary statistics across 1000 SNPs
S_YY = estimateSyy( S_XY = S_XY_full_study1 )

# Viewing the resulting phenotypic correlation matrix
print( S_YY, digit = 3 )

}



\keyword{GenomeWideAssociation}

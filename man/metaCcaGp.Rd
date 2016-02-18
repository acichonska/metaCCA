\name{metaCcaGp}
\alias{metaCcaGp}

\title{
Function to perform genotype-phenotype association analysis according
to metaCCA algorithm.
}

\description{
This function performs genotype-phenotype association analysis according
to metaCCA algorithm (univariate summary statistics-based analysis
of a single or multiple genome-wide association studies (GWAS)
that allows multivariate representation of both genotype and phenotype).

The function accepts a varying number of arguments, depending on the type
of the analysis. By default, single-SNP--multi-trait association analysis
is performed, where each given SNP is tested against all given phenotypic
variables. Other options are to perform single-SNP--multi-trait analysis
of one selected SNP, as well as multi-SNP--multi-trait analysis.
}


\usage{
metaCcaGp( nr_studies, S_XY, std_info, S_YY, N, analysis_type, SNP_id, S_XX )
}



\arguments{

    \item{nr_studies }{
        Number of studies to be analysed.
    }


    \item{S_XY}{

        Univariate summary statistics of the variables to be analysed.

        A list of data frames (one for each study) with row names corresponding
        to SNP IDs (e.g., position or rs_id) and the following columns:

        - \code{allele_0} - string composed of "A", "C", "G" or "T",

        - \code{allele_1} - string composed of "A", "C", "G" or "T",

        - then, two columns for each trait (phenotypic variable) to be included in the
        analysis; in turn:

        1) \code{traitID_b} with linear regression coefficients,

        2) \code{traitID_se} with corresponding standard errors

        ("traitID" in the column name must be an ID of a trait specified by a user;
        do not use underscores "_" in trait IDs outside "_b"/"_se" in order for the
        IDs to be processed correctly).
        }



    \item{std_info}{

        A vector with numerical values \code{0/1} (one value for each study)
        indicating if the univariate analysis has been performed
        on standardised (\code{1}) or non-standardised (\code{0}) data;

        (most likely the data were not standardised - the genotypes were not
         standardised before univariate regression coefficients and standard errors
         were computed - option \code{0} should be used).
    }



    \item{S_YY}{

        A list of phenotypic correlation matrices (one for each study)
        estimated using \code{estimateSyy} function.
    }



    \item{N}{

        A vector with numbers of individuals in each study.
    }


    \strong{\emph{Arguments below are OPTIONAL and depend on the type of the analysis.}}

    \item{analysis_type}{

        Indicator of the analysis type.

        1) Single-SNP--multi-trait analysis of one selected SNP: \code{1}.

        2) Multi-SNP--multi-trait analysis: \code{2}.
    }


    \item{SNP_id}{

    1) Single-SNP--multi-trait analysis of one selected SNP:

    An ID of the SNP of interest.

    2) Multi-SNP--multi-trait analysis:

    A vector with IDs of SNPs to be analysed jointly.
    }



    \item{S_XX}{

        A list of data frames (one for each study) containing correlations
        between SNPs. Row names (and, optionally, column names) must correspond
        to SNP IDs. This argument needs to be given only in case
        of multi-SNP--multi-trait analysis.
    }

}




\value{
    \item{result}{
    Data frame with row names corresponding to SNP IDs.

    Two columns contain:

    1) \code{r_1} - leading canonical correlation value,

    2) \code{-log10(p-val)} - p-value in the -log10 scale.
    }
}



\references{
Cichonska et al. (2016) metaCCA: Summary statistics-based multivariate
meta-analysis of genome-wide association studies using canonical
correlation analysis. Bioinformatics, btw052 (in press, to be updated).
}
\author{
Anna Cichonska
}


\examples{

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#           Analysis of one study according to metaCCA algorithm.           #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# Default single-SNP--multi-trait analysis.
# Here, we will test each of 10 SNPs for an association with a set of 10 traits.
result1 = metaCcaGp( nr_studies = 1,
                     S_XY = list( S_XY_study1 ),
                     std_info = 0,
                     S_YY = list( estimateSyy(S_XY_full_study1 ) ),
                     N = N1 )

# Viewing association results
print( result1, digits = 3 )


# Single-SNP--multi-trait analysis of one selected SNP.
# Here, we will test one of 10 SNPs for an association with a set of 10 traits.
result2 = metaCcaGp( nr_studies = 1,
                     S_XY = list( S_XY_study1 ),
                     std_info = 0,
                     S_YY = list( estimateSyy(S_XY_full_study1) ),
                     N = N1,
                     analysis_type = 1,
                     SNP_id = 'rs80' )

# Viewing association results
print( result2, digits = 3 )


# Multi-SNP--multi-trait analysis.
# Here, we will test a set of 5 SNPs for an association with a set of 10 traits.
result3 = metaCcaGp( nr_studies = 1,
                     S_XY = list( S_XY_study1 ),
                     std_info = 0,
                     S_YY = list( estimateSyy(S_XY_full_study1) ),
                     N = N1,
                     analysis_type = 2,
                     SNP_id = c( 'rs10', 'rs80', 'rs140', 'rs170', 'rs172' ),
                     S_XX = list( S_XX_study1 ) )

# Viewing association results
print( result3, digits = 3 )




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#       Meta-analysis of two studies according to metaCCA algorithm.        #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# Default single-SNP--multi-trait analysis.
# Here, we will test each of 10 SNPs for an association with a set of 10 traits.
meta_result1 = metaCcaGp( nr_studies = 2,
                          S_XY = list( S_XY_study1, S_XY_study2 ),
                          std_info = c( 0, 0 ),
                          S_YY = list( estimateSyy(S_XY_full_study1),
                                       estimateSyy(S_XY_full_study2) ),
                          N = c( N1, N2 ) )

# Viewing association results
print( meta_result1, digits = 3 )


# Single-SNP--multi-trait analysis of one selected SNP.
# Here, we will test one of 10 SNPs for an association with a set of 10 traits.
meta_result2 = metaCcaGp( nr_studies = 2,
                          S_XY = list( S_XY_study1, S_XY_study2 ),
                          std_info = c( 0, 0 ),
                          S_YY = list( estimateSyy(S_XY_full_study1),
                                       estimateSyy(S_XY_full_study2) ),
                          N = c( N1, N2 ),
                          analysis_type = 1,
                          SNP_id = 'rs80' )

# Viewing association results
print( meta_result2, digits = 3 )


# Multi-SNP--multi-trait analysis.
# Here, we will test a set of 5 SNPs for an association with a set of 10 traits.
meta_result3 = metaCcaGp( nr_studies = 2,
                          S_XY = list( S_XY_study1, S_XY_study2 ),
                          std_info = c( 0, 0 ),
                          S_YY = list( estimateSyy(S_XY_full_study1),
                                       estimateSyy(S_XY_full_study2) ),
                          N = c( N1, N2 ),
                          analysis_type = 2,
                          SNP_id = c( 'rs10', 'rs80', 'rs140', 'rs170', 'rs172' ),
                          S_XX = list( S_XX_study1, S_XX_study2 ) )

# Viewing association results
print( meta_result3, digits = 3 )

}



\keyword{GenomeWideAssociation}
\keyword{SNP}
\keyword{Genetics}

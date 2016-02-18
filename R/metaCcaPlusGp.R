# This function performs genotype-phenotype association analysis according
# to metaCCA+ algorithm.


metaCcaPlusGp <- function(nr_studies, S_XY, std_info, S_YY, N, analysis_type, SNP_id, S_XX){

    if ( identical( class(nr_studies), "numeric" ) == FALSE ){
        stop('The first argument (number of studies analysed) needs to be
             numerical. ')
    }

    nr_studies = c( nr_studies, "plus_mode" )

    result  =  metaCcaGp( nr_studies, S_XY, std_info, S_YY, N, analysis_type, SNP_id, S_XX )


    return( result )
}

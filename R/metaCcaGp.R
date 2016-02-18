# This function performs genotype-phenotype association analysis according
# to metaCCA algorithm.


metaCcaGp <- function(nr_studies, S_XY, std_info, S_YY, N, analysis_type, SNP_id, S_XX){

    # define if the "plus" mode was called
    if (  length(nr_studies)==2 && nr_studies[2] == "plus_mode" ){
        # metaCCA+
        plus_mode 	=  1
        nr_studies  =  as.numeric(nr_studies[1])
    } else
        plus_mode 	=  0


    if ( identical( class(nr_studies), "numeric" ) == FALSE ){
        stop('The first argument (number of studies analysed) needs to be
             numerical. ')
    }


    # Identifying analysis type
    if ( missing(analysis_type) ){
        option = 0
    } else {
        option = analysis_type
        if (option != 1 &&  option != 2){
            stop('Wrong indicator of the analysis type. Please provide >>1<<
                 for a single-SNP analysis of the SNP of interest,
                 >>2<< for a multi-SNP analysis.')
        }
    }



    # Validate if the number of elements in the input lists correspond to the
    # number of studies meta-analysed
    if ( length(S_XY) != nr_studies ){
        stop('The number of data frames with univariate summary statistics
             (S_XY) given as input needs to correspond to the number of studies
             analysed (nr_studies).')
    }
    if ( length(std_info) != nr_studies ){
        stop('The length of std_info vector given as input needs to correspond
             to the number of studies analysed (nr_studies).')
    }
    if ( length(S_YY) != nr_studies ){
        stop('The number of phenotypic correlation matrices (S_YY) given
             as input needs to correspond to the number of studies analysed
             (nr_studies).')
    }
    if ( length(N) != nr_studies ){
        stop('The length of N vector with numbers of individuals given as input
             needs to correspond to the number of studies analysed
             (nr_studies).')
    }
    if ( (option == 2)  &&  (length(S_XX) != nr_studies) ){
        stop('The number of data frames containing correlations between SNPs
             (S_XX) given as input needs to correspond to the number of studies
             analysed (nr_studies).')
    }



    trait_id_syy	 		= vector("list", nr_studies)
    trait_ids_betas			= vector("list", nr_studies)
    trait_ids_se			= vector("list", nr_studies)
    SNPid	                = vector("list", nr_studies)

    for (i in 1:nr_studies){

        # Phenotypic correlation matrix - trait IDs
        trait_id_syy[[i]]	=  rownames( S_YY[[i]] )

        # Summary statistics
        # Betas - trait IDs
        trait_ids_betas_temp  =  colnames(S_XY[[i]])[
                                    seq(3,length(S_XY[[i]]),by=2) ]
        trait_ids_betas[[i]]  =  sapply( strsplit(trait_ids_betas_temp, "_") ,
                                    "[[", 1)
        # SE - trait IDs
        trait_ids_se_temp  	  =  colnames(S_XY[[i]])[
                                    seq(4,length(S_XY[[i]]),by=2) ]
        trait_ids_se[[i]]     =  sapply( strsplit(trait_ids_se_temp, "_") ,
                                    "[[", 1)

        # SNP IDs
        SNPid[[i]]	 	=  rownames(S_XY[[i]])


        # Validating if trait IDs of the corresponding regression coefficients
        # and standard errors match
        if ( identical(trait_ids_betas[[i]], trait_ids_se[[i]]) == FALSE ){
            stop('Trait IDs of regression coefficients and standard errors
                 do not match (study ', i, ').')
        }
    }


    # Validating if trait IDs in S_XY match between different studies
    if ( length(unique(trait_ids_betas)) != 1 ){
        stop('Trait IDs in summary statistics of different studies
             do not match.')
    }

    # Validating if trait IDs in S_YY match between different studies
    if ( length(unique(trait_id_syy)) != 1 ){
        stop('Trait IDs in phenotypic correlation structures of different
             studies do not match.')
    }

    # Validating if trait IDs in S_YY and S_XY match
    if ( identical(trait_id_syy[[1]], trait_ids_betas[[1]]) == FALSE ){
        stop('Trait IDs in phenotypic correlation structures and summary
             statistics do not match.')
    }





    # DEFAULT SINGLE-SNP ANALYSIS

    allele0		=	vector("list", nr_studies)
    SXY		    = 	vector("list", nr_studies)
    se			=	vector("list", nr_studies)

    if ( option == 0 ){

        for (i in 1:nr_studies){

            # Validating if alleles are in the correct format
            # allele_0
            out_all0 = valAlleles( levels(S_XY[[i]][,1]) )
            if ( out_all0[[1]] != 0 ) {
                stop("Column 'allele_0' contains at least one invalid
                     character: ", out_all0[[2]] , ", study ", i ,".")
            }
            # allele_1
            out_all1 = valAlleles( levels(S_XY[[i]][,2]) )
            if ( out_all1[[1]] != 0 ) {
                stop("Column 'allele_1' contains at least one invalid
                     character: ", out_all1[[2]] , ", study ", i ,".")
            }

            allele0[[i]] = S_XY[[i]][,1]

            SXY[[i]]  = S_XY[[i]][, seq(3, length(S_XY[[i]]),by=2) ]
            se[[i]]	  = S_XY[[i]][, seq(4, length(S_XY[[i]]),by=2) ]
        }

        #  Validating if SNP IDs match between different studies
        if ( length(unique(SNPid)) != 1 ){
            stop('SNP IDs in summary statistics of different studies
                 do not match.')
        }

        # Validating if allele_0 match between different studies
        if ( length(unique(allele0)) != 1 ){
            stop('Alleles 0 in summary statistics of different studies
                 do not match.')
        }

        # Validating if SNP IDs and trait ids are unique
        if ( length(unique(SNPid[[1]])) != length(SNPid[[1]]) ){
            stop('SNP IDs are not unique! ')
        }
        if ( length(unique(trait_ids_betas[[1]])) != length(trait_ids_betas[[1]]) ){
            stop('Trait IDs are not unique! ')
        }


        # Normalizing regression coefficients (if the univariate analysis has
        # 								been performed on non-standardised data)
        S_XY_norm  =  vector("list", nr_studies)
        for (i in 1:nr_studies){
            if ( std_info[i] == 0 ){
                S_XY_norm[[i]]  =  normalizeSxy(SXY[[i]], se[[i]], N[i]);
            } else if ( std_info[i] == 1 ) {
                S_XY_norm[[i]]  =  SXY[[i]]
            } else
                stop('Wrong indicator of the univariate analysis type
                     (study ', i, '). Please provide >>0<< or >>1<<.')
        }


        # Pooling covariance matrices of the same type
        out_list  	=  poolCov( rep(list(1),nr_studies),  S_YY,  S_XY_norm,  N);
        C_YY 	   	=  out_list[[2]]
        C_XY	  	=  out_list[[3]]
        N_tot     	=  out_list[[4]]


        r_1		=  vector( length = dim(C_XY)[1] )
        p_val	=  vector( length = dim(C_XY)[1] )
        # Analysing one SNP at a time (against all given traits)
        for ( i_snp in 1:dim(C_XY)[1] ){

            t_part    =  matrix( nrow = 1, ncol = dim(C_XY[i_snp,])[2]+1 )
            b_part    =  matrix( nrow = dim(C_XY[i_snp,])[2], ncol = dim(C_YY)[1]+1 )
            full_cov  =  matrix( nrow = dim(t_part)[1]+dim(b_part)[1], ncol = dim(t_part)[2] )
            # Building a full covariance matrix
            t_part  			=  as.matrix( cbind( 1, C_XY[i_snp,] ) )
            colnames(t_part) 	=  NULL; 		rownames(t_part) = NULL
            b_part 				=  as.matrix( cbind( t(C_XY[i_snp,]), C_YY ) )
            colnames(b_part) 	=  NULL; 		rownames(b_part) = NULL
            full_cov			=  rbind( t_part, b_part )

            # Ensuring the PSD property of the full covariance matrix
            if ( plus_mode == 0){
                out_list_m 			=  shrinkPSD(full_cov, 1)
            } else
                out_list_m			=  shrinkPlus(full_cov, 1)
            C_XX_out				=  out_list_m[[1]]
            C_YY_out				=  out_list_m[[2]]
            C_XY_out				=  out_list_m[[3]]

            # Canonical Correlation Analysis (CCA)
            # genotype-phenotype association result
            out_list_cca		=  myCCA(C_XX_out, C_YY_out, C_XY_out, N_tot)
            r_1[i_snp]			=  out_list_cca[[1]]
            p_val[i_snp]		=  -log10(out_list_cca[[2]])

        }


        # Data.frame with the results
        result  			 =  data.frame( r_1, p_val, row.names = SNPid[[1]] )
        colnames(result)[2]  =  '-log10(p-val)'






    # SINGLE-SNP ANALYSIS OF ONE SELECTED SNP
    } else if ( option == 1 ) {

        if ( missing(SNP_id) ){
            stop('Provide an ID of the SNP of interest.')
        }
        if ( length(SNP_id)>1 ){
            stop('An ID of only one SNP should be given.')
        }

        selected_SNPid  =  SNP_id;


        # Summary statistics for a given SNP
        for ( i in 1:nr_studies ){

            # Match SNP ID
            h_id   =  grep( paste("^", selected_SNPid,"$",sep=""), SNPid[[i]] )

            if ( length(h_id) > 1 ) {
                stop('There is more than one SNP with ID "', selected_SNPid,
                     '" in study ',i, '.')
            }
            if ( length(h_id) == 0 ) {
                stop('There is no SNP with ID "', selected_SNPid,
                     '" in study ',i, '.')
            }


            # Validating if alleles are in the correct format
            # allele_0
            out_all0 = valAlleles( levels(S_XY[[i]][h_id, 1]) )
            if ( out_all0[[1]] != 0 ) {
                stop("Column 'allele_0' contains at least one invalid
                     character: ", out_all0[[2]] , ", study ", i ,".")
            }
            # allele_1
            out_all1 = valAlleles( levels(S_XY[[i]][h_id, 2]) )
            if ( out_all1[[1]] != 0 ) {
                stop("Column 'allele_1' contains at least one invalid
                     character: ", out_all1[[2]] , ", study ", i ,".")
            }


            allele0[[i]]   	=  S_XY[[i]][h_id, 1]

            SXY[[i]]		=  S_XY[[i]][h_id, seq(3, length(S_XY[[i]]),by=2) ]
            se[[i]]			=  S_XY[[i]][h_id, seq(4, length(S_XY[[i]]),by=2) ]
        }


        # Validating if allele_0 match between different studies
        if ( length(unique(allele0)) != 1 ){
            stop('Alleles 0 in summary statistics of different studies
                 do not match.')
        }

        # Validating if trait IDs are unique
        if ( length(unique(trait_ids_betas[[1]])) != length(trait_ids_betas[[1]]) ){
            stop('Trait IDs are not unique! ')
        }


        # Normalizing regression coefficients (if the univariate analysis has
        # 							    been performed on non-standardised data)
        S_XY_norm  =  list()
        for (i in 1:nr_studies){
            if ( std_info[i] == 0 ){
                S_XY_norm[[i]]  =  normalizeSxy(SXY[[i]], se[[i]], N[i]);
            } else if ( std_info[i] == 1 ) {
                S_XY_norm[[i]]  =  SXY[[i]]
            } else
                stop('Wrong indicator of the univariate analysis type
                     (study ', i, '). Please provide >>0<< or >>1<<.')
        }


        # Pooling covariance matrices of the same type
        out_list  	=  poolCov( rep(list(1),nr_studies),  S_YY,  S_XY_norm,  N)
        C_YY 	   	=  out_list[[2]]
        C_XY	  	=  out_list[[3]]
        N_tot     	=  out_list[[4]]


        # Building a full covariance matrix
        t_part  			=  as.matrix( cbind( 1, C_XY ) )
        colnames(t_part) 	=  NULL; 		rownames(t_part) = NULL
        b_part 				=  as.matrix( cbind( t(C_XY), C_YY ) )
        colnames(b_part) 	=  NULL; 		rownames(b_part) = NULL
        full_cov			=  rbind( t_part, b_part )


        # Ensuring the PSD property of the full covariance matrix
        if ( plus_mode == 0){
            out_list_m 			=  shrinkPSD(full_cov, 1)
        } else
            out_list_m			=  shrinkPlus(full_cov, 1)
        C_XX_out				=  out_list_m[[1]]
        C_YY_out				=  out_list_m[[2]]
        C_XY_out				=  out_list_m[[3]]


        # Canonical Correlation Analysis (CCA)
        # genotype-phenotype association result
        out_list_cca		=  myCCA(C_XX_out, C_YY_out, C_XY_out, N_tot)
        r_1					=  out_list_cca[[1]]
        p_val				=  -log10(out_list_cca[[2]])


        # Data.frame with the results
        result  		 =  data.frame( r_1, p_val, row.names = selected_SNPid )
        colnames(result)[2]  =  '-log10(p-val)'






    # MULTI-SNP ANALYSIS
    } else if ( option == 2 ) {

        if ( missing(SNP_id) ){
            stop('Provide IDs of SNPs to be analysed jointly.')
        }
        if ( length(SNP_id)<2 ){
            stop('IDs of at least two SNPs should be given in case of multi-SNP
                 analysis.')
        }
        if( missing(S_XX) ){
            stop('Data frames containing correlations between SNPs need
                 to be given.')
        }

        selected_SNPid  = SNP_id

        if ( length(unique(selected_SNPid)) !=  length(selected_SNPid) ){
            stop('IDs of SNPs to be analysed jointly are not unique. ')
        }


        # SNP IDs in S_XX
        S_XX_all		=  vector("list", nr_studies)
        SNPid_inSxx		=  vector("list", nr_studies)
        for ( i in 1:nr_studies ){
            S_XX_all[[i]]		=  as.matrix( S_XX[[i]] )
            SNPid_inSxx[[i]]	=  rownames( S_XX_all[[i]] )
        }


        # Validating if SNP IDs in S_XY match between different studies
        if ( length(unique(SNPid)) != 1 ){
            stop('SNP IDs in summary statistics of different studies
                 do not match.')
        }



        # Summary statistics and genotypic correlation structures corresponding
        # to given SNPs
        SXX  =  vector("list", nr_studies)
        for ( i in 1:nr_studies ){

            # match SNP ids
            h_id_xy = vector( length = length(selected_SNPid) )
            h_id_xx = vector( length = length(selected_SNPid) )
            for ( h in 1:length(selected_SNPid) )	{
                temp_xy  =  grep( paste("^", selected_SNPid[h], "$", sep=""),
                                  SNPid[[i]] )
                temp_xx  =  grep( paste("^", selected_SNPid[h], "$", sep=""),
                                  SNPid_inSxx[[i]] )

                if ( length(temp_xy) == 0 ){
                    stop('SNP "', selected_SNPid[h], '" not found in the summary
                         statistics of study ', i, ' .' )
                } else
                    h_id_xy[h] 	= 	temp_xy

                if ( length(temp_xx) == 0 ){
                    stop('SNP "', selected_SNPid[h], '" not found in the S_XX
                         of study ', i, ' .' )
                } else
                    h_id_xx[h] 	= 	temp_xx
            }


            # Validating if alleles are in the correct format
            # allele_0
            out_all0 = valAlleles( levels(S_XY[[i]][h_id_xy, 1]) )
            if ( out_all0[[1]] != 0 ) {
                stop("Column 'allele_0' contains at least one invalid character: ",
                     out_all0[[2]] , ", study ", i ,".")
            }
            # allele_1
            out_all1 = valAlleles( levels(S_XY[[i]][h_id_xy, 2]) )
            if ( out_all1[[1]] != 0 ) {
                stop("Column 'allele_1' contains at least one invalid character: ",
                     out_all1[[2]] , ", study ", i ,".")
            }


            SXY[[i]]      = S_XY[[i]][h_id_xy, seq(3, length(S_XY[[i]]),by=2) ]
            se[[i]]       = S_XY[[i]][h_id_xy, seq(4, length(S_XY[[i]]),by=2) ]
            allele0[[i]]  = S_XY[[i]][h_id_xy, 1]
            SXX[[i]]      = S_XX_all[[i]][h_id_xx, h_id_xx]
        }


        # Validating if allele_0 match between different studies
        if ( length(unique(allele0)) != 1 ){
            stop('Alleles 0 in summary statistics of different studies
                 do not match.')
        }

        # Validating if trait IDs are unique
        if ( length(unique(trait_ids_betas[[1]])) != length(trait_ids_betas[[1]]) ){
            stop('Trait IDs are not unique! ')
        }


        # Normalizing regression coefficients (if the univariate analysis has
        # 							    been performed on non-standardised data)
        S_XY_norm  =  list()
        for (i in 1:nr_studies){
            if ( std_info[i] == 0 ){
                S_XY_norm[[i]]  =  normalizeSxy(SXY[[i]], se[[i]], N[i]);
            } else if ( std_info[i] == 1 ) {
                S_XY_norm[[i]]  =  SXY[[i]]
            } else
                stop('Wrong indicator of the univariate analysis type
                     (study ', i, '). Please provide >>0<< or >>1<<.')
        }


        # Pooling covariance matrices of the same type
        out_list  	=  poolCov( SXX,  S_YY,  S_XY_norm,  N )
        C_XX		=  out_list[[1]]
        C_YY 	   	=  out_list[[2]]
        C_XY	  	=  out_list[[3]]
        N_tot     	=  out_list[[4]]


        # Building a full covariance matrix
        t_part  			=  as.matrix( cbind( C_XX, C_XY ) )
        colnames(t_part) 	=  NULL; 		rownames(t_part) = NULL
        b_part 				=  as.matrix( cbind( t(C_XY), C_YY ) )
        colnames(b_part) 	=  NULL; 		rownames(b_part) = NULL
        full_cov			=  rbind( t_part, b_part )


        # Ensuring the PSD property of the full covariance matrix
        if ( plus_mode == 0){
            out_list_m 			=  shrinkPSD(full_cov,  length(selected_SNPid))
        } else
            out_list_m			=  shrinkPlus(full_cov, length(selected_SNPid))
        C_XX_out				=  out_list_m[[1]]
        C_YY_out				=  out_list_m[[2]]
        C_XY_out				=  out_list_m[[3]]


        # Canonical Correlation Analysis (CCA)
        # genotype-phenotype association result
        out_list_cca		=  myCCA( C_XX_out, C_YY_out, C_XY_out, N_tot )
        r_1					=  out_list_cca[[1]][1]
        p_val				=  -log10( out_list_cca[[2]][1] )


        # Data.frame with the results
        result 	=  data.frame( r_1, p_val, row.names = list(selected_SNPid) )
        colnames(result)[2]  =  '-log10(p-val)'

    }



    return( result )
}

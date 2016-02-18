# This function uses Pearson correlation to compute phenotypic
# correlation matrix S_YY based on univariate summary statistics S_XY.


estimateSyy <- function(S_XY) {

    # Betas - trait IDs
    trait_ids_betas = colnames(S_XY)[ seq(3, length(colnames(S_XY)), by=2) ]
    trait_ids_betas  = gsub("_b", "", trait_ids_betas)

    # SE - trait IDs
    trait_ids_se = colnames(S_XY)[ seq(4, length(colnames(S_XY)), by=2) ]
    trait_ids_se = gsub("_se", "", trait_ids_se)


    # Validating if trait IDs of the corresponding regression coefficients
    # and standard errors match
    if ( identical(trait_ids_betas, trait_ids_se) == FALSE ) {
        stop("Trait IDs of regression coefficients and standard errors
             do not match")
    }


    # Validating if trait IDs are unique
    if ( length(unique(trait_ids_betas)) != length(trait_ids_betas) ) {
        stop("Trait IDs are not unique")
    }


    # Validating if alleles are in the correct format
    dict = c("A","T","G","C")
    all0 = levels(S_XY[,1])
    all1 = levels(S_XY[,2])
    all0_el = unique(unlist( strsplit(all0, "") ))
    all1_el = unique(unlist( strsplit(all1, "") ))

    if ( length( setdiff( all0_el, dict ) ) != 0 ) {
        stop("Column 'allele_0' contains at least one invalid character: ",
             list(setdiff( all0_el, dict )) , ".")
    }

    if (length( setdiff( all1_el, dict ) ) != 0 ) {
        stop("Column 'allele_1' contains at least one invalid character: ",
             list(setdiff( all1_el, dict )) , ".")
    }



    Betas = S_XY[ seq(3, length(colnames(S_XY)), by=2) ]
    colnames(Betas) = gsub("_b", "", colnames(Betas))

    S_YY = cor(Betas)


    return(S_YY)
}

# This function performs Canonical Correlation Analysis (CCA) based on the
# building blocks of the full covariance matrix (C_XX, C_YY, C_XY),
# as well as significance testing of the canonical correlations.


myCCA <- function(C_XX, C_YY, C_XY, N) {

    C_XX = as.matrix(C_XX);  C_YY = as.matrix(C_YY);

    expM_Cxx = expM(C_XX, -0.5)
    expM_Cyy = expM(C_YY, -0.5)

    K    =  expM_Cxx  %*% C_XY  %*%  expM_Cyy

    # Singular Value Decomposition (SVD)
    svd_out =  svd(K)
    U 		=  svd_out$u
    S		=  as.matrix(svd_out$d)
    V		=  svd_out$v

    # S is a vector, create a diagonal matrix
    temp = diag(length(S))
    diag(temp) = S


    min_r 	=  min( dim(C_XX)[1], dim(C_YY)[1] )
    a 		=  array( dim = c( dim(C_XX)[1], min_r) )  	# canonical weights
    b 		=  array( dim = c( dim(C_YY)[1], min_r) )
    r 		=  array( dim = c(1, min_r) )             	# canonical correlation
    wilks 	=  array( dim = c(1, min_r) )         		# Wilk's lambda
    lambdas	=  array( dim = c(length(S), length(S)) )
    lambdas =  temp^2
    df 		=  array( dim = c(1, min_r) )             	# degrees of freedom

    p 		=  dim(C_XX)[1]
    q 		=  dim(C_YY)[1]
    k		=  0;


    for ( i in 1:min_r ) {

        a[,i]	=  expM_Cxx  %*%  U[,i]
        b[,i] 	=  expM_Cyy  %*%  V[,i]
        r[1,i]	=  ( t(a[,i])%*%C_XY%*%b[,i] ) /
                    ( sqrt( t(a[,i])%*%C_XX%*%a[,i] )  *
                    sqrt( t(b[,i])%*%C_YY%*%b[,i] ) )

        # Wilk's lambda
        prodd = 1
        for (j in i:min_r) {
            prod_temp 	=  1 - lambdas[j,j]
            prodd 	=  prodd * prod_temp
        }

        wilks[i] 	=  prodd
        df[i] 		=  (p-k)*(q-k)       # df for Bartlett Chi-square
        k 			=  k+1
    }


    # Bartlett Chi-square approximation to Wilk's Lambda
    chi 	= 	-( (N-1) - 0.5*(p+q+1) )  *  log(wilks)

    # H0: all canonical correlations = 0
    # p-value
    p_val = 	pchisq(chi, df, lower.tail = FALSE)


    return( list(r, p_val) )
}

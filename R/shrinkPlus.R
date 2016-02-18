# This function shrinks the full covariance matrix beyond the level
# guaranteeing its PSD property.
# It returns building blocks of the resulting shrunken full covariance
# matrix: C_XX, C_XY, C_YY.


shrinkPlus <- function(R, n) {

    # Eigenvalues of the full covariance matrix
    R		=  as.matrix(R)
    lambdas	=  eigen(R)$values

    # Building blocks of the full covariance matrix
    C_XX		=  R[1:n,  1:n]
    C_YY		=  R[(n+1):dim(R)[1],  (n+1):dim(R)[1]]
    C_XY		=  R[1:n,  (n+1):dim(R)[1]]


    d 		=  vector()
    d[1]	=  10^10
    K		=  (expM(C_XX, -0.5))  %*% C_XY  %*%  (expM(C_YY, -0.5))
    d[2]	=  max( svd(K)$d )
    i		=  1


    # Finding an appropriate shrinkage level
    opt_threshold  =  optimalThreshold( d, R, n )


    # Shrinkage
    while(  ( (abs(d[i]-d[i+1]))/d[i]  > opt_threshold )  ||
            ( min(lambdas)<0 ) ){

        i 			=  i+1
        R 			=  0.999*R
        diag(R) 	=  1
        lambdas 	=  eigen(R)$values

        C_XX_new_gen	=  R[1:n,  1:n]
        C_YY_new_gen	=  R[(n+1):dim(R)[1],  (n+1):dim(R)[1]]
        C_XY_new_gen	=  R[1:n,  (n+1):dim(R)[1]]

        K			=  (expM(C_XX_new_gen, -0.5))  %*% C_XY_new_gen  %*%
                        (expM(C_YY_new_gen, -0.5))

        d[i+1]		=  max( svd(K)$d )
    }


    C_XX_shr	=	R[1:n,  1:n]
    C_YY_shr	=	R[(n+1):dim(R)[1],  (n+1):dim(R)[1]]
    C_XY_shr	=	R[1:n,  (n+1):dim(R)[1]]


    return( list(C_XX_shr, C_YY_shr, C_XY_shr) )
}

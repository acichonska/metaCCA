# This function finds an appropriate shrinkage level of the full
# covariance matrix using an elbow criterion.


optimalThreshold <- function( d_full, R, n ) {

    # Monitoring the behaviour of the leading canonical correlation
    # value during shrinkage of the full covariance matrix
    i 	=  1
    while ( ( abs( d_full[length(d_full)-1] - d_full[length(d_full)] )
              / d_full[length(d_full)-1] )  >  0.005 ) {

        i 			=  i+1
        R			=  0.999*R
        diag(R) 	=  1

        C_XX_new_temp	=  R[1:n,  1:n]
        C_YY_new_temp	=  R[(n+1):dim(R)[1],  (n+1):dim(R)[1]]
        C_XY_new_temp	=  R[1:n,  (n+1):dim(R)[1]]

        K			=  (expM(C_XX_new_temp, -0.5))  %*% C_XY_new_temp  %*%
                        (expM(C_YY_new_temp, -0.5))
        d_full[i+1]	=  max( svd(K)$d )
    }



    # Finding an appropriate amount of shrinkage
    d_temp = abs(diff(d_full));  d_temp = d_temp[3:length(d_temp)]
    d_temp = d_temp/d_full[ 3:(length(d_full)-1) ]

    if ( length(d_temp) < 2 ){
        opt_threshold = 0.005
    } else

    x  	=  c( 2:(length(d_temp)+1) )
    y  	=  d_temp

    # Normalizing x and y to [0,1]
    x_n 	=  (x-min(x))/(max(x)-min(x))
    y_n 	=  (y-min(y))/(max(y)-min(y))

    # Computing Euclidean distance of each point to the origin of the plot
    # of the percent change of canonical correlations between subsequent
    # shrinkage iterations versus iteration number
    dist_to_origin  =  vector()
    for ( ic in 1:length(x) ){
        temp				=  c( x_n[ic], y_n[ic] )
        temp_diff			=  temp - c(0,0)
        # 2-norm
        dist_to_origin[ic]	=  sum(abs(temp_diff)^2)^(1/2)
    }

    # Taking the point which has the smallest distance to the origin
    # of the plot
    id_min_dist 	=  which.min(dist_to_origin)

    opt_threshold	=  y[id_min_dist]


    return( opt_threshold )
}

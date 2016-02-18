# This function pools covariance matrices of the same type.
# There are 4 types of inputs: S_XX, S_YY, S_XY, N.
# Each is given in the form of list, and each element of the list corresponds
# to a data.frame.


poolCov <- function(S_XX, S_YY, S_XY, N) {

    nr_s  = length(N);       # number of studies

    counter_XX    =  S_XX[[1]]
    counter_XX[]  =  0
    counter_YY    =  S_YY[[1]]
    counter_YY[]  =  0
    counter_XY    =  S_XY[[1]]
    counter_XY[]  =  0
    denom         =  0
    sum_N         =  0


    for ( i in 1:nr_s ){
        counter_XX  =  counter_XX  +  ( (N[i]-1) * S_XX[[i]] )
        counter_YY  =  counter_YY  +  ( (N[i]-1) * S_YY[[i]] )
        counter_XY  =  counter_XY  +  ( (N[i]-1) * S_XY[[i]] )

        denom  =  denom + (N[i]-1)

        sum_N  =  sum_N + N[i]
    }


    C_XX  =  counter_XX/denom;
    C_YY  =  counter_YY/denom;
    C_XY  =  counter_XY/denom;


    return( list(C_XX, C_YY, C_XY, sum_N) )
}

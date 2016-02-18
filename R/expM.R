# Matrix power


expM <- function(X, e) {

    v	=  La.svd(X)
    res	=  v$u %*% diag(v$d^e) %*% v$vt

    return(res)
}

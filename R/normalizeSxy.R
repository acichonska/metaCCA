# This function normalizes regression coefficients
# using their standard errors 'se' and the number of samples 'N'.


normalizeSxy <- function(S_XY_raw, se, N) {

    S_XY =  (1/sqrt(N)) * (S_XY_raw/se)

    return(S_XY)
}

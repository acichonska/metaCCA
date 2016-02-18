# This function validates if alleles are in the correct format.


valAlleles <- function(al) {

    # allowed characters
    dict = c("A","T","G","C")

    al_el = unique(unlist( strsplit(al, "") ))

    o1 = length( setdiff( al_el, dict ) )
    o2 = list(setdiff( al_el, dict ))
    return( list(o1,o2) )

}

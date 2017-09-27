packIndices <- function(assignments) 
# This function compresses the assignment vectors.
#
# written by Aaron Lun
# created 28 November 2016
{
    out <- .Call(cxx_pack_indices, assignments, TRUE)
    return(out)
}

unpackIndices <- function(assignments) 
# This function decompresses the assignment vectors.
#
# written by Aaron Lun
# created 28 November 2016    
{
    out <- .Call(cxx_pack_indices, assignments, FALSE)
    return(out)
}

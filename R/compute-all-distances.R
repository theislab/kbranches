#'Compute sample distances
#'
#'Computes all SQUARED pairwise distances of data points.
#'
#'@param input_dat: data frame of input data with rows=samles and cols=dimensions
#'
#'@return Dist: a matrix containing all sample distances (square Euclidean)
#'
#'@examples
#'\dontrun{
#'  Dist=compute_all_distances(scdata.3lines)
#'}
#'@export
compute_all_distances=function(input_dat)
{
  #   compute all SQUARED pairwise distances of data points
  #   and return them as matrix 'Dist'
  #
  #   inputs:
  #     -input_dat: data frame of input data with rows=samles and cols=dimensions
  #   return values:
  #     -Dist: matrix containing sample distances
  return(as.matrix(stats::dist(input_dat)^2))#squared Euclidean distance matrix
}

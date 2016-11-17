#' Scale the neighbourhood size
#'
#'scales the number of neighbors (neighbourhood size) depending on the local sample density.
#'keeps increasing K until the distance to the furthest of the K neighbors is above the dmin threshold
#'
#'@param input_dat: data frame of input data with rows=samles and cols=dimensions
#'@param n: position (row) of the sample, in the input_dat data frame
#'@param K_init: initial number of nearest neighbors that needs to be scaled
#'@param Dist: distance matrix of all samples
#'@param dmin: distance threshold (to the furthest neighbor) to be reached
#'
#'@return Kscaled: the scaled neighbourhood size
#'
#'@examples
#'\dontrun{
#' Dist=compute_all_distances(scdata.3lines)
#'
#' #find the 5 nearest neghbours of the second sample
#' neibs=find_K_nearest(n=2,K=5,Dist=Dist)
#'}
#'@export
scale_S_neib=function(input_dat,n,K_init,Dist,dmin)
{
  #   scales the number of neighbors (neighbourhood size) depending on the local sample density
  #   keeps increasing K until the distance to the furthest of the K neighbors is above the dmin threshold
  #
  #   inputs:
  #     -input_dat: data frame of input data with rows=samles and cols=dimensions
  #     -n: position (row) of the sample, in the input_dat data frame
  #     -K_init: initial number of nearest neighbors that needs to be scaled
  #     -Dist: distance matrix of all samples
  #     -dmin: distance threshold (to the furthest neighbor) to be reached
  #
  #   return values:
  #     -Kscaled: numeric, upscaled number of nearest neighbors

  N=nrow(input_dat)#number of samples
  all_neibs=find_K_nearest(n=n,K=N-1,Dist=Dist)#K=N-1, because we exclude the sample itself from the list of neighbors
  #find effective K until the minimum distance (to the furthest neighbor) is reached
  #Kscaled=S_neib
  Kscaled=K_init
  max_dist=all_neibs$dst[Kscaled]
  while(max_dist<dmin)
  {
    Kscaled=Kscaled+1
    max_dist=all_neibs$dst[Kscaled]
  }
  return(Kscaled)#return the new K
}

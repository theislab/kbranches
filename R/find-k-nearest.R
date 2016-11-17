#'Find K nearest neighbours
#'
#'Finds the K nearest neighbours of a sample, given a matrix containing sample distances.
#'
#'@param n: sample for which the nearest neighbors are computed
#'@param K: number of nearest neighbors to be considered
#'@param Dist: matrix containing sample distances, initially computed by compute_all_distances
#'
#'@return a list with elements:
#'\itemize{
#'  \item - pos: the positions (rows) of the N nearest neighbours.
#'  \item - dst: the distances of the N nearest neighbours.
#'}
#'
#'@examples
#' Dist <- compute_all_distances(scdata.3lines.simulated6genes_subsampled)
#'
#' #find the 5 nearest neghbours of the second sample
#' neibs <- find_K_nearest(n=2, K = 5, Dist = Dist)
#'@export
find_K_nearest=function(n,K=3,Dist)
{
  #   select and return the indexes of the K nearest neighbors of sample n
  #
  #   inputs:
  #     -n: sample for which the nearest neighbors are computed
  #     -K: number of nearest neighbors to be considered
  #     -Dist: matrix containing sample distances, initially computed by compute_all_distances
  #   return values:
  #     -pos: positions of the nearest neighbours
  #     -dst: distances of the nearest neighbours

  ndist=Dist[n,]#distances if all neighbors
  neib_ix=sort.int(ndist,decreasing = FALSE, index.return = TRUE)$ix #indices of all neighbors
  neib_ix=neib_ix[-1]#remove the first element, which is the sample itself with distance 0
  neib_ix=neib_ix[1:K]#keep the K nearest neighbors
  neib_dist=ndist[neib_ix]

  return(list(pos=neib_ix, dst=neib_dist))
}

#calculate the original and modified versions of the GAP statistic, given a clustering
calculate_GAP=function(input_dat,clustering_error,Kappa,cluster_labels,
                       nstart_GAP,B_GAP,fixed_center,medoids,Dmat,init_Kmeans,
                       show_plots_GAP=FALSE,silent=TRUE)
{
  N=nrow(input_dat)
  P=ncol(input_dat)

  #Calculate the GAP
  if(is.null(B_GAP)==FALSE)
  {

    #calculate the class centers
    centers=matrix(NA,nrow = Kappa,ncol = P)
    for(k in 1:Kappa){centers[k,]=colMeans(input_dat[cluster_labels==k,,drop=FALSE])}

    dsp=0#dispersion around the class centers
    for(i in 1:N){ dsp=dsp+norm(input_dat[i,]-(centers[cluster_labels[i],]),'2')^2}#calculate the squared Euclidean distance

    #the observed values
    W_k=sum(dsp)
    logW_k=log(W_k)

    # print('initializing')
    err_b=rep(0,B_GAP)
    err_b_log=rep(0,B_GAP)
    # print(err_b)

    W_k_star=rep(0,B_GAP)
    logW_k_star=rep(0,B_GAP)

    for(b in 1:B_GAP)
    {
      if(silent==FALSE){print(paste('GAP iteration:',b,'of',B_GAP))}
      boot_dat=matrix(data = NA, nrow = N, ncol = P)
      feat_max=apply(input_dat,2,max)#create random samples is the range of min and max features
      feat_min=apply(input_dat,2,min)#in order to get a null distribution
      for(i in 1:N)
      {
        for(j in 1:P)
        {
          boot_dat[i,j]=runif(1,min=feat_min[j],max=feat_max[j])
        }
      }
      boot_dat=as.data.frame(boot_dat)
      # plot(x=boot_dat[,1],y=boot_dat[,2])

      if(medoids==TRUE)#compute the matrix of all distances for the bootstrap data
      {
        D_GAP=compute_all_distances(boot_dat)
      }else
      {
        D_GAP=NULL
      }
      clust2=kbranch.global(Kappa=Kappa, input_dat=boot_dat,c0=NULL,Vmat=NULL,show_plots = show_plots_GAP,silent=silent, nstart = nstart_GAP,fixed_center=fixed_center,medoids=medoids,Dmat=D_GAP)

      err_b[b]=clust2$err
      err_b_log[b]=log(clust2$err)

      #calculate the class centers
      centers=matrix(NA,nrow = Kappa,ncol = P)
      for(k in 1:Kappa){centers[k,]=colMeans(boot_dat[clust2$cluster==k,,drop=FALSE])}

      dsp=0#dispersion around the class centers
      for(i in 1:N){ dsp=dsp+norm(boot_dat[i,]-(centers[clust2$cluster[i],]),'2')^2}#calculate the squared Euclidean distance

      #the expected values
      W_k_star[b]=sum(dsp)
      logW_k_star[b]=log(sum(dsp))

    }

    GAP=sum(err_b)/B_GAP-clustering_error
    GAPl=sum(err_b_log)/B_GAP-log(clustering_error)
    EW_k_star=sum(W_k_star)/B_GAP-W_k
    ElogW_k_star=sum(logW_k_star)/B_GAP-logW_k

    #return arguments
    return(list(      GAP=GAP,
                      GAPl=GAPl,
                      GAP_orig=ElogW_k_star,
                      GAP_orig_no_log=EW_k_star,

                      GAP.sd=sqrt(1+1/B_GAP)*sd(err_b),
                      GAPl.sd=sqrt(1+1/B_GAP)*sd(err_b_log),
                      GAP_orig.sd=sqrt(1+1/B_GAP)*sd(logW_k_star),
                      GAP_orig_no_log.sd=sqrt(1+1/B_GAP)*sd(W_k_star)
    )
    )
  }else#do not compute the GAP statistic
  {
    return(list(      GAP=NULL,
                      GAPl=NULL,
                      GAP_orig=NULL,
                      GAP_orig_no_log=NULL,

                      GAP.sd=NULL,
                      GAPl.sd=NULL,
                      GAP_orig.sd=NULL,
                      GAP_orig_no_log.sd=NULL
    )
    )
  }
}

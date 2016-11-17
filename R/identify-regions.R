#'Identify regions of interest
#'
#'Identifies regions of interest (the tips of branches, or branching regions)
#'based on comparing the GAP scores#'acquired by running local clustering by 'kbranch.local'.
#'Performs local filtering to reduce noise in the extracted labels.
#'
#'@param input_dat: data frame of input data with rows=samles and cols=dimensions.
#'@param mode: =c('tip','branch') find either the tips or the branching regions
#'@param tip_mode: =c('2','3','both') find the tips using the GAP statistic of 1 vs 2,3 or both
#'@param gap_scores: list, output of function 'kbranch.local' for the same 'input_dat'
#'@param smoothing_region: number of neighbours in to consider for label filtering
#'@param smoothing_region_thresh: minimum number of thresholds with same label in the neighbourhood required for the sample to keep it's label
#'@param Dist: matrix of sample to sample distances
#'@param dotsize1: size of points of class1 (red), used when plotting in 3D
#'@param dotsize2: size of points of class2 (green), used when plotting in 3D
#'@param nclust: number of clusters (tips/branching regions). If left NULL, will be estimated.
#'@param nclust.max: maximum possible number of clusters to consider. Only used if nclust==NULL
#'@param B.max: maximum number of bootstrap datasets used to calculate the GAP statistic to estimate nclust. Only used if nclust==NULL
#'@param SE.factor: used to estimate nclust, argument to cluster::maxSE. Only used if nclust==NULL
#'@param repeats: number of times to estimate nclust (for stability - the nclust most frequently considered 'best' is finally extracted). Only used if nclust==NULL
#'
#'@return a list with elements:
#'\itemize{
#'  \item - is_in_region: a logical vector indicating which samples are part of the region (TRUE) and which not (FALSE)
#'  \item - is_in_region_filtered: same as 'is_in_region', but after S_neib filtering to reduce noise.
#'  \item - cluster: cluster assignment for each of the data points in is_in_region_filtered
#'  \item - nclust: number of clusters identified, or specified by the user
#'  \item - cluster_frequency: how many times (out of 'repeats') each number of clusters (of the nclust.max) was considered 'best'.
#'  NULL if nclust was provided by the user
#'}
#'
#'@examples
#'
#'see example of kbranch.local
#'
#'@export
identify_regions=function(input_dat,mode='tip',tip_mode='3',gap_scores=NULL,smoothing_region=NULL,smoothing_region_thresh=NULL,Dist,dotsize1=7,dotsize2=7,
                          nclust=NULL,nclust.max=5, B.max=100, SE.factor = 2.58, repeats = 20)#scale the number of neighbors depending on the local density)
{

  #   finds and filters the tips of branches, or branching regions
  #
  #   inputs:
  #     -input_dat: data frame of input data with rows=samles and cols=dimensions
  #     -mode: =c('tip','branch') find either the tips or the branching regions
  #     -tip_mode: =c('2','3','both') find the tips using the GAP statistic of 1 vs 2,3 or both
  #     -gap_scores: output of function 'identify_region' for the same 'input_dat'
  #     -smoothing_region: number of neighbours in to consider for label filtering
  #     -smoothing_region_thresh: minimum number of thresholds with same label in the neighbourhood required for the sample to keep it's label
  #     -Dist: matrix of sample to sample distances
  #     -dotsize1: size of points of class1, used for plotting
  #     -dotsize2: size of points of class2, used for plotting
  #
  #   return values:
  #     -is_in_region: a logical vector indicating which samples are part of the region (TRUE) and which not (FALSE)
  #     -is_in_region_filtered: same as 'is_in_region', but after S_neib filtering to reduce noise.

  N=nrow(input_dat)#number of samples
  P=ncol(input_dat)#number of dimensions

  xlim=c(min(input_dat[,1]),max(input_dat[,1]))
  ylim=c(min(input_dat[,2]),max(input_dat[,2]))
  if(P==3){zlim=c(min(input_dat[,3]),max(input_dat[,3]))}else(zlim=NULL)#if there is a 3rd dimension


  if(mode=='tip')
  {
    gscore=gap_scores$gscore_orig #original gap score
    is_edge_1v2=gscore[,1]>gscore[,2]
    is_edge_1v3=gscore[,1]>gscore[,3]
    if(tip_mode=='2')
    {
      is_edge=is_edge_1v2
    }else if(tip_mode=='3')
    {
      is_edge=is_edge_1v3
    }else
    {
      is_edge=is_edge_1v2|is_edge_1v3#logical OR: is edge if it is edge for 1v2 OR for 1v3
    }
  }else#branching region
  {
    gscore=gap_scores$gscore #modified gap score
    is_edge=gscore[,3]>gscore[,2]
  }

  tmp=is_edge

  filter_internal=function(smoothing_region,smoothing_region_thresh,keep_fixed=TRUE,show_plots=TRUE)
  {
    #filter the color of each samples, based on the majority of it's smoothing_region nearest neighbors
    if(smoothing_region>0)
    {
      is_edge_filtered=rep(F,N)
      for(i in 1:N)
      {
        nnres=find_K_nearest(n=i,K=smoothing_region,Dist=Dist)
        nn=sort(nnres$pos)
        tmp_nn=tmp[nn]

        #count the number of smoothing_region nearest neighbors that have the same label as the sample in question
        #same_neighbors=sum(length(which(tmp_nn==tmp[i])))
        same_neighbors=sum(length(which(tmp_nn==TRUE)))#how many are 'active' tip/branch
        if(same_neighbors>=smoothing_region_thresh)#if above the threshold, the sample keeps it's label
        {
          #is_edge_filtered[i]=tmp[i]
          is_edge_filtered[i]=TRUE
        }else
        {
          #if below the threshold, the sample gets the label of the majority
          rsum=sum(length(which(tmp_nn==T)));#print(gsum)
          gsum=sum(length(which(tmp_nn==F)));#print(rsum)
          if(rsum>gsum)
          {
            is_edge_filtered[i]=T #red
          }else
          {
            is_edge_filtered[i]=F #green
          }
        }
      }
    }else#do not perform filtering, show unfiltered
    {
      is_edge_filtered=is_edge
    }


    if(show_plots==TRUE)
    {
      # main=paste('Branch Tips in red filtered using',smoothing_region_thresh,'of',smoothing_region,'NN')
      main=''
      red_paint=colorRampPalette(RColorBrewer::brewer.pal(3,'Reds'))(3)[3]#colorblind_friendly red
      green_paint=colorRampPalette(RColorBrewer::brewer.pal(3,'Greens'))(3)[3]#colorblind_friendly green
      clr_filtered=ifelse(is_edge_filtered==T,red_paint,green_paint)
      #clr_filtered=ifelse(is_edge_filtered==T,'red','green')
      if(P==2)
      {
        plot(input_dat,main=main)
        points(input_dat,col=clr_filtered)
      }else if(P==3)
      {

        rgl.clear()
        plot3d(x=input_dat[tmp==T,1],y=input_dat[tmp==T,2],z=input_dat[tmp==T,3],xlab='x axis',ylab='y axis',
               zlab='z axis',col=clr_filtered[tmp==T],size=dotsize1,
               main=main,
               xlim=xlim,ylim=ylim,zlim=zlim)
        points3d(x=input_dat[tmp==F,1],y=input_dat[tmp==F,2],z=input_dat[tmp==F,3],size=dotsize2,col=clr_filtered[tmp==F])

      }
    }
    return(list(is_in_region=is_edge,is_in_region_filtered=is_edge_filtered))
    # return(is_edge_filtered)
    # return(list(smoothing_region=smoothing_region,smoothing_region_thresh=smoothing_region_thresh,keep_fixed=keep_fixed,is_edge=is_edge,is_edge_filtered=is_edge_filtered))
  }

  fCallback = function( lastTouched )
  {
    if((lastTouched=="smoothing_region")&(guiGetValue("keep_fixed")==TRUE))
    {
      guiSetValue("smoothing_region_thresh",guiGetValue("smoothing_region"))
    }
  }

  if((is.null(smoothing_region)==FALSE)&(is.null(smoothing_region_thresh)==FALSE))#smoothing_region and smoothing_region_thresh have been provided
  {
    #if(P==3){open3d()}#open a new plot window
    #is_edge_filtered=filter_internal(smoothing_region,smoothing_region_thresh)
    #retval=list(is_edge=is_edge,is_edge_filtered=is_edge_filtered)
    retval=filter_internal(smoothing_region,smoothing_region_thresh,show_plots=FALSE)
  }else
  {
    if(P==2){x11()}else if(P==3){open3d()}#open a new plot window
    filter_res = gui(filter_internal,argSlider = list(smoothing_region=c(0,N,1),smoothing_region_thresh=c(1,N,1)),argOption=list(keep_fixed=c('TRUE','FALSE')),
                     title='Select S_neib',output = NULL,closeOnExec = FALSE, callback = fCallback, argType = list(show_plots='i'))
    #print(y)
    #retval=filter_region(input_dat,mode,tip_mode,gap_scores,smoothing_region=filter_res$smoothing_region,smoothing_region_thresh=filter_res$smoothing_region_thresh,Dist,dotsize1,dotsize2)
    retval=filter_internal(smoothing_region=filter_res$smoothing_region,smoothing_region_thresh=filter_res$smoothing_region_thresh,show_plots = FALSE)
  }
  #print(filter_res)

  if(is.null(nclust))#estimate the number of clusters, if not given
  {
    #now assign samples to specific branches/tips
    K_hist=rep(0,nclust.max)
    for(i in 1:repeats)
    {
      #calculate GAP statistic
      gap=cluster::clusGap(input_dat[retval$is_in_region_filtered,], FUNcluster = kmeans, K.max = nclust.max, B = B.max, d.power = 2,verbose = F)
      #estimate number of clusters
      if(mode=='tip')
      {
        nclust=cluster::maxSE(gap$Tab[,3],gap$Tab[,4],
                              method = "globalSEmax",SE.factor = 2.58)#2.58: 99% conf. int.
      }else
      {
        nclust=cluster::maxSE(gap$Tab[,3],gap$Tab[,4],
                              method = "firstSEmax",SE.factor = 2.58)#2.58: 99% conf. int.
      }
      K_hist[nclust]=K_hist[nclust]+1
    }
    nclust=which.max(K_hist)#which number of clusters is observed most of the time
  }else#the number of clusters was specified by the user
  {
    K_hist=NULL
  }


  # if(mode=='tip')
  # {
  #   gap=cluster::clusGap(input_dat[retval$is_in_region_filtered,],FUNcluster = kmeans, nclust.max = nclust.max, B=B.max,d.power = 2)
  #   #assume that there are >= 2 tips, so remove gap row 1 which corresponds to 1 tip only
  #   nclust=cluster::maxSE(gap$Tab[-1,3],gap$Tab[-1,4],method = "firstSEmax",SE.factor = SE.factor)#suggested number of tips, 2.58: 99% conf. int.
  #   nclust=nclust+1#correct the number of clusters
  # }else#branching region
  # {
  #   gap=cluster::clusGap(input_dat[retval$is_in_region_filtered,],FUNcluster = kmeans, nclust.max = nclust.max, B=B.max,d.power = 2)
  #   nclust=cluster::maxSE(gap$Tab[,3],gap$Tab[,4],method = "firstSEmax",SE.factor = SE.factor)#suggested number of tips, 2.58: 99% conf. int.
  # }

  #now cluster the samples using the 'optimal' number of 'nclust' clusters (tips/branching regions)
  clusters=stats::kmeans(input_dat[retval$is_in_region_filtered,],nclust,nstart = 30)
  cluster=rep(0,nrow(input_dat))
  cluster[retval$is_in_region_filtered]=clusters$cluster
  retval$cluster=cluster
  retval$num_clusters=nclust
  retval$cluster_frequency=nclust_all_repeats=K_hist
  return(retval)
}

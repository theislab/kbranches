#'Refilter already identified regions
#'
#'Refilters the results extracted by identify_regions().
#'Used only in rare cases where there is some noise remaining in the result of identify_regions().
#'
#'@param input_dat: data frame of input data with rows=samles and cols=dimensions.
#'@param sample_labels: vector of TRUE/FALSE labels, one for each sample. If sample_labels[i]==TRUE, then sample i is in region of interest (e.g. branching region)
#'@param smoothing_region: number of neighbours in to consider for label filtering
#'@param smoothing_region_thresh: minimum number of thresholds with same label in the neighbourhood required for the sample to keep it's label
#'@param Dist: matrix of sample to sample distances
#'@param dotsize1: size of points of class1 (red), used when plotting in 3D
#'@param dotsize2: size of points of class2 (green), used when plotting in 3D
#'
#'@return a list with elements:
#'\itemize{
#'  \item - is_in_region: a logical vector indicating which samples are part of the region (TRUE) and which not (FALSE)
#'  \item - is_in_region_filtered: same as 'is_in_region', but after S_neib filtering to reduce noise.
#'}
#'
#'@examples
#'this is an example
#'\dontrun{
#'
#'  example coming soon...
#'
#'}
#'@export
refilter_identified_regions=function(input_dat,sample_labels=NULL,smoothing_region=NULL,smoothing_region_thresh=NULL,Dist,dotsize1=7,dotsize2=7)#scale the number of neighbors depending on the local density)
{
  N=nrow(input_dat)#number of samples
  P=ncol(input_dat)#number of dimensions

  xlim=c(min(input_dat[,1]),max(input_dat[,1]))
  ylim=c(min(input_dat[,2]),max(input_dat[,2]))
  if(P==3){zlim=c(min(input_dat[,3]),max(input_dat[,3]))}else(zlim=NULL)#if there is a 3rd dimension

  tmp=sample_labels
  is_edge=sample_labels

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
        same_neighbors=sum(length(which(tmp_nn==tmp[i])))
        if(same_neighbors>=smoothing_region_thresh)#if above the threshold, the sample keeps it's label
        {
          is_edge_filtered[i]=tmp[i]
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

  if(is.null(smoothing_region==FALSE)&is.null(smoothing_region_thresh==FALSE))#smoothing_region and smoothing_region_thresh have been provided
  {
    if(P==3){open3d()}#open a new plot window
    #is_edge_filtered=filter_internal(smoothing_region,smoothing_region_thresh)
    #retval=list(is_edge=is_edge,is_edge_filtered=is_edge_filtered)
    retval=filter_internal(smoothing_region,smoothing_region_thresh)
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
  return(retval)
}

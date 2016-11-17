#'
#'Local K-Branch clustering for identifying regions of interest!
#'
#'Perform local clustering using kbranch.global in order to generate a GAP score for each sample.
#'The GAP score of each sample can be subsequantly used to identify regions of interest
#'such as tips of branches or branching regions
#'
#'@param input_dat: data frame of input data with rows=samles and cols=dimensions.
#'@param Dmat: matrix containing sample distances input_dat. Should be calculated using Dmat=compute_all_distances(input_dat)
#'@param S_neib: number of neighbours that defines the local neighbourhood. If set to NULL
#'a GUI helper will pop-up to assist in the selection of S_neib, providing a recommended value.
#'@param S_quant: the quantile of cumulative distance used to infer S_neib, the default value should usually suffice.
#'@param S_GUI_helper: If TRUE, a GUI helper will pop-up to aid in the selection of S_neib. If global_S is TRUE, then
#'the GUI will recommend a value for S_neib and visualize the neighbourhood for all data points. If global_S is false, the
#'GUI will only visualize the neighbourhood for every data point, since S_neib is selected automatically for each data point.
#'@param parallel_ncores: number of cores to use if parallel==TRUE. If set to NULL it will use the max number of
#'available cores. Defaults to NULL.
#'@param min_radius_quantile: the percentile to use for the identification of the 'median neighbourhood', defaults to 0.5 (median).
#'@param logfile: logfile to print output in the case of parallel computation.
#'@param nstart: number of initializations for clustering.
#'       Defaults to 5, increase it (e.g. to 10) if the results are too noisy.
#'@param nstart_GAP: number of initializations for clustering when calculating the GAP statistic.
#'       Defaults to 1, increase it (e.g. to 5 or 10) if the results are too noisy.
#'@param B_GAP: number of bootstrap datasets used to compute the GAP statistic.
#'       Defaults to 5, increase it (e.g. to 10 or 100) if the results are too noisy.
#'@param medoids: if TRUE, the medoids version of kbranch will be used (slower).
#'@param init_Kmeans: if TRUE, use K-Means for the initialization of K-haflines, otherwise use random initialization
#'
#'@return a list with elements:
#'\itemize{
#'  \item - gap_scores: list of the four different gap scores for each sample.
#'  \item - call: the call of the function
#'  \item - S_neib: the global value of S_neib used, or 'local' if global_S was false
#'}
#'
#'@examples
#' #this example might take some time to run
#' set.seed(1)
#'
#' #load the data, already in diffusion map format
#' data(scdata.loop.dmap)
#' input_dat <- scdata.loop.dmap[, 1:2] #keep the first 2 diffusion components
#'
#' #if the data are in not in diffusion space then
#' #performing diffusion map dimensionality reduction
#' #is necessary, for example:
#' #load(scdata.loop)
#' #dmap <- destiny::DiffusionMap(scdata.loop, sigma = 1000)
#' #input_dat <- destiny::as.data.frame(dmap)[, 1:2] #keep the first 2 diffusion components
#'
#' #compute the distances among all samples
#' Dmat <- compute_all_distances(input_dat)
#'
#' #perform local clustering to identify regions
#' #if you haven't specified the neighbourhood size S, it will be estimated
#' #set S_GUI_helper to FALSE for manual fine-tuning
#'
#' res <- kbranch.local(input_dat = input_dat, Dmat = Dmat)
#'
#' #identify regions of interest based on the GAP score
#' #of each sample computed by kbranch.local
#'
#' #If smoothing_region and smoothing_region_thresh are NULL, a GUI will
#' #pop-up to aid in their selection. Press OK to update the results.
#' #When you are happy with the filtering press 'x' to close the window.
#'
#' #smoothing: tip cell if at least 5 in 5 neighbors are tip cells
#' tips <- identify_regions(input_dat = input_dat, gap_scores = res$gap_scores, Dist = Dmat,
#'                          smoothing_region = 5, smoothing_region_thresh = 5, mode = 'tip')
#'
#' #plot the separate tips
#' plot(input_dat, pch=21, col = tips$cluster + 1, bg = tips$cluster + 1, main='tip regions')
#'
#' #smoothing: branching region cell if at least 10 in 10 neighbors are branching region cells
#' branch_reg <- identify_regions(input_dat = input_dat, gap_scores = res$gap_scores, Dist = Dmat,
#'                                smoothing_region = 10, smoothing_region_thresh = 10, mode='branch')
#'
#' #plot the branching_region(s)
#' plot(input_dat, pch = 21, col = branch_reg$cluster + 1, bg = branch_reg$cluster + 1, main = 'branching regions')
#'
#' ###################################################
#' #          end of example                         #
#' ###################################################
#'@export
#'
#'@importFrom rgl open3d plot3d points3d rgl.clear
#'@importFrom heplots arrow3d
#'@importFrom doParallel registerDoParallel
#'@importFrom parallel detectCores makeCluster
#'@importFrom foreach foreach getDoParWorkers %dopar%
#'@importFrom fgui gui guiGetValue guiSetValue
kbranch.local=function(input_dat,Dmat=NULL,S_neib=NULL,S_quant=0.1,S_GUI_helper=FALSE,parallel_ncores=NULL,min_radius_quantile=0.5,logfile='log.kbranch.local.txt',
                       nstart=5,nstart_GAP=1,B_GAP=5,medoids=FALSE,init_Kmeans=TRUE)
{
  # kbranch.local: identify regions based on the gap statistic produced by local clustering of k halflines

  ###these variables are required for compatibility with future releases###
  input_dat_orig=NULL
  global_S=TRUE
  #########################################################################

  input_dat=unique(input_dat)#to avoid errors later in the local clustering
  N=nrow(input_dat)
  P=ncol(input_dat)

  #print(paste(lsf.str()))

  #model selection will not work for a large number of diffusion components
  # if(P>3)
  # {
  #   print('Only using the first 3 diffusion components for model selection')
  #   P=3
  #   input_dat=input_dat[,1:P]
  #   Dmat=compute_all_distances(input_dat)
  # }

  cumul_dist = function(dst)
  {
    #calculates the cumulative distance for all neighbours
    #dst: vector of distance to all neighbours of a single sample
    #     it is the $dst field of the return value of find_K_nearest()
    #     for a given data point.
    retval=rep(0,length(dst))
    for( i in 1:length(dst))
    {
      retval[i]=sum(dst[1:i])
    }
    return(retval)
  }

  #cumulative distance used to suggest a value of S_neib
  #only necessary if S_neib was not provided by the user
  if(is.null(S_neib)==TRUE)
  {
    cdist=matrix(NA,N,N-1)
    #calculate the cumulative distance of each data point to all others
    for(i in 1:N)
    {
      neibs=find_K_nearest(n=i,K=N-1,Dist=Dmat)
      #cdist[i,j]: cumulative distance of point i to point j.
      #if point j is the Kst nearest neighbour (NN) of point i, then
      #cdist[i,j]=distance to j + distance to (K-1)st NN+...+dist to 1st NN of point i.
      cdist[i,]=cumul_dist(neibs$dst)
    }
  }

  # S_quant=0.1#quantile used to calculate the recommended S_neib

  if(global_S==TRUE)
  {
    #calculate the average density: mean distance of S_neib neighbors for each sample
    if(is.null(S_neib)==TRUE)#use a GUI helper to select S_neib
    {
      # print('got here!')
      #calculate a recommended value for the neighbourhood size S_neib
      #based on the average of all cumulative distances
      cdist_avg=apply(cdist,2,mean)

      #S_neib is the number of S_neib nearest neighbours that corerspond to the "S_quant"-th
      #percentile of average cumulative distance to all neighbours
      S_neib_recommended=which.min(abs(cdist_avg-quantile(cdist_avg,S_quant)))
      S_neib=S_neib_recommended #use the recommended value, unless a different value is selected using the GUI
      if(S_GUI_helper==TRUE)
      {
        #visualize the average cumulative distance
        x11();plot(cdist_avg,xlab='S_neib',ylab='Average cumulative distance')

        #gui_median_neighbourhood=function(S_neib,center.x=NULL,center.y=NULL,center.z=NULL)
        gui_median_neighbourhood=function(S_neib,minimum_radius=0.5,center)
        {
          center=center[1]
          # dens=rep(0,N)
          # dst=rep(0,N)
          # for(i in 1:N)
          # {
          #   # dens[i]=1/mean(find_K_nearest(n=i,K=S_neib,Dist=Dmat)$dst)
          #   dst[i]=mean(find_K_nearest(n=i,K=S_neib,Dist=Dmat)$dst)
          #   dens[i]=1/dst[i]
          #   # dens[i]=max(dst[i])
          # }
          # #find the threshold for the minimum distance
          # pos_median=which(abs(dens-quantile(dens,min_radius_quantile))==min(abs(dens-quantile(dens,min_radius_quantile))))[1]
          # min_radius=max(find_K_nearest(n=pos_median,K=S_neib,Dist=Dmat)$dst)

          S_neib_radius=rep(0,N)#neighbourhood radius for each sample, used to calculate min neighbourhood size for scaling
          for(i in 1:N)
          {
            S_neib_radius[i]=max(find_K_nearest(n=i,K=S_neib,Dist=Dmat)$dst)
          }

          # S_neib_radius=rep(0,N)#neighbourhood radius for each sample
          # for(i in 1:N)
          # {
          #   S_neib_radius[i]=max(find_K_nearest(n=i,K=S_neib,Dist=Dmat)$dst)
          # }
          #find the threshold for the minimum distance
          #pos_median=1
          # min_radius_quantile=min_radius_quantile
          min_radius_quantile_internal=minimum_radius
          min_radius=quantile(S_neib_radius,min_radius_quantile_internal)#if a neighbourhood radius is smaller than min_radius, it is scaled up to min_radius before clustering

          # #for VISUALIZATION ONLY:
          # #check whether to use a sample (other than the median) as the center of the neighbourhood
          # if((P==2)&(is.null(center.x)==FALSE)&(is.null(center.y)==FALSE))
          # {
          #   #center the neighbourhood on the sample closest to the specified x,y coords, instead of the median
          #   pos_median=which.min(t(apply(input_dat,1,function(x){norm(x-c(center.x,center.y),'2')^2})))
          #
          # }else if((P==3)&(is.null(center.x)==FALSE)&(is.null(center.y)==FALSE)&(is.null(center.z)==FALSE))
          # {
          #   #center the neighbourhood on the sample closest to the specified x,y coords, instead of the median
          #   pos_median=which.min(t(apply(input_dat,1,function(x){norm(x-c(center.x,center.y,center.z),'2')^2})))
          # }

          #get all the samples in the neighbourhood
          nnres=find_K_nearest(n=center,K=S_neib,Dist=Dmat)
          nn=sort(nnres$pos)
          #scale the number of nearest neighbors in dense regions
          Kscaled=scale_S_neib(input_dat=input_dat,n=center,K_init=S_neib,Dist=Dmat,dmin=min_radius)
          nnres=find_K_nearest(n=center,K=Kscaled,Dist=Dmat)
          nn_scaled=sort(nnres$pos)
          nn_extra=setdiff(nn_scaled,nn)

          red_paint=colorRampPalette(RColorBrewer::brewer.pal(3,'Reds'))(3)[3]#colorblind_friendly red
          green_paint=colorRampPalette(RColorBrewer::brewer.pal(3,'Greens'))(3)[3]#colorblind_friendly green

          if(P==2)
          {
            plot(input_dat,main=paste('Median neighbourhood, recommended S_neib:',S_neib_recommended))
            #points(x=input_dat[pos_median,1],y=input_dat[pos_median,2],col='magenta',pch=23,bg='magenta',cex=2)
            points(x=input_dat[nn,1],y=input_dat[nn,2],col=green_paint,pch=21,bg=green_paint)
            if(length(nn_extra)>0){points(x=input_dat[nn_extra,1],y=input_dat[nn_extra,2],col=red_paint,pch=21,bg=red_paint)}#scaled neighbourhood
            text(x=input_dat[center,1],y=input_dat[center,2],col='magenta',labels = 'X',cex=2)
          }else if(P==3)
          {
            plot3d(x=input_dat[,1],y=input_dat[,2],z=input_dat[,3],xlab='x axis',ylab='y axis',zlab='z axis',size=7,main=paste('Median neighbourhood, recommended S_neib:',S_neib_recommended))
            points3d(x=input_dat[nn,1],y=input_dat[nn,2],z=input_dat[nn,3],size=10,col=green_paint)#plot the neighbourhood of the median
            if(length(nn_extra)>0){points3d(x=input_dat[nn_extra,1],y=input_dat[nn_extra,2],z=input_dat[nn_extra,3],size=10,col=red_paint)}#scaled neighbourhood
            #points3d(x=input_dat[pos_median,1],y=input_dat[pos_median,2],z=input_dat[pos_median,3],size=15,col='magenta')#plot center
            rgl::rgl.texts(x=input_dat[center,1],y=input_dat[center,2],z=input_dat[center,3],size=25,col='magenta',text='X')
          }
          #return(paste(' S_neib =',S_neib,'\n center.x=',center.x,'\n center.y=',center.y,'\n center.z=',center.z))
        }

        #{x11()}else if(P==3){open3d()}#open plotting windows for the GUI
        if(P==2)
        {
          x11()
          plot(input_dat,main='Select points to visualize neighbourhood')
          pts=identify(input_dat)
          dev.off();x11();#close the plotting window and start a new one
          argSlider = list(S_neib=c(1,N-1,1),minimum_radius=c(0.1,1,0.1))
          argList = list(center=pts) #point(s) to center the neighbourhood
          gui_retval = gui(gui_median_neighbourhood,argSlider = argSlider, argList=argList,
                             title='Select local neighbourhood size: S_neib',output = NULL,closeOnExec = FALSE)
          # argSlider = list(S_neib=c(1,N-1,1),
          #                  center.x=c(min(input_dat[,1]),max(input_dat[,1]),1/N),#1/N step size
          #                  center.y=c(min(input_dat[,2]),max(input_dat[,2]),1/N),
          #                  center.z=c(0,0,0))
          # S_neib_selected = gui(gui_median_neighbourhood,argSlider = argSlider,title='Select local neighbourhood size: S_neib',output = NULL,closeOnExec = FALSE)


        }else if(P==3)
        {
          open3d()
          plot3d(x=input_dat[,1],y=input_dat[,2],z=input_dat[,3],xlab='x axis',ylab='y axis',zlab='z axis',size=7,main='Select points to visualize neighbourhood')
          pts=rgl::identify3d(x=input_dat[,1],y=input_dat[,2],z=input_dat[,3])
          rgl.clear()#clear the plot
          argSlider = list(S_neib=c(1,N-1,1),minimum_radius=c(0.1,1,0.1))
          argList = list(center=pts) #point(s) to center the neighbourhood
          gui_retval = gui(gui_median_neighbourhood,argSlider = argSlider, argList=argList,
                             title='Select local neighbourhood size: S_neib',output = NULL,closeOnExec = FALSE)
          # argSlider = list(S_neib=c(1,N-1,1),
          #                  center.x=c(min(input_dat[,1]),max(input_dat[,1]),1/N),
          #                  center.y=c(min(input_dat[,2]),max(input_dat[,2]),1/N),
          #                  center.z=c(min(input_dat[,3]),max(input_dat[,3]),1/N))
          # S_neib_selected = gui(gui_median_neighbourhood,argSlider = argSlider,title='Select local neighbourhood size: S_neib',output = NULL,closeOnExec = FALSE)
        }
        #S_neib_selected = gui(gui_median_neighbourhood,argSlider = argSlider,title='Select local neighbourhood size: S_neib',output = NULL,closeOnExec = FALSE)
        S_neib=gui_retval$S_neib
        min_radius_quantile=gui_retval$minimum_radius
      }
    }
    S_neib_radius=rep(0,N)#neighbourhood radius for each sample, used to calculate min neighbourhood size for scaling
    for(i in 1:N)
    {
      S_neib_radius[i]=max(find_K_nearest(n=i,K=S_neib,Dist=Dmat)$dst)
    }
    min_radius=quantile(S_neib_radius,min_radius_quantile)
    # #calculate the median neighbourhood distance
    # dens=rep(0,N)
    # dst=rep(0,N)
    # for(i in 1:N)
    # {
    #   # dens[i]=1/mean(find_K_nearest(n=i,K=S_neib,Dist=Dmat)$dst)
    #   dst[i]=mean(find_K_nearest(n=i,K=S_neib,Dist=Dmat)$dst)
    #   dens[i]=1/dst[i]
    #   # dens[i]=max(dst[i])
    # }
    # #find the threshold for the minimum distance
    # pos_median=which(abs(dens-quantile(dens,min_radius_quantile))==min(abs(dens-quantile(dens,min_radius_quantile))))[1]
    # min_radius=max(find_K_nearest(n=pos_median,K=S_neib,Dist=Dmat)$dst)

    # min_radius=max(find_K_nearest(n=pos_median,K=S_neib,Dist=Dmat)$dst)
  }else#global_S==FALSE
  {
    #use a GUI to visualize the neighbourhood
    if(S_GUI_helper==TRUE)
    {
      gui_local_neighbourhood=function(center.x=NULL,center.y=NULL,center.z=NULL)
      {
        #for VISUALIZATION ONLY:
        #check whether to use a sample (other than the median) as the center of the neighbourhood
        if((P==2)&(is.null(center.x)==FALSE)&(is.null(center.y)==FALSE))
        {
          #center the neighbourhood on the sample closest to the specified x,y coords, instead of the median
          pos_median=which.min(t(apply(input_dat,1,function(x){norm(x-c(center.x,center.y),'2')^2})))

        }else if((P==3)&(is.null(center.x)==FALSE)&(is.null(center.y)==FALSE)&(is.null(center.z)==FALSE))
        {
          #center the neighbourhood on the sample closest to the specified x,y coords, instead of the median
          pos_median=which.min(t(apply(input_dat,1,function(x){norm(x-c(center.x,center.y,center.z),'2')^2})))
        }

        #get all the samples in the neighbourhood
        S_neib_recommended=which.min(abs(cdist[pos_median,]-quantile(cdist[pos_median,],S_quant)))#local S_neib of given data point
        nnres=find_K_nearest(n=pos_median,K=S_neib_recommended,Dist=Dmat)
        nn=sort(nnres$pos)
        if(P==2)
        {
          plot(input_dat,main=paste('Local neighbourhood, recommended S_neib:',S_neib_recommended))
          points(x=input_dat[pos_median,1],y=input_dat[pos_median,2],col='red',pch=23,bg='magenta',cex=2)
          points(x=input_dat[nn,1],y=input_dat[nn,2],col='green',pch=21,bg='green')
        }else if(P==3)
        {
          plot3d(x=input_dat[,1],y=input_dat[,2],z=input_dat[,3],xlab='x axis',ylab='y axis',zlab='z axis',size=7,main=paste('Local neighbourhood, recommended S_neib:',S_neib_recommended))
          points3d(x=input_dat[nn,1],y=input_dat[nn,2],z=input_dat[nn,3],size=10,col='green')#plot the neighbourhood of the median
          points3d(x=input_dat[pos_median,1],y=input_dat[pos_median,2],z=input_dat[pos_median,3],size=15,col='magenta')#plot median point
        }
        return(paste(' S_neib =',S_neib,'\n center.x=',center.x,'\n center.y=',center.y,'\n center.z=',center.z))
      }

      #{x11()}else if(P==3){open3d()}#open plotting windows for the GUI
      if(P==2)
      {
        x11()
        argSlider = list(S_neib=c(1,N-1,1),
                         center.x=c(min(input_dat[,1]),max(input_dat[,1]),1/N),#1/N step size
                         center.y=c(min(input_dat[,2]),max(input_dat[,2]),1/N),
                         center.z=c(0,0,0))
        #argSlider = list(S_neib=c(1,N-1,1),center.x=c(min(input_dat[,1]),max(input_dat[,1]),1),center.y=c(min(input_dat[,2]),max(input_dat[,2]),1))
        #argText=list(center.z='unavailable')
        S_neib_selected = gui(gui_local_neighbourhood,argSlider = argSlider,title='local neighbourhood',output = NULL,closeOnExec = FALSE)
      }else if(P==3)
      {
        open3d()
        argSlider = list(S_neib=c(1,N-1,1),
                         center.x=c(min(input_dat[,1]),max(input_dat[,1]),1/N),
                         center.y=c(min(input_dat[,2]),max(input_dat[,2]),1/N),
                         center.z=c(min(input_dat[,3]),max(input_dat[,3]),1/N))
        S_neib_selected = gui(gui_local_neighbourhood,argSlider = argSlider,title='local neighbourhood',output = NULL,closeOnExec = FALSE)
      }
    }
  }

  K_all=c(1,2,3)
  gscore=matrix(data=NA,nrow = N,ncol = length(K_all))
  gscorel=matrix(data=NA,nrow = N,ncol = length(K_all))
  gscore_orig=matrix(data=NA,nrow = N,ncol = length(K_all))
  gscore_orig_nl=matrix(data=NA,nrow = N,ncol = length(K_all))
  #clr=rep(0,N)

  #if(global_S==T){print(paste('global S_neib:',S_neib))}else{print('local S_neib used')}
  print(paste('S_neib:',S_neib,', min_radius_quantile: ',min_radius_quantile))
  # if(parallel==TRUE)#if parallel, register a parallel backend
  # {
    #set the backend for parallel computation
    numcores=parallel::detectCores()#see how many cores there are
    if(is.null(parallel_ncores))#use the maximum number of available cores
    {
      cl = makeCluster(numcores-1)#create the cluster, -1 for stability
      print(paste('Parrallel execution using max num of ',numcores,'-1 =',numcores-1,'cores (-1 for stability)'))
    }else# a number of cores has been specified
    {
      if((parallel_ncores>=numcores)|(parallel_ncores<1))#can't use more than the maximum, or less than 1
      {
        cl = makeCluster(numcores-1)#create the cluster, -1 for stability
        print(paste('Parrallel execution using max num of ',numcores,'-1 =',numcores-1,'cores (-1 for stability)'))
      }else
      {
        cl = makeCluster(parallel_ncores)#create the cluster, -1 for stability
        print(paste('Using',parallel_ncores,'core(s)'))
      }
    }
    print(paste('refer to',logfile,'for progress'))
    doParallel::registerDoParallel(cl)#register the cluster
    # doParallel::registerDoParallel(cores=numcores-1)#register the cluster
    #print(getDoParWorkers())#see how many workers there are
  # }

  writeLines(c(''), logfile)
  #res_foreach <- foreach(i=1:N,.export = curr_env) %dopar%
  res_foreach <- foreach(i=1:N,.packages = 'kbranch') %dopar%
  {
    #print the progress in a logfile
    sink(logfile,append = T)
    print(paste('sample',i,'of',N))
    sink()
    if(global_S==T)#upscale the S_neib number in high density regions
    {
      Kscaled=scale_S_neib(input_dat=input_dat,n=i,K_init=S_neib,Dist=Dmat,dmin=min_radius)#scale the number of nearest neighbors
    }else
    {
      Kscaled=which.min(abs(cdist[i,]-quantile(cdist[i,],S_quant)))
    }
    if(Kscaled<3){Kscaled=3}#Cannot be less than the number of clusters
    # print(Kscaled)
    nnres=find_K_nearest(n=i,K=Kscaled,Dist=Dmat)
    nn=sort(nnres$pos)
    if(is.null(input_dat_orig))#cluster in diffusion map dimensions
    {
      local_dat=rbind(input_dat[i,],input_dat[nn,])
      D_local=Dmat[c(i,nn),c(i,nn)]

      tmp_gap=rep(NA,length(K_all))
      tmp_gapl=rep(NA,length(K_all))
      tmp_gap_orig=rep(NA,length(K_all))
      tmp_gap_orig_nl=rep(NA,length(K_all))
      for(k in 1:length(K_all))
      {
        Kappa=K_all[k]
        clust=kbranch.global(Kappa=Kappa,input_dat=local_dat,c0=NULL,Vmat=NULL,silent=T,nstart=nstart,nstart_GAP=nstart_GAP,B_GAP=B_GAP,fixed_center=1,medoids=medoids,Dmat=D_local,init_Kmeans = init_Kmeans)
        tmp_gap[k]=clust$GAP
        tmp_gapl[k]=clust$GAPl
        tmp_gap_orig[k]=clust$GAP_orig
        tmp_gap_orig_nl[k]=clust$GAP_orig_no_log
      }
    }else#cluster in the original Dimensions
    {
      local_dat_orig=rbind(input_dat_orig[i,],input_dat_orig[nn,])#use the original space for clustering
      D_local_orig=Dmat_orig[c(i,nn),c(i,nn)]
      local_dat=rbind(input_dat[i,],input_dat[nn,])#use the dmap space for the gap statistic
      D_local=Dmat[c(i,nn),c(i,nn)]

      tmp_gap=rep(NA,length(K_all))
      tmp_gapl=rep(NA,length(K_all))
      tmp_gap_orig=rep(NA,length(K_all))
      tmp_gap_orig_nl=rep(NA,length(K_all))
      for(k in 1:length(K_all))
      {
        Kappa=K_all[k]
        clust=kbranch.global(Kappa=Kappa,input_dat=local_dat_orig,c0=NULL,Vmat=NULL,silent=T,nstart=nstart,nstart_GAP=NULL,B_GAP=NULL,fixed_center=1,medoids=medoids,Dmat=D_local_orig,init_Kmeans = init_Kmeans)

        #calculate the class centers
        centers=matrix(NA,nrow = Kappa,ncol = P)
        for(k in 1:Kappa){centers[k,]=colMeans(local_dat[clust$cluster==k,,drop=FALSE])}

        #fist estimate the kbranch.global error in dmap space, using the cluster labels from the original space
        #Vmat is equal to the centers matrix calculated above, centered at the data point being examined
        c0_local=as.numeric(local_dat[1,])#c0 is always row 1 in the local data
        Vmat_local=t(t(centers)-c0_local)
        kbranch.global_err=0#clustering error
        for(i in 1:nrow(local_dat))#compute the error (distance from halfline object) for every data point
        {
          pos=clust$cluster[i]#select the direction vector according to the class of the data point
          v=Vmat_local[pos,]

          x=as.numeric(local_dat[i,])#otherwise matrix operations won't work
          x_hat=x-c0_local#take the C1 vector as the starting point of the axes of the vector space
          if(sign(as.numeric(x_hat%*%v))>=0)#check dot product positive -> must be close to half-line, not line
          {
            kbranch.global_err=kbranch.global_err+(norm(as.matrix((x-c0_local)-((v%o%v)/as.numeric(v%*%v))%*%(x-c0_local)),'2')^2)#line_dist
          }else#point on the wrong side, measure distance from the center
          {
            kbranch.global_err=kbranch.global_err+(norm(as.matrix(x-c0_local),'2')^2)#point_dist
          }
        }

        #calculate the GAP
        resGAP=calculate_GAP(input_dat=local_dat,clustering_error=kbranch.global_err,Kappa=Kappa,cluster_labels=clust$cluster,
                             nstart_GAP=nstart_GAP,B_GAP=B_GAP,fixed_center=1,medoids=medoids,Dmat=D_local,init_Kmeans = init_Kmeans)
        tmp_gap[k]=resGAP$GAP
        tmp_gapl[k]=resGAP$GAPl
        tmp_gap_orig[k]=resGAP$GAP_orig
        tmp_gap_orig_nl=resGAP$GAP_orig_no_log

      }
    }

    #this list is 'saved' for each iteration of the foreach loop
    lst=list(gscore=tmp_gap,gscorel=tmp_gapl,gscore_orig=tmp_gap_orig,gscore_orig_nl=tmp_gap_orig_nl)

  }
  # doParallel::stopImplicitCluster()#stop the cluster in case something went wrong
  parallel::stopCluster(cl)
  for(i in 1:N)
  {
    gscore[i,]=res_foreach[[i]]$gscore
    gscorel[i,]=res_foreach[[i]]$gscorel
    gscore_orig[i,]=res_foreach[[i]]$gscore_orig
    gscore_orig_nl[i,]=res_foreach[[i]]$gscore_orig_nl
  }
  retval=list()
  retval$gap_scores=list(gscore=gscore,gscorel=gscorel,gscore_orig=gscore_orig,gscore_orig_nl=gscore_orig_nl)
  retval$call=match.call()
  if(global_S==T){retval$S_neib=S_neib}else{retval$S_neib='local'}

  return(retval)
}

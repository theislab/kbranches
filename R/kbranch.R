#'
#'Local K-Star clustering for identifying regions of interest!
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
#'@param medoids: if TRUE, the medoids version of khalflines will be used (slower).
#'@param init_Kmeans: if TRUE, use K-Means for the initialization of K-haflines, otherwise use random initialization
#'
#'@return a list with elements:
#'\itemize
#'{
#'  \item - gap_scores: list of the four different gap scores for each sample.
#'  \item - call: the call of the function
#'  \item - S_neib: the global value of S_neib used, or 'local' if global_S was false
#'}
#'
#'@examples
#'  #this example might take some time to run
#'  set.seed(1)
#'
#'  #load the data, already in diffusion map format
#'  input_dat=scdata.loop.map[,1:2] #keep the first 2 diffusion components
#'
#'  #if the data are in not in diffusion map then
#'  #performing diffusion map dimensionality reduction
#'  #is necessary, for example:
#'  #raw_dat=scdata.loop
#'  #dmap=destiny::DiffusionMap(raw_dat,sigma = 1000)
#'  #input_dat=destiny::as.data.frame(dmap)[,1:2] #keep the first 2 diffusion components
#'
#'  #compute the distances among all samples
#'  Dmat=compute_all_distances(input_dat)
#'
#'  #perform local clustering to identify regions
#'  #if you haven't specified the neighbourhood size S, it will be estimated
#'  #set S_GUI_helper to FALSE for manual fine-tuning
#'
#'  res=kbranch.local(input_dat=input_dat,Dmat=Dmat)
#'
#'  #identify regions of interest based on the GAP score
#'  #of each sample computed by kbranch.local
#'
#'  #If smoothing_region and smoothing_region_thresh are NULL, a GUI will
#'  #pop-up to aid in their selection. Press OK to update the results.
#'  #When you are happy with the filtering press 'x' to close the window.
#'
#'  #smoothing: tip cell if at least 2 in 5 neighbors are tip cells
#'  tips=identify_regions(input_dat=input_dat,gap_scores=res$gap_scores,Dist=Dmat,
#'                        smoothing_region = 5,smoothing_region_thresh = 2,mode='tip')
#'
#'  #plot the separate tips
#'  plot(input_dat,pch=21,col=tips$cluster+1,bg=tips$cluster+1,main='tip regions')
#'
#'  #smoothing: branching region cell if at least 5 in 5 neighbors are branching region cells
#'  branch_reg=identify_regions(input_dat=input_dat,gap_scores=res$gap_scores,
#'              Dist=Dmat,smoothing_region = 5,smoothing_region_thresh = 5,mode='branch')
#'
#'  #plot the branching_region(s)
#'  plot(input_dat,pch=21,col=branch_reg$cluster+1,bg=branch_reg$cluster+1,main='branching regions')
#'\dontrun
#'{
#'  ###################################################
#'  #          end of example                         #
#'  ###################################################
#'
#'
#'}
#'@export
#'
#'@importFrom rgl open3d plot3d points3d rgl.clear
#'@importFrom heplots arrow3d
#'@importFrom doParallel registerDoParallel
#'@importFrom parallel detectCores makeCluster
#'@importFrom foreach foreach getDoParWorkers %dopar%
#'@importFrom fgui gui guiGetValue guiSetValue
kbranch.local=function(input_dat,Dmat=NULL,S_neib=NULL,S_quant=0.1,S_GUI_helper=FALSE,parallel_ncores=NULL,min_radius_quantile=0.5,logfile='log.kbranch.global.txt',
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
  res_foreach <- foreach(i=1:N,.packages = 'khalflines') %dopar%
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
#'\itemize
#'{
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
#'\itemize
#'{
#'  \item - is_in_region: a logical vector indicating which samples are part of the region (TRUE) and which not (FALSE)
#'  \item - is_in_region_filtered: same as 'is_in_region', but after S_neib filtering to reduce noise.
#'}
#'
#'@examples
#'this is an example
#'\dontrun
#'{
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
#'\dontrun
#'{
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


#'
#'Clustering on K-Halflines
#'
#'Clusters data on K-Halflines with a common center
#'and calculates the corresponding GAP statistic
#'
#'@param input_dat: data frame of input data with rows=samles and cols=dimensions.
#'@param Kappa: number of clusters (halflines)
#'@param Dmat: matrix containing sample distances
#'@param init_Kmeans: if TRUE: initialize directions v1,...,vk using K-Means. FALSE: use directions of randomly selected samples
#'@param c0: initial value for the center of all half-lines
#'@param Vmat: matrix whose K rows are the direction vectors
#'@param nstart_GAP: number of initializations for clustering when calculating the GAP statistic
#'@param nstart_kmeans: number of initializations for Kmeans (when using Kmeans to initialize khalflines)
#'@param B_GAP: number of bootstrap datasets used to compute the GAP statistic, if NULL (default), it won't be computed
#'@param fixed_center: if not NULL, then K-halflines will run with the given center fixed
#'@param medoids: if TRUE, the medoids version of khalflines will be used (slower)
#'@param silent: set to FALSE to display messages (for debugging)
#'@param silent_internal: set to TRUE to display messages and plots of internal clustering functions (for debugging)
#'@param show_plots: if TRUE, the clustering will result be plotted
#'@param show_lines: if TRUE, show the halflines in the plot
#'@param show_plots_GAP: if TRUE, show the plots when performing clustering under the null distribution to calculate the GAP statistic (for debugging)
#'
#'@return a list with elements:
#'\itemize
#'{
#'  \item - cluster: cluster assignment for each sample (numeric)
#'  \item - Kappa: number of clusters (halflines)
#'  \item - err: total clustering cost
#'  \item - iters: total iterations of the algorithm
#'  \item - c0: position (row index in input_dat) of the center sample
#'  \item - Vmat: positions (row indices in input_dat) of the direction samples
#'  \item - clust_counts: number (count) of samples in each of the clusters
#'  \item - all_clustering_errors: vector of total clustering error for each of the nstart different initializations
#'  \item - all_clusterings: total results for each of the nstart different initializations
#'  \item - GAP: value of the modified GAP statistic for the given Kappa
#'  \item - GAPl: value of the modified GAP statistic for the given Kappa using the logarithm of the expected dispersion
#'  \item - GAP_orig: value of the oroginal GAP statistic for the given Kappa (using the logarithm of the expected dispersion)
#'  \item - GAP_orig_no_log: value of the oroginal GAP statistic for the given Kappa (without using the logarithm of the expected dispersion)
#'  \item - GAP.sd: standard deviation of GAP
#'  \item - GAPl.sd: standard deviation of GAPl
#'  \item - GAP_orig.sd: standard deviation of GAP_orig
#'  \item - GAP_orig_no_log.sd: standard deviation of GAP_orig_no_log
#'  \item - call: function call
#'}
#'
#'@examples
#'  #cluster the 2D data on three haflines
#'  set.seed(1)
#'
#'  #laod the data
#'  raw_dat=scdata.3lines.simulated6genes_subsampled
#'
#'  #perform diffusion map dimensionality reduction
#'  dmap=destiny::DiffusionMap(raw_dat,sigma = 1000)
#'
#'  #keep the first 2 diffusion components
#'  input_dat=destiny::as.data.frame(dmap)[,1:2]
#'
#'  #cluster into a K-Star with K=3
#'  clust=kbranch.global(input_dat,Kappa=3)
#'
#'  #plot the clustering results
#'  plot(input_dat,pch=21,col=clust$cluster,bg=clust$cluster,main='K-Star clustering')
#'
#'
#'@export
kbranch.global=function(input_dat,Kappa,Dmat=NULL,init_Kmeans=TRUE,c0=NULL,Vmat=NULL,nstart=20,nstart_GAP=20,nstart_kmeans=20,B_GAP=NULL,
                 fixed_center=NULL,medoids=FALSE,silent=TRUE,silent_internal=TRUE,show_plots=FALSE,show_lines = TRUE,show_plots_GAP=FALSE)
{
  #   external function for clustering k halflines
  #   serves as wrapper to kbranch.global_internal and also initializes c0 and Vmat if values have not been provided
  #   and calculates the GAP statistic for the given number of clusters
  #
  #   inputs:
  #     -input_dat: data frame of input data with rows=samles and cols=dimensions
  #     -Kappa: number of clusters (halflines)
  #     -Dmat: matrix containing sample distances
  #     -init_Kmeans: if TRUE: initialize directions v1,...,vk using K-Means. FALSE: use directions of randomly selected samples
  #     -c0: initial value for the center of all half-lines
  #     -Vmat: matrix whose K rows are the direction vectors
  #     -nstart_GAP: number of initializations for clustering when calculating the GAP statistic
  #     -nstart_kmeans: number of initializations for Kmeans (when using Kmeans to initialize khalflines)
  #     -B_GAP: number of bootstrap datasets used to compute the GAP statistic
  #     -fixed_center: if not NULL, then K-halflines will run with the given center fixed
  #     -medoids: if FALSE, the medoids version of khalflines will be used (slower)
  #     -silent: set to FALSE to display messages (for debugging)
  #     -silent_internal: set to TRUE to display messages and plots of internal clustering functions (for debugging)
  #     -show_plots: if TRUE, the clustering will result be plotted
  #     -show_lines: if TRUE, show the halflines in the plot
  #     -show_plots_GAP: if TRUE, show the plots when performing clustering under the null distribution to calculate the GAP statistic (for debugging)
  #
  #   return values:
  #     -cluster: cluster assignment for each sample (numeric)
  #     -Kappa: number of clusters (halflines)
  #     -err: total clustering cost
  #     -iters: total iterations of the algorithm
  #     -c0: position (row index in input_dat) of the center sample
  #     -Vmat: positions (row indices in input_dat) of the direction samples
  #     -clust_counts: number (count) of samples in each of the clusters
  #     -all_clustering_errors: vector of total clustering error for each of the nstart different initializations
  #     -all_clusterings: total results for each of the nstart different initializations
  #     -GAP: value of the modified GAP statistic for the given Kappa
  #     -GAPl: value of the modified GAP statistic for the given Kappa using the logarithm of the expected dispersion
  #     -GAP_orig: value of the oroginal GAP statistic for the given Kappa (using the logarithm of the expected dispersion)
  #     -GAP_orig_no_log: value of the oroginal GAP statistic for the given Kappa (without using the logarithm of the expected dispersion)
  #     -GAP.sd: standard deviation of GAP
  #     -GAPl.sd: standard deviation of GAPl
  #     -GAP_orig.sd: standard deviation of GAP_orig
  #     -GAP_orig_no_log.sd: standard deviation of GAP_orig_no_log
  #     -call: function call

  #######################################################
  # Internal functions of kbranch.global
  #######################################################

  line_dist=function(x,c0,v_outer=NULL)
  {
    #   calculate point-line squared Euclidian distance
    #
    #   inputs:
    #     -x: numeric vector
    #     -c0: numeric vector
    #
    #   return values:
    #     -point-line distance

    #d=(x-c0)-((v%o%v)/as.numeric(v%*%v))%*%(x-c0)
    x_hat=x-c0
    # if(is.null(v_outer))
    # {
    #   d=x_hat-(v%o%v)%*%x_hat#no need to normalize since v already as unit length
    # }else
    # {
    #   d=x_hat-v_outer%*%x_hat#no need to normalize since v already as unit length
    # }
    d=x_hat-v_outer%*%x_hat#no need to normalize since v already as unit length
    return(norm(as.matrix(d),'2')^2) #return the 2-norm (Euclidian) squared
  }

  point_dist=function(x,c0)
  {
    #   calculate point-point squared Euclidian distance
    #
    #   inputs:
    #     -x: numeric vector
    #     -c0: numeric vector
    #
    #   return values:
    #     -point-point distance

    return(norm(as.matrix(x-c0),'2')^2) #return the 2-norm (Euclidian) squared
  }

  assign_to_closest_hline=function(input_dat,Kappa,c0,Vmat)#assign data points to the closest half line
  {
    #   assign each sample to it's closest cluster (halfline)
    #
    #   inputs:
    #     -input_dat: data frame of input data with rows=samles and cols=dimensions
    #     -Kappa: number of clusters (halflines)
    #     -c0: current value for the center of all half-lines
    #     -Vmat: matrix whose K rows are the current direction vectors
    #
    #   return values:
    #     -dclass: data frame of class label for each sample

    xdata_class=data.frame(class=rep(0,nrow(input_dat)),wrong_x=rep(NA,nrow(input_dat)),wrong_y=rep(NA,nrow(input_dat)))
    for(i in 1:nrow(input_dat))
    {
      x=as.numeric(input_dat[i,])#otherwise matrix operations won't work
      nc=-1;#which is the nearest class
      nc_dist=Inf
      x_hat=x-c0#take the C1 vector as the starting point of the axes of the vector space
      # print('Vmat in assign, before')
      # print(Vmat)
      for(k in 1:Kappa)
      {
        v_k=Vmat[k,]#current direction vector
        v_outer=v_k%o%v_k
        if(anyNA(Vmat)){
          print(paste('x_hat',as.numeric(x_hat)))
          print(paste('v_k',as.numeric(v_k)))
          print(Vmat)}

        if(sign(as.numeric(x_hat%*%v_k))>=0)#check dot product positive -> must be close to half-line, not line
        {
          #nc_dist_new=line_dist(x,c0,v_k)
          nc_dist_new=line_dist(x,c0,v_outer=v_outer)
          if(nc_dist_new<nc_dist)
          {
            nc=k;
            nc_dist=nc_dist_new
          }
        }else#point on the wrong side, measure distance from the center
        {
          nc_dist_new=point_dist(x,c0)
          if(nc_dist_new<nc_dist)
          {
            nc=k;
            nc_dist=nc_dist_new
          }
        }
      }

      xdata_class$class[i]=nc;#assign closest cluster
      # View(xdata_class)
      if(nc==-1)#none of the three clusters was selected, all had negative dot-products
      {
        print(paste('------- ERROR! vector with x=',x[1],' y=',x[2],'-------'))
      }
    }

    return(dclass=xdata_class) #return the class for each sample
  }

  plot_hline_clust2D=function(input_dat,dclass,c0,Vmat,plot_title,show_lines=TRUE,len=10^10,xlim=NULL,ylim=NULL,zlim=NULL)
  {
    #plot the result of the k halflines method
    #note, if the halflines are not visible or short, increase len
    #input_dat: input data
    #c0: center of all the half-lines
    #Vmat: matrix of direction vectors
    #dclass: vector of cluster assignments (cluster labels)
    #add_on: add to existing plot

    #   plot the result of the k halflines method in 2D
    #   note, if the halflines are not visible or short, increase len
    #
    #   inputs:
    #     -input_dat: data frame of input data with rows=samles and cols=dimensions
    #     -dclass: data frame of class label for each sample
    #     -c0: current value for the center of all half-lines
    #     -Vmat: matrix whose K rows are the current direction vectors
    #     -plot_title: string, title of the plot
    #     -show_lines: if TRUE, show the halflines on the plot, if they still aren't visible increase len
    #     -len: length of the halflines
    #     -xlim: x axlis limits
    #     -ylim: y axlis limits
    #     -zlim: z axlis limits
    #
    #   return values:
    #     -no return values

    # plot.new()
    plot(x=input_dat[,1],y=input_dat[,2],col=dclass$class+1,xlab = 'x', ylab = 'y',xlim=xlim,ylim=ylim,zlim=zlim) #plot the data points
    title(plot_title)

    if(show_lines==TRUE)
    {
      points(x = c0[1], y = c0[2],pch=3, cex=2, col='blue3') #plot the center of all clusters
      K=length(unique(dclass$class))
      for(k in 1:K)#print all half-lines of length len
      {
        v1=as.numeric(Vmat[k,])
        points(x=c(c0[1],len*v1[1]+c0[1]),y=c(c0[2],len*v1[2]+c0[2]),type = 'l', col=k+1) #direction vectors v_i starting from c0
      }
      return()
    }
  }

  plot_hline_clust3D=function(input_dat,dclass,c0,Vmat,plot_title,show_lines=TRUE,len=10^10,xlim=NULL,ylim=NULL,zlim=NULL)#to plot data with 3 centers c1,c2,c3 instead of c0
  {
    #   plot the result of the k halflines method in 3D
    #   note, if the halflines are not visible or short, increase len
    #
    #   inputs:
    #     -input_dat: data frame of input data with rows=samles and cols=dimensions
    #     -dclass: data frame of class label for each sample
    #     -c0: current value for the center of all half-lines
    #     -Vmat: matrix whose K rows are the current direction vectors
    #     -plot_title: string, title of the plot
    #     -show_lines: if TRUE, show the halflines on the plot, if they still aren't visible increase len
    #     -len: length of the halflines
    #     -xlim: x axlis limits
    #     -ylim: y axlis limits
    #     -zlim: z axlis limits
    #
    #   return values:
    #     -no return values

    open3d()
    plot3d(x=input_dat[,1],y=input_dat[,2],z=input_dat[,3],xlab='x axis',ylab='y axis',zlab='z axis',col=dclass$class+1,size=7,main=plot_title,
           xlim=xlim,ylim=ylim,zlim=zlim)
    if(show_lines==TRUE)
    {
      points3d(x=c0[1],y=c0[2],z=c0[3],size=10)#plot the center of all half-lines in black
      K=length(unique(dclass$class))
      for(k in 1:K)#print all half-lines with length from c0 to the center of each half-line cluster
      {
        clust1=colMeans(input_dat[dclass$class==k,,drop=FALSE])
        arrow3d(p0=as.numeric(c0),p1=clust1,barblen=0,color=k+1)#as numeric c0 because if data frame it crashes
      }
    }
    return()
  }

  update_Vmat=function(input_dat,Kappa,dclass,c0)#updates all direction vectors
  {
    #   update Vmat, the directions of the halflines
    #
    #   inputs:
    #     -input_dat: data frame of input data with rows=samles and cols=dimensions
    #     -Kappa: number of clusters (halflines)
    #     -dclass: data frame of class label for each sample
    #     -c0: current value for the center of all half-lines
    #
    #   return values:
    #     -Vmat: matrix whose K rows are the new direction vectors

    N=nrow(input_dat)#number of samples
    P=ncol(input_dat)
    Vmat=matrix(data = 0,nrow = Kappa,ncol = P)
    for(k in 1:Kappa)#update each direction vector
    {
      Vmat[k,]=update_v(input_dat=input_dat,dclass=dclass,cluster_selected=k,c0=c0)
    }
    return(Vmat)
  }

  update_v=function(input_dat,dclass,cluster_selected,c0)
  {
    #   update the direction vector v as the average direction of the selected cluster
    #
    #   inputs:
    #     -input_dat: data frame of input data with rows=samles and cols=dimensions
    #     -dclass: data frame of class label for each sample
    #     -cluster_selected: numeric, corresponding to the cluster label for which the direction is updated
    #     -c0: current value for the center of all half-lines
    #
    #   return values:
    #     -v_cm: updated direction vector

    data_centered=input_dat[dclass$class==cluster_selected,,drop=FALSE]
    data_centered=t(t(data_centered)-c0)
    v_cm=as.numeric(colMeans(data_centered))
    #v_cm=v_cm/sqrt(v_cm%*%v_cm)#normalize to unit length
    v_cm=v_cm/sqrt(sum(v_cm^2))#normalize to unit length
    return(v_cm)
  }

  update_c=function(input_dat,Kappa,dclass,c_old,Vmat)#updates the center/beginning of all halflines c, given data dat and direction vectors v1,v2,v3
  {
    #   update c, the center of the halflines
    #
    #   inputs:
    #     -input_dat: data frame of input data with rows=samles and cols=dimensions
    #     -Kappa: number of clusters (halflines)
    #     -dclass: data frame of class label for each sample
    #     -c_old: current (old) value for the center of all half-lines
    #     -Vmat: matrix whose K rows are the current direction vectors
    #
    #   return values:
    #     -c_new: new value for the center of all half-lines

    N=nrow(input_dat)#number of samples
    P=ncol(input_dat)#number of dimensions

    r_nk=matrix(0,nrow = N,ncol = Kappa)
    for(k in 1:Kappa)
    {
      r_nk[dclass$class==k,k]=1
    }
    Id=diag(P)

    V_lst=list()
    A_lst=list()
    for(k in 1:Kappa)
    {
      #V_lst[[k]]=Id-(Vmat[k,]%o%Vmat[k,]/as.numeric(Vmat[k,]%*%Vmat[k,]))
      V_lst[[k]]=Id-(Vmat[k,]%o%Vmat[k,])#no normalization since v_k has unit length
      A_lst[[k]]=t(V_lst[[k]])%*%V_lst[[k]];
    }

    #calculate the new center Co
    left=0;
    N_k=rep(0,Kappa)#Nk: number of elemets in cluster k with negative dot product (angle>90)
    M_k=rep(0,Kappa)#Mk: number of elemets in cluster k with positive dot product (angle<=90)

    for(i in 1:N)
    {
      x=as.numeric(input_dat[i,])#otherwise matrix operations won't work
      x_hat=x-c_old#take the C1 vector as the starting point of the axes of the vector space
      for(k in 1:Kappa)
      {
        if(r_nk[i,k]==1)#the data point belongs to cluster 1
        {
          # left=left+x%*%A1
          v=Vmat[k,]
          if(sign(as.numeric(x_hat%*%v))>=0)#check dot product positive -> must be close to half-line, not line
          {
            left=left+x%*%A_lst[[k]]
            M_k[k]=M_k[k]+1
          }else#point on the wrong side, measure distance from the center
          {
            left=left+x
            N_k[k]=N_k[k]+1
          }
        }
      }
    }

    # right=lambda*(N1+N2+N3)*Id+M1*A1+M2*A2+M3*A3
    N_all=sum(N_k)
    M_times_A=matrix(0,nrow = P, ncol = P)

    for(k in 1:Kappa)
    {
      M_times_A=M_times_A+M_k[k]*A_lst[[k]]
    }
    right=N_all*Id+M_times_A
    # right=solve(right)
    c_new=left%*%solve(right)
    return(as.numeric(c_new))
  }

  total_clust_error=function(input_dat,dclass,c0,Vmat)#computes the overall error of the clustering of all data points
  {
    #   calculate the overall clustering cost
    #
    #   inputs:
    #     -input_dat: data frame of input data with rows=samles and cols=dimensions
    #     -Kappa: number of clusters (halflines)
    #     -dclass: data frame of class label for each sample
    #     -c0: current value for the center of all half-lines
    #     -Vmat: matrix whose K rows are the current direction vectors
    #
    #   return values:
    #     -err: total clustering cost

    data=input_dat
    c0=as.numeric(c0)
    err=0
    for(i in 1:nrow(data))#compute the error (distance from halfline object) for every data point
    {
      pos=dclass$class[i]#select the direction vector according to the class of the data point
      v=Vmat[pos,]
      v_outer=v%o%v
      x=as.numeric(data[i,])#otherwise matrix operations won't work
      # print('x');print(x);print('c0');print(c0);print(x-c0)
      x_hat=x-c0#take the C1 vector as the starting point of the axes of the vector space
      # print('x-c0 happened')
      if(sign(as.numeric(x_hat%*%v))>=0)#check dot product positive -> must be close to half-line, not line
      {
        err=err+line_dist(x,c0,v_outer=v_outer)
        # pdprod[i]=1#the data point has a positive dot product with the half line
      }else#point on the wrong side, measure distance from the center
      {
        err=err+point_dist(x,c0)
      }
    }

    return(err)
  }

  kbranch.global_internal=function(input_dat,Kappa,c0=NULL,Vmat=NULL,silent=FALSE,keep_c_fixed=FALSE)
  {
    #   internal function for clustering k halflines
    #   called by the wrapper function kbranch.global()
    #
    #   inputs:
    #     -input_dat: data frame of input data with rows=samles and cols=dimensions
    #     -Kappa: number of clusters (halflines)
    #     -dclass: data frame of class label for each sample
    #     -c0: initial value for the center of all half-lines
    #     -Vmat: matrix whose K rows are the initial direction vectors
    #     -silent: if TRUE, generate output to show progress of the method as it runs (for debugging)
    #     -keep_c_fixed: if TRUE, keep the center of the halflines fixed
    #
    #   return values:
    #     -dclass: cluster assignment for each sample
    #     -Kappa: number of clusters (halflines)
    #     -err: total clustering cost
    #     -iters: total iterations of the algorithm
    #     -c0: final value for the center of all half-lines
    #     -Vmat: matrix whose K rows are the final direction vectors
    #     -clust_counts: number (count) of samples in each of the clusters

    N=nrow(input_dat)#number of samples
    P=ncol(input_dat)#the number of dimensions, the data are in R^d (d-dimensional)

    #normalize each direction vector to unit length
    Vmat=t(apply(Vmat,1,function(x){x/sqrt(sum(x^2))}))

    # print('assigning...')
    # print(Vmat)
    # dclass=assign_to_closest_hline2(xdata=input_dat,c1,c2,c3,v1,v2,v3) #initialize class labels
    dclass=assign_to_closest_hline(Kappa=Kappa,input_dat=input_dat,c0=c0,Vmat=Vmat) #initialize class labels

    # print('classes assigned')

    #perfrom check that all clusters have samples
    if(length(unique(dclass$class))<Kappa)
    {
      if(silent==FALSE)
      {
        print(paste('NOT all clusters have been assigned samples, only',length(unique(dclass$class)),'of',Kappa))
        print(unique(dclass$class))
        print('Qutting the clustering process...')
      }
      #crash detected, return NULL error and zero iterations
      retval=list(dclass=dclass,Kappa=Kappa,err=NULL,iters=0,c0=c0,Vmat=Vmat)
      return(retval)
    }else
    {
      if(silent==FALSE){print(paste('OK, all clusters have been assigned samples:',length(unique(dclass$class)),'of',Kappa))}
    }

    if((P==2)&&(silent==FALSE))
    {
      plot_hline_clust2D(input_dat,c0=c0,Vmat=Vmat,dclass=dclass,plot_title=paste('Initialization'))
    }

    c0_start=c0
    dclass_old=rep(0,length(dclass$class))
    i=1
    while(sum(dclass_old-dclass$class)!=0)#stop when no samples change clusters
    {
      if(silent==FALSE){print(paste('round',i))}

      if(keep_c_fixed==FALSE)#update the center, unless it is fixed
      {
        # print('updating c...')
        #update the center of the half-lines
        c0=update_c(Kappa=Kappa,input_dat = input_dat,dclass = dclass,Vmat,c_old = c0)#;print(c0)
      }
      #update the direction vectors of the half-lines
      Vmat_old=Vmat
      Vmat=update_Vmat(Kappa=Kappa,input_dat=input_dat,dclass=dclass,c0=c0)
      if(anyNA(Vmat))#NA rows indicate directions with no samples assigned to them after updating the cluster labels in the previous round
      {
        #print('@@@@@@@@@@@@@@');print(unique(dclass$class));plot_hline_clust2D(input_dat,c0=c0,Vmat=Vmat_old,dclass=dclass,plot_title='Problem')
        for(rcount in 1:nrow(Vmat))
        {
          if(anyNA(Vmat[rcount,]))#keep the direction from the old Vmat
          {
            Vmat[rcount,]=Vmat_old[rcount,]
          }
        }
      }
      # print('Vmat ok')

      #re-assign samples to clusters
      dclass_old=dclass$class#save old values to check for convergence
      dclass=assign_to_closest_hline(Kappa=Kappa,input_dat=input_dat,c0=c0,Vmat=Vmat)
      #perfrom check that all clusters have sumples
      if(length(unique(dclass$class))<Kappa)
      {
        if(silent==FALSE)
        {
          print(paste('NOT all clusters have been assigned samples, only',length(unique(dclass$class)),'of',Kappa))
          print('Quitting the clustering process...')
        }
        #crash detected, return NULL error and zero iterations
        retval=list(dclass=dclass,Kappa=Kappa,err=NULL,iters=0,c0=c0,Vmat=Vmat)
        return(retval)
      }else
      {
        if(silent==FALSE){print(paste('OK, all clusters have been assigned samples:',length(unique(dclass$class)),'of',Kappa))}
      }

      if((P==2)&&(silent==FALSE))
      {
        #plot_hline_clust2D(input_dat,c0=c0,Vmat=Vmat,dclass=dclass,plot_title=paste('after iteration',i,'| with starting c0:',round(c0_start[1],2),',',round(c0_start[2],2)))
        plot_hline_clust2D(input_dat,c0=c0,Vmat=Vmat,dclass=dclass,plot_title=paste('after iteration',i))
      }
      i=i+1
    }
    i=i-1;#correct the value of i from the last iteration
    err=total_clust_error(input_dat,dclass,c0,Vmat)
    if(silent==FALSE)
    {
      print(paste('center',c0))
      print(paste('error',err))
      if(P==2)
      {
        #plot_hline_clust2D(input_dat,c0=c0,Vmat=Vmat,dclass=dclass,plot_title=paste('Iterations',i,'| starting c0:(',round(c0_start[1],2),',',round(c0_start[2],2),') | error:',round(err,2)))
        plot_hline_clust2D(input_dat,c0=c0,Vmat=Vmat,dclass=dclass,plot_title=paste('Iterations',i))
      }
    }
    #compute the number of elements in each cluster
    clust_counts=rep(-1,Kappa)
    for(k in 1:Kappa){clust_counts[k]=sum(dclass$class==k)}
    #return a list containing all classification data
    #retval=list(dclass=dclass,Kappa=Kappa,lambda=lambda,err=err,iters=i,c0=c0,Vmat=Vmat,clust_counts=clust_counts)
    retval=list(dclass=dclass,Kappa=Kappa,err=err,iters=i,c0=c0,Vmat=Vmat,clust_counts=clust_counts)
    return(retval)
  }

  kbranch.global_internal_PAM=function(input_dat,Dmat=NULL,Kappa,c0_ix=NULL,Vmat_ix=NULL,silent=FALSE,keep_c_fixed=FALSE)
  {
    #   internal function for clustering k halflines
    #   called by the wrapper function kbranch.global()
    #   uses a PAM like process to identify the K halflines
    #
    #   inputs:
    #     -input_dat: data frame of input data with rows=samles and cols=dimensions
    #     -Dmat: matrix containing sample distances
    #     -Kappa: number of clusters (halflines)
    #     -c0_ix: position (row index in input_dat) of the center sample
    #     -Vmat_ix: positions (row indices in input_dat) of the direction samples
    #     -silent: if TRUE, generate output to show progress of the method as it runs (for debugging)
    #     -keep_c_fixed: if TRUE, keep the center of the halflines fixed
    #
    #   return values:
    #     -dclass: cluster assignment for each sample
    #     -Kappa: number of clusters (halflines)
    #     -err: total clustering cost
    #     -iters: total iterations of the algorithm
    #     -c0_ix: position (row index in input_dat) of the center sample
    #     -Vmat_ix: positions (row indices in input_dat) of the direction samples
    #     -clust_counts: number (count) of samples in each of the clusters

    N=nrow(input_dat)#number of samples
    P=ncol(input_dat)#the number of dimensions, the data are in R^d (d-dimensional)
    if(is.null(Dmat)==TRUE)
    {
      if(silent==FALSE){'the distance matrix for all samples wasn\'t provided, calculating...'}
      Dmat=compute_all_distances(input_dat)
      if(silent==FALSE){'...all distances were calculated.'}
    }

    medoid_ix=c(c0_ix,Vmat_ix)#the indexes (rows) of all medoid samples in the data set
    assgn=assign_to_closest_hline_PAM(Dmat=Dmat,Kappa=Kappa,input_dat=input_dat,c0_ix=c0_ix,Vmat_ix=Vmat_ix)#assign data points to the closest half line
    err_old=assgn$total_error
    dclass=assgn$dclass#just in case the initialization is the best case

    if(silent==FALSE)
    {
      Vmat_temp=input_dat[Vmat_ix,]
      Vmmat=t(t(as.matrix(Vmat_temp))-as.numeric(input_dat[c0_ix,]))
      if(P==2)
      {
        plot_hline_clust2D(input_dat=input_dat,c0=input_dat[c0_ix,],Vmat=Vmmat,dclass=dclass,plot_title=paste('initialization'),len=10^10,show_lines=TRUE)
        points(input_dat[Vmat_ix,],pch=4,cex=3)
      }else
      {
        plot_hline_clust3D(input_dat=input_dat,c0=input_dat[c0_ix,],Vmat=Vmmat,dclass=dclass,plot_title=paste('initialization'),len=10^10,show_lines=TRUE)
        points3d(input_dat[Vmat_ix,],pch=4,cex=3)
      }
    }

    iter=1
    error_decreases=TRUE

    Vmat_ix_new=Vmat_ix #initialize the direction vectors
    while(error_decreases==TRUE)
    {
      if(silent==FALSE){print(paste('iteration',iter))}

      err_new=err_old #initialize the error or the new iteration


      if(keep_c_fixed==FALSE)
      {
        if(silent==FALSE){print('exploring the center')}
        #search all neighboring states for a half-line center
        for(i in 1:N)#search for a new center in the samples
        {
          if(is.element(i,medoid_ix)==FALSE)#only check data points that are not medoids already
          {
            assgn=assign_to_closest_hline_PAM(Dmat=Dmat,Kappa=Kappa,input_dat=input_dat,c0_ix=i,Vmat_ix=Vmat_ix)#assign data points to the closest half line
            err_cur=assgn$total_error
            if(err_cur<err_new)#new center has better error, select it
            {
              err_new=err_cur #save the new error (currently best)
              c0_ix_new=i #replace the direction vector
              assgn_new=assgn #replace the cluster assignments
              change_center=TRUE #flag: change the center
            }
          }
        }
      }else{if(silent==FALSE){print('the center stays fixed...')}}


      if(silent==FALSE){print('exploring the direction vectors')}
      #search all neighboring states for direction vectors
      for(k in 1:Kappa)
      {
        #search for a new direction vector in the samples
        for(i in 1:N)
        {
          if(is.element(i,medoid_ix)==FALSE)#only check data points that are not medoids already
          {
            Vmat_ix_tmp=Vmat_ix_new
            Vmat_ix_tmp[k]=i
            assgn=assign_to_closest_hline_PAM(Dmat=Dmat,Kappa=Kappa,input_dat=input_dat,c0_ix=c0_ix,Vmat_ix=Vmat_ix_tmp)#assign data points to the closest half line
            err_cur=assgn$total_error
            if(err_cur<err_new)#new direction has better error, select it
            {
              err_new=err_cur #save the new error (currently best)
              Vmat_ix_new[k]=i #replace the direction vector
              assgn_new=assgn #replace the cluster assignments
              change_center=FALSE #flag: change the direction vectors, NOT the center
            }
          }
        }
      }

      if(err_new>=err_old)#reached local minimum, exit the loop
      {
        if(silent==FALSE){print('reached local minimum, exiting...')}
        #do not update Vmat_ix, the value of the previous iteration is better
        error_decreases=FALSE
      }else#the error was decreased, save/update the improved parameters
      {
        err_old=err_new #update the old value for the error
        if(change_center==TRUE)
        {
          c0_ix=c0_ix_new #save the improved center
          changed_var='center'
        }else
        {
          Vmat_ix=Vmat_ix_new #save the improved direction vectors
          changed_var='directions'
        }
        medoid_ix=c(c0_ix,Vmat_ix)
        # assignment=assgn#assignment of samples to clusters + clustering error
        dclass=assgn_new$dclass#save the 'best' cluster assignments so far
        if(silent==FALSE)
        {
          print('error decreased, continue...')

          Vmat_temp=input_dat[Vmat_ix,]
          Vmmat=t(t(as.matrix(Vmat_temp))-as.numeric(input_dat[c0_ix,]))
          if(P==2)
          {
            plot_hline_clust2D(input_dat=input_dat,c0=input_dat[c0_ix,],Vmat=Vmmat,dclass=dclass,plot_title=paste('iteration',iter,'changed',changed_var),len=10^10,show_lines=TRUE)
            points(input_dat[Vmat_ix,],pch=4,cex=3)
          }else
          {
            print('plotting 3D')
            plot_hline_clust3D(input_dat=input_dat,c0=input_dat[c0_ix,],Vmat=Vmmat,dclass=dclass,plot_title=paste('iteration',iter,'changed',changed_var),len=10^10,show_lines=TRUE)
            points3d(input_dat[Vmat_ix,],pch=4,cex=3)
          }

        }
        iter=iter+1
      }
    }
    clust_counts=rep(-1,Kappa)
    for(k in 1:Kappa){clust_counts[k]=sum(dclass$class==k)}
    #return a list containing all classification data
    # retval=list(dclass=dclass,Kappa=Kappa,lambda=lambda,err=err_old,iters=iter,c0_ix=c0_ix,Vmat_ix=Vmat_ix,clust_counts=clust_counts)
    retval=list(dclass=dclass,Kappa=Kappa,err=err_old,iters=iter,c0_ix=c0_ix,Vmat_ix=Vmat_ix,clust_counts=clust_counts)
    return(retval)

  }

  assign_to_closest_hline_PAM=function(Dmat,Kappa,input_dat,c0_ix,Vmat_ix)#assign data points to the closest half line
  {
    #   Assingn each sample to it's closest halfline in the case of PAMlike k halflines clustering
    #
    #   inputs:
    #     -Dmat: matrix containing sample distances
    #     -Kappa: number of clusters (halflines)
    #     -input_dat: data frame of input data with rows=samles and cols=dimensions
    #     -c0_ix: position (row index in input_dat) of the center sample
    #     -Vmat_ix: positions (row indices in input_dat) of the direction samples
    #
    #   return values:
    #     -dclass: cluster assignment for each sample
    #     -total_error: total clustering cost

    if(missing(input_dat))
    {
      print('missing data....')
    }
    # print(paste('input_dat dim',dim(input_dat)))
    c0=as.numeric(input_dat[c0_ix,])
    Vmat=as.matrix(input_dat[Vmat_ix,])

    total_error=0#the overall clustering error
    xdata_class=data.frame(class=rep(0,nrow(input_dat)))
    # xdata_class=data.frame(class=rep(0,nrow(input_dat)),wrong_x=rep(NA,nrow(input_dat)),wrong_y=rep(NA,nrow(input_dat)))
    for(i in 1:nrow(input_dat))
    {
      # print(paste('i is',i))
      x=as.numeric(input_dat[i,])#otherwise matrix operations won't work
      nc=-1;#which is the nearest class
      nc_dist=Inf
      nc_error=-1
      x_hat=x-c0#take the C1 vector as the starting point of the axes of the vector space
      for(k in 1:Kappa)
      {
        v_k=Vmat[k,]-c0#current direction vector

        if(anyNA(Vmat)){
          print(paste('x_hat',as.numeric(x_hat)))
          print(paste('v_k',as.numeric(v_k)))
          print(Vmat)}

        if(sign(as.numeric(x_hat%*%v_k))>=0)#check dot product positive -> must be close to half-line, not line
        {
          d_ic_squared=Dmat[c0_ix,i];#print(paste('d_ic_squared',d_ic_squared))
          d_iv_squared=Dmat[Vmat_ix[k],i];#print(paste('d_iv_squared',d_iv_squared))
          d_cv_squared=Dmat[c0_ix,Vmat_ix[k]];#print(paste('d_cv_squared',d_cv_squared))

          nc_dist_new=abs(d_ic_squared- ((d_ic_squared+d_cv_squared-d_iv_squared)/(2*sqrt(d_cv_squared)))^2)#abs for numerical stabilitytability

          if(nc_dist_new<nc_dist)
          {
            nc=k;
            nc_dist=nc_dist_new
          }
        }else#point on the wrong side, measure distance from the center
        {
          nc_dist_new=Dmat[c0_ix,i]#Dmat contains the squared Euclidean distance
          if(nc_dist_new<nc_dist)
          {
            nc=k;
            nc_dist=nc_dist_new
          }
        }
      }

      xdata_class$class[i]=nc;#assign closest cluster
      total_error=total_error+nc_dist#add the distance to the clostest half-line to the overall error

      # View(xdata_class)
      if(nc==-1)#none of the three clusters was selected, all had negative dot-products
      {
        print(paste('------- ERROR! vector with x=',x[1],' y=',x[2],'-------'))
      }
    }

    return(list(dclass=xdata_class,total_error=total_error)) #return the class for each sample
  }

  #######################################################
  # end of internal functions of kbranch.global
  #######################################################

  N=nrow(input_dat)#the number of samples
  P=ncol(input_dat)#the number of dimensions, the data are in R^d (d-dimensional)
  #if no c0 or Vmat have been provided, initialize using K-means
  if( ( (is.null(c0))||(is.null(Vmat)) ))#&(medoids==FALSE) )#if the medoid version is run, no need to initialize manually
  {
    # print('got in the init phase')
    if(medoids==FALSE)#perform initialization of the usual K-halflines method
    {
      if(init_Kmeans==FALSE)
      {
        if(silent==FALSE){print('Initializing using random directions')}
      }else #initialize using the directions of K-means
      {
        if(silent==FALSE){print('Initializing the directions using K-means...')}
        # P=ncol(input_dat)#the number of dimensions, the data are in R^d (d-dimensional)
        Vmat=matrix(data = 0,nrow = Kappa,ncol = P)
        Cmat=matrix(data = 0,nrow = Kappa,ncol = P)
        km.out=kmeans(x=input_dat, centers=Kappa, nstart = nstart_kmeans)

        for(k in 1:Kappa)
        {
          kdata=input_dat[km.out$cluster==k,,drop=FALSE]#select the k-th cluster of kmeans
          # print(colMeans(kdata))
          Cmat[k,]=colMeans(kdata)
          # points(Cmat[k,])
        }
        # print(Cmat)
        c0_kmeans=colMeans(Cmat)
        dclass=data.frame(class=km.out$cluster)
        if(silent==FALSE){print('K-means run OK...')}

        # if(Kappa==2)
        if(Kappa<3)
        {
          #Kappa == 2 is a special case
          #take 2 random points to initialize the 2 direction vectors
          #because the usual initialization through K-Means returns 2
          #directions on the SAME line and the clustering fails
          for(k in 1:Kappa) #in 1:2
          {
            Vmat[k,]=runif(P)
          }
          Vmat_kmeans=Vmat

          #       Vmat_kmeans=update_Vmat(Kappa,input_dat,dclass=dclass,c0=Cmat[1,])
          #       Vmat_kmeans=Vmat_kmeans[1:2,]
          #       # print(Vmat_kmeans)

        }else #Kappa>=3, set the directions to the K-means directions centered at the K-Means center
        {
          Vmat_kmeans=update_Vmat(Kappa=Kappa,input_dat=input_dat,dclass=dclass,c0=c0_kmeans)
        }

        if(silent==FALSE)#print the k-Means initialization, if the results are 2 or 3 dimensional
        {
          # if(P==2){plot_hline_clust2D(input_dat,c0=c0_kmeans,Vmat=Vmat_kmeans,dclass=dclass,plot_title=paste('K-Means'))}
          if(P==2){plot(x=input_dat[,1],y=input_dat[,2],col=km.out$cluster+1,main='K-Means',xlab='x',ylab='y');points(x=c0_kmeans[1],y=c0_kmeans[2],pch=3, cex=2, col='blue3')}
          if(P==3){open3d();plot3d(x=input_dat[,1],y=input_dat[,2],z=input_dat[,3],xlab='x axis',ylab='y axis',zlab='z axis',col=dclass$class+1,size=7,main='K-Means')}
        }
      }
    }
    all_clusterings=vector("list", nstart)#pre-allocate memory for the list of all classification results
    all_clustering_errors=rep(0,nstart)#the clustering error of every successful clustering performed
    reps=1 #count the number of successful clusterings performed
    iterations_performed=1 #count the number of total iterations, to account for unsucessful clusterings, as well.
    while(reps<=nstart)
    {
      if(silent==FALSE){print(paste('######### External Iteration #########',iterations_performed))}

      if(is.null(fixed_center)==TRUE)#run the usual K-Halflines, without a fixed center
      {
        if(medoids==FALSE)#run the usual clustering algorithm (no the medoids version)
        {
          if(init_Kmeans==FALSE)#perform random initialization
          {
            #initialize c0 to a random sample
            #initialize v1,...,vk to random samples (centered at c0): e.g.vk=ck-c0

            # print('running as usual')
            random_center_pos=sample(x=NROW(input_dat),size=1)
            c0=as.numeric(input_dat[random_center_pos,])#select a random data point as the center of the half-lines

            input_dat2=input_dat[-random_center_pos,]
            Vmat=as.matrix(input_dat2[sample(nrow(input_dat2),Kappa),])
            #         print('Vmat before')
            #         print(Vmat)
            for(i in 1:nrow(Vmat))
            {
              Vmat[i,]=Vmat[i,]-c0
            }
            Vmat_kmeans=Vmat

            retval=kbranch.global_internal(Kappa=Kappa,input_dat=input_dat,c0=c0,Vmat=Vmat_kmeans,silent=silent_internal)
          }else#initialize using k-means
          {
            #initialize c0 to a random sample
            #initialize v1,...,vk to the k-means directions

            # print('running as usual')
            c0=as.numeric(input_dat[sample(x=NROW(input_dat),size=1),])#select a random data point as the center of the half-lines
            # Vmat_kmeans=update_Vmat(Kappa,input_dat,dclass=dclass,c0=c0)
            retval=kbranch.global_internal(Kappa=Kappa,input_dat=input_dat,c0=c0,Vmat=Vmat_kmeans,silent=silent_internal)
          }
        }else #run the medoids version (with random initialization)
        {
          init_medoids_ix=sample(1:nrow(input_dat),Kappa+1)
          c0=init_medoids_ix[1]#the first sample is the center of the halflines
          Vmat=init_medoids_ix[-1]#the rest Kappa samples are the direction vectors
          # print('starting internal')
          retval=kbranch.global_internal_PAM(Dmat=Dmat,Kappa=Kappa,input_dat=input_dat,c0_ix=c0,Vmat_ix=Vmat,silent=silent_internal)
        }
      }else #run K-Halflines with a fixed center and only update the direction vectors
      {
        if(medoids==FALSE)
        {
          if(NROW(fixed_center)==1)#it is only one value: refers to the row of the sample to use as a fixed center
          {
            # print('NROW 1')
            # if(Kappa==1)
            c0=as.numeric(input_dat[fixed_center,])

            # c0=input_dat[fixed_center,]
            input_dat2=input_dat[-fixed_center,]
            # View(input_dat2)
            Vmat=as.matrix(input_dat2[sample(nrow(input_dat2),Kappa),])
            # View(Vmat)
            # if((Kappa==1)&(NCOL(input_dat)<3)){Vmat=t(Vmat)}
            # print(paste('NCOL(Vmat)',NCOL(Vmat)))
            # print(paste('NROW(c0)',NROW(c0)))
            if(NCOL(Vmat)!=NROW(c0)){Vmat=t(Vmat)}#need to transpose Vmat to get it right...
          }else#use the vector fixed center as the fixed center of the halflines
          {
            # print('NROW NOT 1')

            c0=as.numeric(fixed_center)
            Vmat=as.matrix(input_dat[sample(nrow(input_dat),Kappa),])
          }
          # print(paste('c0 is',c0))
          # print(NROW(c0))
          # print(NCOL(c0))
          # print('and vmat...')
          # print(Vmat)
          #print(dim(Vmat))
          for(i in 1:nrow(Vmat))
          {
            Vmat[i,]=Vmat[i,]-c0
          }
          Vmat_kmeans=Vmat
          retval=kbranch.global_internal(Kappa=Kappa,input_dat=input_dat,c0=c0,Vmat=Vmat_kmeans,silent=silent_internal,keep_c_fixed=T)
        }else#run the medoids version of the algorithm
        {
          c0=fixed_center
          set_ix=1:N
          set_ix=set_ix[-which(set_ix==c0)]#remove the index of the fixed center from the set of all available indices
          Vmat=sample(set_ix,Kappa)
          retval=kbranch.global_internal_PAM(Dmat=Dmat,Kappa=Kappa,input_dat=input_dat,c0_ix=c0,Vmat_ix=Vmat,silent=silent_internal,keep_c_fixed=T)
        }
      }
      # print('got out of the run medoids')
      if(is.null(retval$err)==FALSE)
      {
        all_clusterings[[reps]]=retval
        all_clustering_errors[reps]=retval$err
        reps=reps+1
      }
      iterations_performed=iterations_performed+1
    }
    if(silent==FALSE){print('clustering errors:');print(all_clustering_errors)}
    if(silent==FALSE){print('min_error:');print(min(all_clustering_errors))}
    pos_min=which(all_clustering_errors==min(all_clustering_errors))[1]
    retval=all_clusterings[[pos_min]]#keep only the best clustering
    retval$all_clustering_errors=all_clustering_errors#return the errors for all clusterings performed
    retval$all_clusterings=all_clusterings#return data for all clusterings performed

    #calculate the GAP statistic
    resGAP=calculate_GAP(input_dat=input_dat,clustering_error=retval$err,Kappa=Kappa,cluster_labels=retval$dclass$class,
                                  nstart_GAP=nstart_GAP,B_GAP=B_GAP,fixed_center=fixed_center,medoids=medoids,Dmat=Dmat,init_Kmeans = init_Kmeans)
    retval$GAP=resGAP$GAP
    retval$GAPl=resGAP$GAPl
    retval$GAP_orig=resGAP$GAP_orig
    retval$GAP_orig_no_log=resGAP$GAP_orig_no_log

    retval$GAP.sd=resGAP$GAP.sd
    retval$GAPl.sd=resGAP$GAPl.sd
    retval$GAP_orig.sd=resGAP$GAP_orig.sd
    retval$GAP_orig_no_log.sd=resGAP$GAP_orig_no_log.sd

    if(show_plots==TRUE)#print the best clustering results, if the results are 2 dimensional
    {

      if(medoids==FALSE)
      {
        if(P==2){plot_hline_clust2D(input_dat,c0=retval$c0,Vmat=retval$Vmat,dclass=retval$dclass,plot_title='K-Halflines')}
        if(P==3)
        {plot_hline_clust3D(input_dat,c0=retval$c0,Vmat=retval$Vmat,dclass=retval$dclass,plot_title='K-Halflines',show_lines=show_lines)}
      }else
      {
        Vmat_temp=input_dat[retval$Vmat,]
        Vmmat=t(t(as.matrix(Vmat_temp))-as.numeric(input_dat[retval$c0,]))
        if(P==2)
        {
          plot_hline_clust2D(input_dat,c0=input_dat[retval$c0,],Vmat=Vmmat,dclass=retval$dclass,plot_title='K-Halflines')
          if(ncol(as.matrix(input_dat[retval$c0,]))==2)
          {
            points(as.matrix(input_dat[retval$c0,]),pch=3,cex=3)#mark the center
          }else
          {
            points(t(as.matrix(input_dat[retval$c0,])),pch=3,cex=3)#mark the center
          }
          points(input_dat[retval$Vmat,],pch=4,cex=3)#mark the direction vectors
        }
        if(P==3)
        {
          plot_hline_clust3D(input_dat,c0=input_dat[retval$c0,],Vmat=Vmmat,dclass=retval$dclass,plot_title='K-Halflines',show_lines=show_lines)
        }
      }

      # rgl.set() #to select the proper rgl window
      # rgl.snapshot( filename = '3d_final.png', fmt = "png") #to print the window
      # rgl.postscript(filename = 'local_guo_unfiltered_scaled_largeK_100.pdf', fmt='pdf')
    }


  }else
  {
    if(silent==FALSE){print('c0 and Vmat manually provided')}
    #if c0 and Vmat have been provided, simply run the function
    retval=kbranch.global_internal(Kappa,input_dat,c0,Vmat,silent)
  }

  retval$cluster=retval$dclass$class
  retval$call=match.call()
  retval$dclass=NULL
  return(retval)

}

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

  # N=nrow(input_dat)#number of samples
  # P=ncol(input_dat)#the number of dimensions, the data are in R^d (d-dimensional)
  # # print(paste('N',N,'P',P))
  # Dist=matrix(data = NA,nrow = N,ncol = N)
  # for(i in 1:N)
  # {
  #   for(j in 1:N)
  #   {
  #     if(i==j)#distance of sample to itself is zero
  #     {
  #       Dist[i,j]=0
  #     }else
  #     {
  #       if(is.na(Dist[i,j])==TRUE)#has not been computed yet
  #       {
  #         x1=as.numeric(input_dat[i,])#the i-th sample
  #         x2=as.numeric(input_dat[j,])#the j-th sample
  #         d=norm(as.matrix(x1-x2),'2')^2#compute the squared Euclidian distance
  #         # d=abs(norm(as.matrix(x1-x2),'2'))#compute the Euclidian distance, abs for numerical stability
  #         Dist[i,j]=d;
  #         Dist[j,i]=d;#the distance is symmetric, so the matrix Dist is symmetric
  #       }
  #     }
  #
  #   }
  # }
  # return(Dist)
  return(as.matrix(stats::dist(input_dat)^2))#squared Euclidean distance matrix
}

#'Find K nearest neighbours
#'
#'Finds the K nearest neighbours of a sample, given a matrix containing sample distances.
#'
#'@param n: sample for which the nearest neighbors are computed
#'@param K: number of nearest neighbors to be considered
#'@param Dist: matrix containing sample distances, initially computed by compute_all_distances
#'
#'@return a list with elements:
#'\itemize
#'{
#'  \item - pos: the positions (rows) of the N nearest neighbours.
#'  \item - dst: the distances of the N nearest neighbours.
#'}
#'
#'@examples
#'\dontrun{
#'  Dist=compute_all_distances(scdata.3lines.simulated6genes_subsampled)
#'
#'  #find the 5 nearest neghbours of the second sample
#'  neibs=find_K_nearest(n=2,K=5,Dist=Dist)
#' }
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

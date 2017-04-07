#'
#'Clustering on K-Branches
#'
#'Clusters data on K-Branches (halflines) with a common center
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
#'\itemize{
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
#' #cluster the 2D data on three halflines
#' set.seed(1)
#'
#' #load the data
#' data(scdata.3lines.simulated6genes_subsampled)
#' raw_dat <- scdata.3lines.simulated6genes_subsampled
#'
#' #perform diffusion map dimensionality reduction
#' dmap <- destiny::DiffusionMap(raw_dat, sigma = 1000)
#'
#' #keep the first 2 diffusion components
#' input_dat <- destiny::as.data.frame(dmap)[, 1:2]
#'
#' #cluster with K=3
#' clust <- kbranches.global(input_dat, Kappa = 3)
#'
#' #plot the clustering results
#' plot(input_dat, pch=21, col=clust$cluster, bg=clust$cluster, main = 'K-Branch clustering')
#'
#'
#'@export
kbranches.global=function(input_dat,Kappa,Dmat=NULL,init_Kmeans=TRUE,c0=NULL,Vmat=NULL,nstart=20,nstart_GAP=20,nstart_kmeans=20,B_GAP=NULL,
                        fixed_center=NULL,medoids=FALSE,silent=TRUE,silent_internal=TRUE,show_plots=FALSE,show_lines = TRUE,show_plots_GAP=FALSE)
{
  #   external function for clustering k halflines
  #   serves as wrapper to kbranches.global_internal and also initializes c0 and Vmat if values have not been provided
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
  #     -cluster: cluster assignment for each sample (integer)
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
  # Internal functions of kbranches.global
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
    
    xdata_class=data.frame(class=rep(0,nrow(input_dat)),wrong_x=rep(NA,nrow(input_dat)),wrong_y=rep(NA,nrow(input_dat)),wrong_side=rep(NA,nrow(input_dat)))
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
            xdata_class$wrong_side[i]=F
            nc=k;
            nc_dist=nc_dist_new
          }
        }else#point on the wrong side, measure distance from the center
        {
          nc_dist_new=point_dist(x,c0)
          if(nc_dist_new<nc_dist)
          {
            xdata_class$wrong_side[i]=T
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

  update_Vmat=function(input_dat,Kappa,dclass,c0,Vmat_old=NULL)#updates all direction vectors
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
      Vmat[k,]=update_v(input_dat=input_dat,dclass=dclass,cluster_selected=k,c0=c0,v_old=Vmat_old[k,])
    }
    return(Vmat)
  }

  update_v=function(input_dat,dclass,cluster_selected,c0,v_old=NULL)
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

    # data_centered=input_dat[dclass$class==cluster_selected,,drop=FALSE]
    # data_centered=t(t(data_centered)-c0)
    # v_cm=as.numeric(colMeans(data_centered))
    # #v_cm=v_cm/sqrt(v_cm%*%v_cm)#normalize to unit length
    # v_cm=v_cm/sqrt(sum(v_cm^2))#normalize to unit length
    # return(v_cm)
    
    ix=(dclass$class==cluster_selected)&(dclass$wrong_side==F)
    # print(paste('cluster',cluster_selected,'right:',sum(ix),'of',sum(dclass$class==cluster_selected)))
    # print(sum(ix))
    # plot(input_dat,main='test')
    # points(input_dat[ix,],col=cluster_selected+1,pch=24,bg=cluster_selected+1)
    if(sum(ix)==0)
    {
      print('no points on the right side, using all data points...')
      ix=dclass$class==cluster_selected
    }
    data_centered=input_dat[ix,,drop=FALSE]
    
    # points(data_centered,pch=19,col=cluster_selected+1)
    
    # print(dim(data_centered))
    
    data_centered=t(t(data_centered)-c0)
    data_centered=as.matrix(data_centered)
    XtX=t(data_centered)%*%data_centered
    # XtX=data_centered%*%t(data_centered)
    # print(dim(XtX))
    
    # v_cm=as.numeric(colMeans(data_centered))
    v_cm=as.numeric(eigen(XtX)$vectors[,1])
    if(v_cm%*%v_old<0)
    {
      v_cm=-v_cm
    }
    
    # print(NROW(v_cm))
    # v_cm=v_cm/sqrt(v_cm%*%v_cm)#normalize to unit length
    v_cm=v_cm/sqrt(sum(v_cm^2))#normalize to unit length
    return(v_cm)
    
  }

  update_VmatInit=function(input_dat,Kappa,dclass,c0)#updates all direction vectors
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
      Vmat[k,]=update_vInit(input_dat=input_dat,dclass=dclass,cluster_selected=k,c0=c0)
    }
    return(Vmat)
  }
  
  update_vInit=function(input_dat,dclass,cluster_selected,c0)
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
      x_hat=x-c0#take the C1 vector as the starting point of the axes of the vector space
      if(sign(as.numeric(x_hat%*%v))>=0)#check dot product positive -> must be close to half-line, not line
      {
        err=err+line_dist(x,c0,v_outer=v_outer)
      }else#point on the wrong side, measure distance from the center
      {
        err=err+point_dist(x,c0)
      }
    }

    return(err)
  }

  kbranches.global_internal=function(input_dat,Kappa,c0=NULL,Vmat=NULL,silent=FALSE,keep_c_fixed=FALSE)
  {
    #   internal function for clustering k halflines
    #   called by the wrapper function kbranches.global()
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
    # errlist=list()
    #normalize each direction vector to unit length
    Vmat=t(apply(Vmat,1,function(x){x/sqrt(sum(x^2))}))

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
    c0_old=c0
    dclass_old=rep(0,length(dclass$class))
    i=1
    while(sum(dclass_old-dclass$class)!=0)#stop when no samples change clusters
    {
      if(silent==FALSE)
      {
        print(paste('round',i))
        if(P==2){plot_hline_clust2D(input_dat,c0=c0,Vmat=Vmat,dclass=dclass,plot_title=paste('iter',i,'start'))}
      }
      
      if(keep_c_fixed==FALSE)#update the center, unless it is fixed
      {
        c0_old=c0
        #update the center of the half-lines
        c0=update_c(Kappa=Kappa,input_dat = input_dat,dclass = dclass,Vmat,c_old = c0)#;print(c0)
        
        # plot_hline_clust2D(input_dat,c0=c0,Vmat=Vmat,dclass=dclass,plot_title=paste('iter',i,'start'))
        # points(c0[1],c0[2],pch=24,bg='black',col='black')
        
        # #re-assign samples to clusters
        # dclass_old=dclass$class#save old values to check for convergence
        # dclass=assign_to_closest_hline(Kappa=Kappa,input_dat=input_dat,c0=c0,Vmat=Vmat)
        # #perfrom check that all clusters have sumples
        # if(length(unique(dclass$class))<Kappa)
        # {
        #   if(silent==FALSE)
        #   {
        #     print(paste('NOT all clusters have been assigned samples, only',length(unique(dclass$class)),'of',Kappa))
        #     print('Quitting the clustering process...')
        #   }
        #   #crash detected, return NULL error and zero iterations
        #   retval=list(dclass=dclass,Kappa=Kappa,err=NULL,iters=0,c0=c0,Vmat=Vmat)
        #   return(retval)
        # }else
        # {
        #   if(silent==FALSE){print(paste('OK, all clusters have been assigned samples:',length(unique(dclass$class)),'of',Kappa))}
        # }
        
      }
      
      #update the direction vectors of the half-lines
      Vmat_old=Vmat
      Vmat=update_Vmat(Kappa=Kappa,input_dat=input_dat,dclass=dclass,c0=c0_old,Vmat_old=Vmat_old)
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
        points(Vmat[1,1],Vmat[1,2],pch=24,bg=2,col=2)
        points(Vmat[2,1],Vmat[2,2],pch=24,bg=3,col=3)
        points(Vmat[3,1],Vmat[3,2],pch=24,bg=4,col=4)
        
        points(c(0,Vmat[1,1]),c(0,Vmat[1,2]),type='l',col=2,lty='dashed')
        points(c(Vmat[1,1]),c(Vmat[1,2]),pch=24,bg=2,col=2)
        
        points(c(0,Vmat[2,1]),c(0,Vmat[2,2]),type='l',col=3,lty='dashed')
        points(c(Vmat[2,1]),c(Vmat[2,2]),pch=24,bg=3,col=3)
        
        points(c(0,Vmat[3,1]),c(0,Vmat[3,2]),type='l',col=4,lty='dashed')
        points(c(Vmat[3,1]),c(Vmat[3,2]),pch=24,bg=4,col=4)
      }
      # errlist[i]=total_clust_error(input_dat,dclass,c0,Vmat)
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

  kbranches.global_internal_PAM=function(input_dat,Dmat=NULL,Kappa,c0_ix=NULL,Vmat_ix=NULL,silent=FALSE,keep_c_fixed=FALSE)
  {
    #   internal function for clustering k halflines
    #   called by the wrapper function kbranches.global()
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
  # end of internal functions of kbranches.global
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

        }else #Kappa>=3, set the directions to the K-means directions centered at the K-Means center
        {
          Vmat_kmeans=update_VmatInit(Kappa=Kappa,input_dat=input_dat,dclass=dclass,c0=c0_kmeans)
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
            random_center_pos=sample(x=NROW(input_dat),size=1)
            c0=as.numeric(input_dat[random_center_pos,])#select a random data point as the center of the half-lines

            input_dat2=input_dat[-random_center_pos,]
            Vmat=as.matrix(input_dat2[sample(nrow(input_dat2),Kappa),])
            for(i in 1:nrow(Vmat))
            {
              Vmat[i,]=Vmat[i,]-c0
            }
            Vmat_kmeans=Vmat

            retval=kbranches.global_internal(Kappa=Kappa,input_dat=input_dat,c0=c0,Vmat=Vmat_kmeans,silent=silent_internal)
          }else#initialize using k-means
          {
            #initialize c0 to a random sample
            #initialize v1,...,vk to the k-means directions

            c0=as.numeric(input_dat[sample(x=NROW(input_dat),size=1),])#select a random data point as the center of the half-lines
            retval=kbranches.global_internal(Kappa=Kappa,input_dat=input_dat,c0=c0,Vmat=Vmat_kmeans,silent=silent_internal)
          }
        }else #run the medoids version (with random initialization)
        {
          init_medoids_ix=sample(1:nrow(input_dat),Kappa+1)
          c0=init_medoids_ix[1]#the first sample is the center of the halflines
          Vmat=init_medoids_ix[-1]#the rest Kappa samples are the direction vectors
          # print('starting internal')
          retval=kbranches.global_internal_PAM(Dmat=Dmat,Kappa=Kappa,input_dat=input_dat,c0_ix=c0,Vmat_ix=Vmat,silent=silent_internal)
        }
      }else #run K-Halflines with a fixed center and only update the direction vectors
      {
        if(medoids==FALSE)
        {
          if(NROW(fixed_center)==1)#it is only one value: refers to the row of the sample to use as a fixed center
          {
            c0=as.numeric(input_dat[fixed_center,])

            input_dat2=input_dat[-fixed_center,]
            Vmat=as.matrix(input_dat2[sample(nrow(input_dat2),Kappa),])
            if(NCOL(Vmat)!=NROW(c0)){Vmat=t(Vmat)}#need to transpose Vmat to get it right...
          }else#use the vector fixed center as the fixed center of the halflines
          {
            # print('NROW NOT 1')

            c0=as.numeric(fixed_center)
            Vmat=as.matrix(input_dat[sample(nrow(input_dat),Kappa),])
          }
          for(i in 1:nrow(Vmat))
          {
            Vmat[i,]=Vmat[i,]-c0
          }
          Vmat_kmeans=Vmat
          retval=kbranches.global_internal(Kappa=Kappa,input_dat=input_dat,c0=c0,Vmat=Vmat_kmeans,silent=silent_internal,keep_c_fixed=T)
        }else#run the medoids version of the algorithm
        {
          c0=fixed_center
          set_ix=1:N
          set_ix=set_ix[-which(set_ix==c0)]#remove the index of the fixed center from the set of all available indices
          Vmat=sample(set_ix,Kappa)
          retval=kbranches.global_internal_PAM(Dmat=Dmat,Kappa=Kappa,input_dat=input_dat,c0_ix=c0,Vmat_ix=Vmat,silent=silent_internal,keep_c_fixed=T)
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
    retval=kbranches.global_internal(Kappa,input_dat,c0,Vmat,silent)
  }

  retval$cluster=as.integer(retval$dclass$class)
  retval$call=match.call()
  retval$dclass=NULL
  return(retval)

}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kbranches-local.R
\name{kbranches.local}
\alias{kbranches.local}
\title{Local K-Branch clustering for identifying regions of interest!}
\usage{
kbranches.local(input_dat, Dmat = NULL, S_neib = NULL, S_quant = 0.1,
  S_GUI_helper = FALSE, parallel_ncores = NULL, min_radius_quantile = 0.5,
  logfile = "log.kbranches.local.txt", nstart = 5, nstart_GAP = 1,
  B_GAP = 5, medoids = FALSE, init_Kmeans = TRUE)
}
\arguments{
\item{input_dat:}{data frame of input data with rows=samles and cols=dimensions.}

\item{Dmat:}{matrix containing sample distances input_dat. Should be calculated using Dmat=compute_all_distances(input_dat)}

\item{S_neib:}{number of neighbours that defines the local neighbourhood. If set to NULL
a GUI helper will pop-up to assist in the selection of S_neib, providing a recommended value.}

\item{S_quant:}{the quantile of cumulative distance used to infer S_neib, the default value should usually suffice.}

\item{S_GUI_helper:}{If TRUE, a GUI helper will pop-up to aid in the selection of S_neib. If global_S is TRUE, then
the GUI will recommend a value for S_neib and visualize the neighbourhood for all data points. If global_S is false, the
GUI will only visualize the neighbourhood for every data point, since S_neib is selected automatically for each data point.}

\item{parallel_ncores:}{number of cores to use. If set to NULL (default) it will use the max number of
available cores. Defaults to NULL.}

\item{min_radius_quantile:}{the percentile to use for the identification of the 'median neighbourhood', defaults to 0.5 (median).}

\item{logfile:}{logfile to print output in the case of parallel computation.}

\item{nstart:}{number of initializations for clustering.
Defaults to 5, increase it (e.g. to 10) if the results are too noisy.}

\item{nstart_GAP:}{number of initializations for clustering when calculating the GAP statistic.
Defaults to 1, increase it (e.g. to 5 or 10) if the results are too noisy.}

\item{B_GAP:}{number of bootstrap datasets used to compute the GAP statistic.
Defaults to 5, increase it (e.g. to 10 or 100) if the results are too noisy.}

\item{medoids:}{if TRUE, the medoids version of kbranch will be used (slower).}

\item{init_Kmeans:}{if TRUE, use K-Means for the initialization of K-haflines, otherwise use random initialization}
}
\value{
a list with elements:
\itemize{
 \item - gap_scores: list of the four different gap scores for each sample.
 \item - call: the call of the function
 \item - S_neib: the global value of S_neib used, or 'local' if global_S was false
}
}
\description{
Perform local clustering using kbranches.global in order to generate a GAP score for each sample.
The GAP score of each sample can be subsequantly used to identify regions of interest
such as tips of branches or branching regions
}
\examples{
#this example might take some time to run
set.seed(1)

#load the data, already in diffusion map format
data(scdata.loop.dmap)
input_dat <- scdata.loop.dmap[, 1:2] #keep the first 2 diffusion components

#if the data are in not in diffusion space then
#performing diffusion map dimensionality reduction
#is necessary, for example:
#load(scdata.loop)
#dmap <- destiny::DiffusionMap(scdata.loop, sigma = 1000)
#input_dat <- destiny::as.data.frame(dmap)[, 1:2] #keep the first 2 diffusion components

#compute the distances among all samples
Dmat <- compute_all_distances(input_dat)

#perform local clustering to identify regions
#if you haven't specified the neighbourhood size S, it will be estimated
#set S_GUI_helper to FALSE for manual fine-tuning

res <- kbranches.local(input_dat = input_dat, Dmat = Dmat)

#identify regions of interest based on the GAP score
#of each sample computed by kbranches.local

#If smoothing_region and smoothing_region_thresh are NULL, a GUI will
#pop-up to aid in their selection. Press OK to update the results.
#When you are happy with the filtering press 'x' to close the window.

#smoothing: tip cell if at least 5 in 5 neighbors are tip cells
tips <- identify_regions(input_dat = input_dat, gap_scores = res$gap_scores, Dist = Dmat,
                         smoothing_region = 5, smoothing_region_thresh = 5, mode = 'tip')

#plot the separate tips
plot(input_dat, pch=21, col = tips$cluster + 1, bg = tips$cluster + 1, main='tip regions')

#smoothing: branching region cell if at least 10 in 10 neighbors are branching region cells
branch_reg <- identify_regions(input_dat = input_dat, gap_scores = res$gap_scores, Dist = Dmat,
                               smoothing_region = 10, smoothing_region_thresh = 10, mode='branch')

#plot the branching_region(s)
plot(input_dat, pch = 21, col = branch_reg$cluster + 1, bg = branch_reg$cluster + 1, main = 'branching regions')

###################################################
#          end of example                         #
###################################################
}


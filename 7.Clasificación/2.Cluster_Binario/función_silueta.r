fviz_nbclust <- function (x, FUNcluster = NULL, method = c("silhouette", "wss", "gap_stat"),
                          diss = NULL, k.max = 10, nboot = 100, verbose = interactive(),
                          barfill="steelblue", barcolor="steelblue", 
                          linecolor = "steelblue", print.summary = TRUE,  ...) 
{
  set.seed(123)
  if(k.max < 2) stop("k.max must bet > = 2")
  method = match.arg(method)
  if(!inherits(x, c("data.frame", "matrix")) & !("Best.nc" %in% names(x)))
    stop("x should be an object of class matrix/data.frame or ",
         "an object created by the function NbClust() [NbClust package].")
  
  # x is an object created by the function NbClust() [NbClust package]
  if(inherits(x, "list") & "Best.nc" %in% names(x)){
    best_nc <- x$Best.nc
    if(class(best_nc) == "numeric") print(best_nc)
    else if(class(best_nc) == "matrix") 
      .viz_NbClust(x, print.summary, barfill, barcolor)
  }
  else if(is.null(FUNcluster)) stop("The argument FUNcluster is required. ",
                                    "Possible values are kmeans, pam, hcut, clara, ...")
  else if(!is.function(FUNcluster)){
    stop(
      "The argument FUNcluster should be a function. ",
      "Check if you're not overriding the specified function name somewhere."
    )
  }
  else if(method %in% c("silhouette", "wss")) {
    
    if (is.data.frame(x)) x <- as.matrix(x)
    if(is.null(diss)) diss <- stats::dist(x)
    
    v <- rep(0, k.max)
    if(method == "silhouette"){
      for(i in 2:k.max){
        clust <- FUNcluster(x, i, ...)
        v[i] <- .get_ave_sil_width(diss, clust$cluster)
      }
    }
    else if(method == "wss"){
      for(i in 1:k.max){
        clust <- FUNcluster(x, i, ...)
        v[i] <- .get_withinSS(diss, clust$cluster)
      }
      
    }
    
    df <- data.frame(clusters = as.factor(1:k.max), y = v, stringsAsFactors = TRUE)
    
    ylab <- "Suma de cuadrados totales"
    if(method == "silhouette") ylab <- "Amplitud Promedio de Silueta"
    
    p <- ggpubr::ggline(df, x = "clusters", y = "y", group = 1,
                        color = linecolor, ylab = ylab,
                        xlab = "Número de k clúster",
                        main = "Número optimo de clúster"
    )
    if(method == "silhouette") 
      p <- p + geom_vline(xintercept = which.max(v), linetype=2, color = linecolor)
    
    return(p) 
  }
  
  else if(method == "gap_stat"){
    extra_args <- list(...)
    gap_stat <- cluster::clusGap(x, FUNcluster, K.max = k.max,  B = nboot, 
                                 verbose = verbose, ...)
    if(!is.null(extra_args$maxSE)) maxSE <- extra_args$maxSE
    else maxSE <- list(method = "firstSEmax", SE.factor = 1)
    p <- fviz_gap_stat(gap_stat,  linecolor = linecolor, maxSE = maxSE)
    return(p) 
  }
  
  
}

#' @rdname fviz_nbclust
#' @param gap_stat an object of class "clusGap" returned by the function
#'   clusGap() [in cluster package]
#' @param maxSE a list containing the parameters (method and SE.factor) for
#'   determining the location of the maximum of the gap statistic (Read the
#'   documentation ?cluster::maxSE). Allowed values for maxSE$method include: 
#'   \itemize{ \item "globalmax": simply corresponds to the global maximum,
#'   i.e., is which.max(gap) \item "firstmax": gives the location of the first
#'   local maximum \item "Tibs2001SEmax": uses the criterion, Tibshirani et al
#'   (2001) proposed: "the smallest k such that gap(k) >= gap(k+1) - s_{k+1}". 
#'   It's also possible to use "the smallest k such that gap(k) >= gap(k+1) -
#'   SE.factor*s_{k+1}" where SE.factor is a numeric value which can be 1
#'   (default), 2, 3, etc. \item "firstSEmax": location of the first f() value
#'   which is not larger than the first local maximum minus SE.factor * SE.f[],
#'   i.e, within an "f S.E." range of that maximum. \item
#'   see ?cluster::maxSE for more options }
#'   
#'   
#' @export
fviz_gap_stat <- function(gap_stat,  linecolor = "steelblue",
                          maxSE = list(method = "firstSEmax", SE.factor = 1)){
  if(!inherits(gap_stat, "clusGap"))
    stop("Only an object of class clusGap is allowed. (cluster package)")
  if(is.list(maxSE)){
    if(is.null(maxSE$method)) maxSE$method = "firstmax"
    if(is.null(maxSE$SE.factor)) maxSE$SE.factor = 1
  }
  else stop("The argument maxSE must be a list containing the parameters method and SE.factor")
  
  # first local max
  gap <- gap_stat$Tab[, "gap"]
  se <- gap_stat$Tab[, "SE.sim"]
  decr <- diff(gap) <= 0
  k <- .maxSE(gap, se, method = maxSE$method, SE.factor = maxSE$SE.factor)
  
  #k <- length(gap)
  #k = if (any(decr)) which.max(decr) else k
  
  df <- as.data.frame(gap_stat$Tab, stringsAsFactors = TRUE)
  df$clusters <- as.factor(1:nrow(df))
  df$ymin <- gap-se
  df$ymax <- gap + se
  p <- ggpubr::ggline(df, x = "clusters", y = "gap", group = 1, color = linecolor)+
    ggplot2::geom_errorbar(aes_string(ymin="ymin", ymax="ymax"), width=.2, color = linecolor)+
    geom_vline(xintercept = k, linetype=2, color = linecolor)+
    labs(y = "Estad?stico de Intervalo", x = "Número de k clúster",
         title = "Número optimo de clúster")
  p
}




# Get the Amplitud Promedio de Silueta
# ++++++++++++++++++++++++++
# Cluster package required
# d: dist object
# cluster: cluster number of observation
.get_ave_sil_width <- function(d, cluster){
  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("cluster package needed for this function to work. Please install it.")
  }
  ss <- cluster::silhouette(cluster, d)
  mean(ss[, 3])
}

# Get total within sum of square
# +++++++++++++++++++++++++++++
# d: dist object
# cluster: cluster number of observation
.get_withinSS <- function(d, cluster){
  d <- stats::as.dist(d)
  cn <- max(cluster)
  clusterf <- as.factor(cluster)
  clusterl <- levels(clusterf)
  cnn <- length(clusterl)
  
  if (cn != cnn) {
    warning("cluster renumbered because maximum != number of clusters")
    for (i in 1:cnn) cluster[clusterf == clusterl[i]] <- i
    cn <- cnn
  }
  cwn <- cn
  # Compute total within sum of square
  dmat <- as.matrix(d)
  within.cluster.ss <- 0
  for (i in 1:cn) {
    cluster.size <- sum(cluster == i)
    di <- as.dist(dmat[cluster == i, cluster == i])
    within.cluster.ss <- within.cluster.ss + sum(di^2)/cluster.size
  }
  within.cluster.ss
}





# Visualization of the output returned by the function
# NbClust()
# x : an object generated by the function NbClust()
.viz_NbClust <- function(x, print.summary = TRUE,
                         barfill = "steelblue", barcolor = "steelblue")
{
  best_nc <- x$Best.nc
  if(class(best_nc) == "numeric") print(best_nc)
  else if(class(best_nc) == "matrix"){
    best_nc <- as.data.frame(t(best_nc), stringsAsFactors = TRUE)
    best_nc$Number_clusters <- as.factor(best_nc$Number_clusters)
    
    # Summary
    if(print.summary){
      ss <- summary(best_nc$Number_clusters)
      cat ("Among all indices: \n===================\n")
      for(i in 1 :length(ss)){
        cat("*", ss[i], "proposed ", names(ss)[i], "as the best number of clusters\n" )
      }
      cat("\nConclusion\n=========================\n")
      cat("* According to the majority rule, the best number of clusters is ",
          names(which.max(ss)),  ".\n\n")
    }
    df <- data.frame(Number_clusters = names(ss), freq = ss, stringsAsFactors = TRUE )
    p <- ggpubr::ggbarplot(df,  x = "Number_clusters", y = "freq", fill = barfill, color = barcolor)+
      labs(x = "Número de k clúster", y = "Frequency among all indices",
           title = paste0("Número optimo de clúster - k = ", names(which.max(ss)) ))
    
    return(p)
  }
}


#  Determines the location of the maximum of f see ?cluster::maxSE
# +++++++++++++++++++++++++++++++++++++++++++
# f: numeric vector containing the gap statistic
# SE.f : standard error of the gap statistic
# method : character string indicating how the "optimal" number of clusters, k, 
# is computed from the gap statistics (and their standard deviations), 
# or more generally how the location k^ of the maximum of f[k] should be determined.
# SE.factor:  Determining the Número optimo de clúster, Tibshirani et al. proposed the "1 S.E."-rule.
.maxSE <- function (f, SE.f, method = c("firstSEmax", "Tibs2001SEmax", 
                                        "globalSEmax", "firstmax", "globalmax"), SE.factor = 1) 
{
  method <- match.arg(method)
  stopifnot((K <- length(f)) >= 1, K == length(SE.f), SE.f >= 
              0, SE.factor >= 0)
  fSE <- SE.factor * SE.f
  switch(method, firstmax = {
    decr <- diff(f) <= 0
    if (any(decr)) which.max(decr) else K
  }, globalmax = {
    which.max(f)
  }, Tibs2001SEmax = {
    g.s <- f - fSE
    if (any(mp <- f[-K] >= g.s[-1])) which.max(mp) else K
  }, firstSEmax = {
    decr <- diff(f) <= 0
    nc <- if (any(decr)) which.max(decr) else K
    if (any(mp <- f[seq_len(nc - 1)] >= f[nc] - fSE[nc])) which(mp)[1] else nc
  }, globalSEmax = {
    nc <- which.max(f)
    if (any(mp <- f[seq_len(nc - 1)] >= f[nc] - fSE[nc])) which(mp)[1] else nc
  })
}


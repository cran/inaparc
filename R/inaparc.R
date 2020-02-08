aldaoud <- function(x, k){ 
  if(missing(x)) 
    stop("Missing data set")
  if(is.vector(x) || is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x)) 
    stop("Data set must be a numeric vector, data frame or matrix")
  n <- nrow(x); p <- ncol(x)
  if(missing(k))
    stop("Missing input argument k")
  if(!is.numeric(k))
    stop("Argument k must be a positive integer")
  if(k < 1 || k > n)
    stop(paste0("Argument k should be between 1 and ", n, " for the data set being processed."))
  v <- matrix(nrow=k, ncol=p, 0)
  fmax <- which.max(apply(x,2,var))
  x <- as.matrix(x[order(x[,fmax]),])
  r <- floor(n/k)
  ridx <- 1
  for(i in 1:k){
    for(j in 1:p)
      v[i,j] <- median(x[ridx:(ridx+r-1),j])
    ridx <- ridx + r
  }
  colnames(v) <- colnames(x)
  rownames(v) <- paste0("Cl.", 1:k)
  result <- list()
    result$v <- v
    result$ctype <- "obj"
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}

ballhall <- function(x, k, tv){
  if(missing(x)) 
    stop("Missing data set")
  if(is.vector(x) || is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x)) 
    stop("Data set must be a numeric vector, data frame or matrix")
  n <- nrow(x); p <- ncol(x)
  if(missing(k))
    stop("Missing input argument k")
  if(!is.numeric(k))
    stop("Argument k must be a positive integer")
  if(k < 1 || k > n)
    stop(paste0("Argument k should be between 1 and ", n, " for the data set being processed."))
  if(missing(tv))
    tv = "md"
  if(!is.numeric(tv)){
    if(!is.element(tv, c("cd1", "cd2", "md", "mm")))
      stop("Argument tv must be a user specified number or 'cd1', 'cd2', 'md' or 'mm' for internal computation.")
    if(tv=="md"){
      distx <- numeric(n-1)
      for(i in 1:n-1)
        distx[i] <- .sqeucdist(x[i, 1:p], x[i+1, 1:p])
      T <- 0.5*mean(distx)/k
    }
    if(tv=="mm"){
      distx <- numeric(n-1)
      for(i in 1:n-1)
        distx[i] <- .sqeucdist(x[i, 1:p], x[i+1, 1:p])
      T <- 0.5*(max(distx)-min(distx))/k
    }
    if(tv=="cd1")
      T <- mean(diff(x[,1:p])^2)
    if(tv=="cd2")
      T <- min(diff(x[,1:p])^2 > 0)
  }
  else{
    T <- tv
  }
  v <- matrix(nrow=k, ncol=p, 0)
  xrows <- which.min(x-mean(x))/p
  x <- cbind(x, 1:n)
  if(k > 1){
    for(i in 2:k){
      xs <- x[-xrows,]
      for(j in 1:nrow(xs)){
        ridx <- 0
        for(l in xrows){
           if(.sqeucdist(xs[j,1:p], x[l, 1:p]) > T)
             ridx <- ridx + 1
        }
        if(ridx == length(xrows)){
          xrows <- c(xrows, xs[j, p+1])
          break
        }
      }
    }
  }
  v <- x[xrows, 1:p] 
  if(p==1)
    v <- as.matrix(x[xrows, 1:p]) 
  if(k > 1)
    rownames(v) <- paste0("Cl.", 1:nrow(v))
  result <- list()
    result$v <- v
    result$ctype <- "obj"
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}

crsamp <- function(x, k, r, ctype){
  if(missing(x)) 
    stop("Missing data set")
  if(is.vector(x) || is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x)) 
    stop("Data set must be a numeric vector, data frame or matrix")
  n <- nrow(x); p <- ncol(x)
  if(missing(k))
    stop("Missing input argument k")
  if(!is.numeric(k))
    stop("Argument k must be a positive integer")
  if(k < 2 || k > n)
    stop(paste0("Argument k should be between 2 and ", n, " for the data set being processed."))
  if(missing(r)) 
    r <- 2
  if(!is.numeric(r))
    stop(paste0("Argument r must be a positive integer between 2 and ", k,"."))
  if(r < 2 || r > floor(n/k))
    stop(paste0("Argument r must be a positive integer between 2 and ", floor(n/k),"."))
  if(missing(ctype)) 
    ctype <- "obj"
  if(!is.element(ctype, c("avg","med", "obj")))
    stop(paste0("'", ctype, "' is an invalid option for the centroid selection. The valid options are 'avg', 'med' or 'obj'."))
  v <- matrix(nrow=k, ncol=p, 0)
  xrows <- c()
  if(ctype=="avg"){
    for(i in 1:k){
      xrows <- sample(n, r, replace=FALSE)
      if(p > 1)
        v[i,] <- apply(x[xrows,], 2, mean)
      else
        v[i,] <- mean(x[xrows,])
    }
  }
  if(ctype=="med"){
    for(i in 1:k){
      xrows <- sample(n, r, replace=FALSE)
      if(p > 1)
        v[i,] <- apply(x[xrows,], 2, median)
      else
        v[i,] <- median(x[xrows,])
    }
  }
  if(ctype=="obj"){
    xrows <- c()
    for(i in 1:k){
      repeat{
        sampled <- sample(1:n, r, replace=FALSE)     
        if(p > 1)
          sampmean <- apply(x[sampled,], 2, mean)
        else
          sampmean <- mean(x[sampled,])
        distx <- numeric(n)
        for(j in 1:n)
          distx[j] <- .sqeucdist(x[j,], sampmean) 
        ridx <- which.min(distx)
        if(!is.element(ridx, xrows)){
          xrows <- c(xrows, ridx)
          break
        }
      }    
    }
    v <- as.matrix(x[xrows,])
  }
  rownames(v) <- paste0("Cl.", 1:k)
  colnames(v) <- colnames(x)
  result <- list()
    result$v <- v
    result$ctype <- ctype
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}

firstk <- function(x, k){ 
  if(missing(x)) 
    stop("Missing data set")
  if(is.vector(x) || is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x)) 
    stop("Data set must be a numeric vector, data frame or matrix")
  n <- nrow(x); p <- ncol(x)
  if(missing(k))
    stop("Missing input argument k")
  if(!is.numeric(k))
    stop("Argument k must be a positive integer")
  if(k < 1 || k > n)
    stop(paste0("Argument k should be between 1 and ", n, " for the data set being processed."))
  v <- matrix(nrow=k, ncol=ncol(x), x[1:k,])
  colnames(v) <- colnames(x)
  rownames(v) <- paste0("Cl.", 1:k)
  result <- list()
    result$v <- v
    result$ctype <- "obj"
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}

forgy <- function(x, k){
  if(missing(x)) 
    stop("Missing data set")
  if(is.vector(x) || is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x)) 
    stop("Data set must be a numeric vector, data frame or matrix")
  n <- nrow(x); p <- ncol(x)
  if(missing(k))
    stop("Missing input argument k")
  if(!is.numeric(k))
    stop("Argument k must be a positive integer")
  if(k < 1 || k > nrow(x))
    stop(paste0("Argument k should be between 1 and ", n, " for the data set being processed."))
  v <- matrix(nrow=k, ncol=p, 0)
  cidx <- sample(rep(seq_len(k), each = ceiling(n/k)), size=n)
  for(i in 1:k)
    for(j in 1:p)
      v[i,j] <- mean(x[cidx == i,j])
  colnames(v) <- colnames(x)
  rownames(v) <- paste0("Cl.", 1:k)
  result <- list()
    result$v <- v
    result$ctype <- "avg"
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}

hartiganwong <- function(x, k){ 
  if(missing(x)) 
    stop("Missing data set")
  if(is.vector(x) || is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x)) 
    stop("Data set must be a numeric vector, data frame or matrix")
  n <- nrow(x); p <- ncol(x)
  if(missing(k))
    stop("Missing input argument k")
  if(!is.numeric(k))
    stop("Argument k must be a positive integer")
  if(k < 1 || k > n)
    stop(paste0("Argument k should be between 1 and ", n, " for the data set being processed."))
  v <- matrix(nrow=k, ncol=p, 0)
  xmean <- apply(x, 2, mean) 
  x <- cbind(x, rowSums(abs(x-xmean)))
  x <- x[order(x[,p+1]),]
  r <- floor(n/k)
  for(i in 1:k){
    v[i,] <- x[1+(i-1)*r, 1:p]
  }
  colnames(v) <- colnames(x[,1:p])
  rownames(v) <- paste0("Cl.", 1:k)
  result <- list()
    result$v <- v
    result$ctype <- "obj"
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}

inofrep <- function(x, k, sfidx, sfpm, binrule, nbins, tcmethod, tc){
  if(missing(x)) 
    stop("Missing data set")
  if(is.vector(x) || is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x)) 
    stop("Data set must be a numeric vector, data frame or matrix")
  n <- nrow(x); p <- ncol(x)
  if(missing(k))
    stop("Input argument k is missing")
  if(!is.numeric(k))
    stop("Input argument k must be a positive integer")
  else
    k <- as.integer(k)
  if(k<1 || k > n)
    stop(paste0("k should be between 1 and ", n, "."))
  if(missing(binrule) & !missing(nbins)) 
    binrule <- "usr"
  if(missing(binrule) & missing(nbins)) 
    binrule <- "sqr"
  if(!binrule %in% c("sqr", "scott", "sturges", "huntsberger", "bc", "cencov", "fd", "doane", "cebeci","usr")) 
    stop("The valid rules for class generation are 'sqr', 'scott', 'sturges', 'huntsberger', 'bc', 'cencov', 'fd', 'doane', 'cebeci','usr'")
  if(binrule=="usr" & missing(nbins))
    nbins <- 20
  if(missing(tcmethod)) 
    tcmethod <- "min"
  if(missing(tcmethod) & !missing(tc)) 
    tcmethod <- "usr"
  if(!tcmethod %in% c("sd1", "sd2", "q1", "iqr", "avg", "min", "min2", "log2","usr")) 
    stop("The valid of threshold frequency are 'sd1', 'sd2', 'q1', 'iqr', 'avg', 'min', 'min2', 'log2' and 'usr'")
  if(!missing(tc)){  
    if(!is.numeric(tc))
      stop("tc, threshold frequency value should be 0 or a positive integer")
    else
      tc <- as.integer(tc)
    if(tc<0 || tc>nrow(x))
      stop("tc, threshold frequency value should be 0 or a positive integer less than the number of rows of data set")
  }else{
    if(tcmethod=="usr" & missing(tc))
      tc <- 1
  }
  if(missing(sfidx)){
    npeaks <- 0
    for(j in 1:p){
      hvals <- kpeaks::genpolygon(x[,j], binrule=binrule, nbins=nbins)
      resfpp <- kpeaks::findpolypeaks(hvals$mids, hvals$freqs, tcmethod=tcmethod, tc=tc)
      if(npeaks < resfpp$np){
        npeaks <- resfpp$np
        sfidx <- j
      }
    }
  }
  if(!is.numeric(sfidx))
    stop("The column index of selected feature must be a positive integer")
  sfidx <- as.integer(sfidx)
  if(sfidx < 1 || sfidx > p)
    stop("The column index of selected feature cannot be less than 1 or greater than ", p, " for this data set" )
  nuniqx <- length(unique(x[,sfidx]))
  if(k > nuniqx)
    stop(paste0("k should be less than ", nuniqx, " because there are only ", nuniqx, " distinct values of the selected feature for avoiding the coincided clusters."))
  if(!missing(sfpm)){
    if(is(sfpm, "data.frame"))
      sfpm <- as.matrix(sfpm)
    if(!is(sfpm, "matrix")) 
     stop("Peak matrix of the selected feature must be a numeric data frame or matrix")
  } 
  else{
    hvals <- kpeaks::genpolygon(x[,sfidx], binrule=binrule, nbins=nbins)
    sfpm <- kpeaks::findpolypeaks(hvals$mids, hvals$freqs, tcmethod=tcmethod, tc=tc)$pm
  }
  nsampled <- nrow(sfpm)
  if(nsampled > 1){
    sfpm <- sfpm[order(sfpm[,2], decreasing = TRUE),]
    if(nsampled > k )
      nsampled <- k
  }
  v <- matrix(nrow=k, ncol=p, 0)
  xsampled <- c()
  for(i in 1:nsampled){
    ridx <- which.min(abs(x[,sfidx] - sfpm[i,1]))
    xsampled <- c(xsampled, ridx)
  }
  if(k > nsampled){
    g <- 1:n
    xremained <- g[-xsampled]
    xsampled <- c(xsampled, sample(xremained, k-nsampled, replace=FALSE))
  }
  v <- x[xsampled,]
  if(p==1){
    v <- as.matrix(v)
    rownames(v) <- paste0("Cl.", 1:k)
  } 
  else if(k > 1){
    colnames(v) <- colnames(x)
    rownames(v) <- paste0("Cl.", 1:k)
  }
  result <- list()
    result$v <- v
    result$sfidx <- sfidx
    result$ctype <- "obj"
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}

inscsf <- function(x, k, sfidx, ctype){
  if(missing(x)) 
    stop("Missing data set")
  if(is.vector(x) || is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x)) 
    stop("Data set must be a numeric vector, data frame or matrix")
  n <- nrow(x); p <- ncol(x)
  if(missing(k))
    stop("Missing k, number of clusters")
  if(!is.numeric(k))
    stop("k must be a positive integer")
  else
    k <- as.integer(k)
  if(k<1 || k > nrow(x))
    stop(paste0("k should be between 1 and ", n, " for the data set being processed."))
  if(!missing(sfidx)){
    if(!is.numeric(sfidx))
      stop(paste0("Selected feature index must be an integer between 1 and ", p, " for the data set being processed."))
    else
      sfidx <- as.integer(sfidx)
   }else{
    sfidx <- 1
    if(p > 1){
      nuniqx <- length(unique(x[,1]))
      for(j in 2:p){
        if(length(unique(x[,j])) > nuniqx){
          nuniqx <- length(unique(x[,j]))
          sfidx <- j
        }
      }
    }
  }
  if(sfidx < 1 |sfidx > p)
    stop(paste0("Selected feature index must be an integer between 1 and ", p, " for the data set being processed."))
  xuniques <- unique(x[,sfidx])
  modex <- xuniques[which.max(tabulate(match(x[,sfidx], xuniques)))]
  meanx <- mean(x[,sfidx])
  medx <- median(x[,sfidx])
  if(missing(ctype)){
    if(modex > medx) 
      ctype <- "mod"
    else if(medx > meanx) 
        ctype <- "med"
      else 
        ctype <- "avg"
  }
  if(!is.element(ctype, c("avg","med", "mod")))
    stop(paste0("'", ctype, "' is an invalid option for the center selection. Enter 'avg', 'med' or 'mod' for type selection."))
  if(ctype=="avg") center <- meanx
  if(ctype=="med") center <- medx
  if(ctype=="mod") center <- modex
  r1 <- (max(x[,sfidx])-center)/(k-1)/2
  r2 <- (center-min(x[,sfidx]))/(k-1)/2
  nuniqx <- length(unique(x[,sfidx]))
  if(k > nuniqx)
    stop(paste0("For avoiding the coincided clusters, k should be less than ", nuniqx, " because there are only ", nuniqx, " distinct values of the selected feature."))
  v <- matrix(nrow=k, ncol=p, 0)
  if(p > 1){
    ridx <- which.min(abs(x[,sfidx] - center))
    v[1,] <- x[ridx,]
    if(k > 1){
      for(i in 2:k){
        x <- x[-ridx,]
        if(i%%2==1) 
          ridx <- which.min(abs(x[,sfidx] - (center + (i-1) * r1)))
        else
          ridx <- which.min(abs(x[,sfidx] - (center - (i-1) * r2)))
        v[i,] <- x[ridx,] 
      }
    }
  }else{
    ridx <- which.min(abs(x - center))
    v[1,] <- x[ridx]
    if(k > 1){
      for(i in 2:k){
        x <- x[-ridx]
        if(i%%2==1) 
          ridx <- which.min(abs(x - (center + (i-1) * r1)))
        else
          ridx <- which.min(abs(x - (center - (i-1) * r2)))
        v[i,] <- x[ridx] 
      }
    }
  }
  if(p==1){
    v <- as.matrix(v)
    rownames(v) <- paste0("Cl.", 1:k)
  } 
  else if(k > 1){
    colnames(v) <- colnames(x)
    rownames(v) <- paste0("Cl.", 1:k)
  }
  result <- list()
    result$v <- v
    result$sfidx <- sfidx
    result$ctype <- "obj"
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}	

insdev <- function(x, k, sfidx){
  if(missing(x)) 
    stop("Missing data set")
  if(is.vector(x) || is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x)) 
    stop("Data set must be a numeric vector, data frame or matrix")
  n <- nrow(x); p <- ncol(x)
  if(missing(k))
    stop("Missing k, number of clusters")
  if(!is.numeric(k))
    stop("Argument k must be a positive integer")
  if(k<1 || k > n)
    stop(paste0("k should be between 1 and ", nrow(x), " for the data set being processed."))
  if(!missing(sfidx)){
    if(!is.numeric(sfidx))
      stop(paste0("Selected feature index must be an integer between 1 and ", p, " for the data set being processed."))
    else 
      sfidx <- round(sfidx)
  }else{
    sfidx <- 1
    if(p > 1){
      cvx <- sqrt(var(x[,1]))/mean(x[,1])
      for(j in 2:p){
        if(sqrt(var(x[,j]))/mean(x[,j]) > cvx){
          cvx <- sqrt(var(x[,j]))/mean(x[,j])
          sfidx <- j
        }
      }
    }  
  }
  if(sfidx < 1 || sfidx > p)
    stop(paste0("Selected feature index must be an integer between 1 and ", p, " for the data set being processed."))
  v <- matrix(nrow=k, ncol=p, 0)
  ridx <- which.min(abs(x[,sfidx] - mean(x[,sfidx])))
  xinc <- sqrt(var(x[,sfidx]))/k/2
  v[1,] <- x[ridx,]
  if(k > 1){
    for(i in 2:k){
      if(p > 1)
        x <- x[-ridx,]
      else
        x <- x[-ridx]
      if(p > 1){
        if(i%%2==1) 
          ridx <- which.min(abs(x[,sfidx] - (v[1,sfidx] + (i-1) * xinc)))
        else
          ridx <- which.min(abs(x[,sfidx] - (v[1,sfidx] - i * xinc)))
        v[i,] <- x[ridx,] 
      }else{
        if(i%%2==1) 
          ridx <- which.min(abs(x - (v[1,sfidx] + (i-1) * xinc)))
        else
          ridx <- which.min(abs(x - (v[1,sfidx] - i * xinc)))
        v[i,] <- x[ridx] 
      }
    }
  }
  if(p==1){
    v <- as.matrix(v)
    rownames(v) <- paste0("Cl.", 1:k)
  } 
  else if(k > 1){
    colnames(v) <- colnames(x)
    rownames(v) <- paste0("Cl.", 1:k)
  }
  result <- list()
    result$v <- v
    result$sfidx <- sfidx
    result$ctype <- "obj"
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}	

kkz <- function(x, k){
  if(missing(x)) 
    stop("Missing data set")
  if(is.vector(x) || is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x)) 
    stop("Data set must be a numeric vector, data frame or matrix")
  n <- nrow(x); p <- ncol(x)
  if(missing(k))
    stop("Missing input argument k")
  if(!is.numeric(k))
    stop("Argument k must be a positive integer")
  if(k < 1 || k > nrow(x))
    stop(paste0("Argument k should be between 1 and ", n, " for the data set being processed."))
  v <- matrix(nrow=k, ncol=p, 0)
  ix <- numeric(k); dx <- numeric(k)
  ridx <- k + 1
  if(k == n){
    xrows <- 1:k
  }
  else{
    norms <- numeric(n)
    for(i in 1:n){
      norms[i] <- sum(.sqeucdist(x[i,], x[-x[i,],]))
    }
    xrows <- which.max(norms)
    v[1,] <- x[xrows,]
    if(k >= 2 && k < n){
      ridx <- which.max(.sqeucdist(x[-xrows,], v[1,]))
      v[2,] <- x[ridx,]
      xrows <- c(xrows, ridx)
      ridx <- 3
      x <- cbind(x, 1:n)
    }
  }
  while(ridx <= k){
    ix <- rep(0, k); dx <- rep(0, k)
    xs <- x[-xrows,]
    for(j in 1:length(xrows)){
      mindist <- Inf
      for(i in 1:nrow(xs)){
        distx <- .sqeucdist(xs[i, 1:p], x[xrows[j], 1:p])
        if(distx < mindist){
          mindist <- distx
          ix[j] <- xs[i, p+1]
          dx[j] <- mindist
          if(mindist == 0)
            break
        }
      }
    }
    if(sum(ix) == 0)
      break
    xrows <- c(xrows, ix[which.max(dx)])
    ridx <- ridx + 1
  }
  if(k > 1)
    v <- as.matrix(x[xrows, 1:p]) 
  colnames(v) <- colnames(x[,1:p])
  rownames(v) <- paste0("Cl.", 1:k)
  result <- list()
    result$v <- v
    result$ctype <- "obj"
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}

kmpp <- function(x, k){
  if(missing(x)) 
    stop("Missing data set")
  if(is.vector(x) || is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x)) 
    stop("Data set must be a numeric vector, data frame or matrix")
  n <- nrow(x); p <- ncol(x)
  if(missing(k))
    stop("Missing input argument k")
  if(!is.numeric(k))
    stop("Argument k must be a positive integer")
  if(k < 1 || k > n)
    stop(paste0("Argument k should be between 1 and ", n, " for the data set being processed."))
  v <- matrix(nrow=k, ncol=p, 0)
  xrows <- numeric(k) 
  distmat <- matrix(nrow=n, ncol=k-1, 0) 
  prbx <- rep(1, n)
  if(k == n)
    xrows <- 1:k
  else{ 
    if(k > 1){
      for(i in 1:(k-1)){
        xrows[i] <- sample.int(n, 1, prob = prbx) 
        distmat[,i] <- colSums((t(x)-x[xrows[i],])^2) 
        prbx <- distmat[cbind(1:n, max.col(-distmat[, 1:i, drop = FALSE]))]
      }
    }
    xrows[k] <- sample.int(n, 1, prob = prbx)
  }
  v <- x[xrows, ]
  if(k > 1){
    v <- as.matrix(x[xrows, ])
    colnames(v) <- colnames(x)
    rownames(v) <- paste0("Cl.", 1:k)
  }
  result <- list()
    result$v <- v
    result$ctype <- "obj"
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}

ksegments <- function(x, k, ctype){ 
  if(missing(x)) 
    stop("Missing data set")
  if(is.vector(x) || is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x)) 
    stop("Data set must be a numeric vector, data frame or matrix")
  n <- nrow(x); p <- ncol(x)
  if(missing(k))
    stop("Missing input argument k")
  if(!is.numeric(k))
    stop("Argument k must be a positive integer")
  if(k < 1 || k > nrow(x))
    stop(paste0("Argument k should be between 1 and ", nrow(x), " for the data set being processed."))
  if(missing(ctype)) 
    ctype <- "avg"
  if(!is.element(ctype, c("avg","med")))
    stop(paste0("'", ctype, "' is an invalid option for the centroid selection. The valid options are 'avg' or 'med'."))
  v <- matrix(nrow=k, ncol=p, 0)
  ns <- floor(n/k)
  step <- 1
  for(i in 1:k){
    sridx <- step
    eridx <- sridx+ns-1
    if(ctype=="avg"){
      if(ns >= 2){
        if(p > 1)
          v[i,] <- apply(x[sridx:eridx,], 2, mean)
        else 
          v[i,] <- mean(x[sridx:eridx,])
      }else
        v[i,] <- x[sridx:eridx,]
    }
    else{
       if(ns >= 2){
        if(p > 1)
          v[i,] <- apply(x[sridx:eridx,], 2, median)
        else 
          v[i,] <- median(x[sridx:eridx,])
      }else
        v[i,] <- x[sridx:eridx,]
    }
    step <- eridx + 1
  }
  colnames(v) <- colnames(x)
  rownames(v) <- paste0("Cl.", 1:k)
  result <- list()
    result$v <- v
    result$ctype <- ctype
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}

ksteps <- function(x, k, ctype){ 
  if(missing(x)) 
    stop("Missing data set")
  if(is.vector(x) || is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x)) 
    stop("Data set must be a numeric vector, data frame or matrix")
  n <- nrow(x); p <- ncol(x)
  if(missing(k))
    stop("Missing input argument k")
  if(!is.numeric(k))
    stop("Argument k must be a positive integer")
  if(k < 1 || k > as.integer(sqrt(nrow(x))))
    stop(paste0("Argument k should be between 1 and ", as.integer(sqrt(nrow(x))), " for the data set being processed."))
  if(missing(ctype)) 
    ctype <- "avg"
  if(!is.element(ctype, c("avg","med")))
    stop(paste0("'", ctype, "' is an invalid option for the centroid selection. The valid options are 'avg' or 'med'."))
  v <- matrix(nrow=k, ncol=p, 0)
  xrows <- c()
  for(i in 1:k){
    for(j in 1:k){
      xrows <- c(xrows, i+(j-1)*k)
    }
    if(ctype=="avg"){
      if(p > 1){
        if(k > 1)
          v[i,] <- apply(x[xrows,], 2, mean)
        else
          v[i,] <- x[xrows,]
      }
      else 
        v[i,] <- mean(x[xrows,])
    }
    if(ctype=="med"){
      if(p > 1){
        if(k > 1)
          v[i,] <- apply(x[xrows,], 2, median)
        else
          v[i,] <- x[xrows,]
      }
      else 
        v[i,] <- median(x[xrows,])
    }
  }
  if(p==1){
    v <- as.matrix(v)
    rownames(v) <- paste0("Cl.", 1:k)
  } 
  else if(k > 1){
    colnames(v) <- colnames(x)
    rownames(v) <- paste0("Cl.", 1:k)
  }
  result <- list()
    result$v <- v
    result$ctype <- ctype
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}

lastk <- function(x, k){ 
  if(missing(x)) 
    stop("Missing data set")
  if(is.vector(x) || is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x)) 
    stop("Data set must be a numeric vector, data frame or matrix")
  n <- nrow(x); p <- ncol(x)
  if(missing(k))
    stop("Missing input argument k")
  if(!is.numeric(k))
    stop("Argument k must be a positive integer")
  if(k < 1 || k > nrow(x))
    stop(paste0("Argument k should be between 1 and ", n, " for the data set being processed."))
  v <- matrix(nrow=k, ncol=p, x[(n-k+1):n,])
  colnames(v) <- colnames(x)
  rownames(v) <- paste0("Cl.", 1:k)
  result <- list()
    result$v <- v
    result$ctype <- "obj"
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}

lhsmaximin <- function(x, k, ncp){ 
  if(missing(x)) 
    stop("Missing data set")
  if(is.vector(x) || is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x)) 
    stop("Data set must be a numeric vector, data frame or matrix")
  n <- nrow(x); p <- ncol(x)
  if(missing(k))
    stop("Missing input argument k")
  if(!is.numeric(k))
    stop("Argument k must be a positive integer")
  if(k < 1 || k > nrow(x))
    stop(paste0("Argument k should be between 1 and ", n, " for the data set being processed."))
  if(missing(ncp))
    ncp <- 1
  if(!is.numeric(ncp))
    stop("Argument dup must be a positive integer")
  if(ncp < 1 || ncp > k)
    stop(paste0("Argument ncp should be between 1 and ", k, "."))
  v <- matrix(nrow=k, ncol=p, 0)
  xmean <- numeric(p)
  xsd <- numeric(p)
  for(j in 1:p){
    xmean[j] <- mean(x[,j])
    xsd[j] <- sqrt(var(x[,j]))
  }
  z <- lhs::maximinLHS(k, p, dup=ncp)
  for(j in 1:p)
    v[,j] <- qnorm(z[,j], mean=xmean[j], sd=xsd[j])
  if(p==1){
    v <- as.matrix(v)
    rownames(v) <- paste0("Cl.", 1:k)
  } 
  else if(k > 1){
    colnames(v) <- colnames(x)
    rownames(v) <- paste0("Cl.", 1:k)
  }
  result <- list()
    result$v <- v
    result$ctype <- "obj"
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}

lhsrandom <- function(x, k){ 
  if(missing(x)) 
    stop("Missing data set")
  if(is.vector(x) || is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x)) 
    stop("Data set must be a numeric vector, data frame or matrix")
  n <- nrow(x); p <- ncol(x)
  if(missing(k))
    stop("Missing input argument k")
  if(!is.numeric(k))
    stop("Argument k must be a positive integer")
  if(k < 1 || k > nrow(x))
    stop(paste0("Argument k should be between 1 and ", n, " for the data set being processed."))
  v <- matrix(nrow=k, ncol=p, 0)
  xmean <- numeric(p)
  xsd <- numeric(p)
  for(j in 1:p){
    xmean[j] <- mean(x[,j])
    xsd[j] <- sqrt(var(x[,j]))
  }
  z <- lhs::randomLHS(k, p)
  for(j in 1:p)
    v[,j] <- qnorm(z[,j], mean=xmean[j], sd=xsd[j])
  if(p==1){
    v <- as.matrix(v)
    rownames(v) <- paste0("Cl.", 1:k)
  } 
  else if(k > 1){
    colnames(v) <- colnames(x)
    rownames(v) <- paste0("Cl.", 1:k)
  }
  result <- list()
    result$v <- v
    result$ctype <- "obj"
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}

maximin <- function(x, k){
  if(missing(x)) 
    stop("Missing data set")
  if(is.vector(x) || is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x)) 
    stop("Data set must be a numeric vector, data frame or matrix")
  n <- nrow(x); p <- ncol(x)
  if(missing(k))
    stop("Missing input argument k")
  if(!is.numeric(k))
    stop("Argument k must be a positive integer")
  if(k < 1 || k > n)
    stop(paste0("Argument k should be between 1 and ", n, " for the data set being processed."))
  v <- matrix(nrow=k, ncol=p, 0)
  ix <- numeric(k); dx <- numeric(k)
  ridx <- k + 1
  if(k == n)
    xrows <- 1:k
  else{
    xrows <-  sample(n, 1)
    v[1,] <- x[xrows,]
    if(k >= 2 && k < n){
      ridx <- which.max(.sqeucdist(x[-xrows,], v[1,]))
      v[2,] <- x[ridx,]
      xrows <- c(xrows, ridx)
      ridx <- 3
      x <- cbind(x, 1:n)
    }
  }
  while(ridx <= k){
    ix <- rep(0, k); dx <- rep(0, k)
    xs <- x[-xrows,]
    for(j in 1:length(xrows)){
      mindist <- Inf
      for(i in 1:nrow(xs)){
        distx <- .sqeucdist(xs[i, 1:p], x[xrows[j],1:p])
        if(distx < mindist){
          mindist <- distx
          ix[j] <- xs[i, p+1]
          dx[j] <- mindist
          if(mindist == 0)
            break
        }
      }
    }
    if(sum(ix) == 0)
      break
    xrows <- c(xrows, ix[which.max(dx)])
    ridx <- ridx + 1
  }
  v <- x[xrows, 1:p] 
  if(p==1){
    v <- as.matrix(v)
    rownames(v) <- paste0("Cl.", 1:k)
  } 
  else if(k > 1){
    colnames(v) <- colnames(x[,1:p])
    rownames(v) <- paste0("Cl.", 1:k)
  }
  result <- list()
    result$v <- v
    result$ctype <- "obj"
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}

mscseek <- function(x, k, tv){
  if(missing(x)) 
    stop("Missing data set")
  if(is.vector(x) || is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x)) 
    stop("Data set must be a numeric vector, data frame or matrix")
  n <- nrow(x); p <- ncol(x)
  if(missing(k))
    stop("Missing input argument k")
  if(!is.numeric(k))
    stop("Argument k must be a positive integer")
  if(k < 1 || k > nrow(x))
    stop(paste0("Argument k should be between 1 and ", n, " for the data set being processed."))
  if(missing(tv))
    tv = "cd1"
  if(!is.numeric(tv)){
    if(!is.element(tv, c("cd1", "cd2", "md", "mm")))
      stop("Argument tv must be a user specified number or 'cd1', 'cd2', 'md' or 'mm' for internal computation.")
    if(tv=="md"){
      distx <- numeric(n-1)
      for(i in 1:n-1)
        distx[i] <- .sqeucdist(x[i, 1:p], x[i+1, 1:p])
      T <- 0.5*mean(distx)/k
    }
    if(tv=="mm"){
      distx <- numeric(n-1)
      for(i in 1:n-1)
        distx[i] <- .sqeucdist(x[i, 1:p], x[i+1, 1:p])
      T <- 0.5*(max(distx)-min(distx))/k
    }
    if(tv=="cd1")
      T <- mean(diff(x[,1:p])^2)
    if(tv=="cd2")
      T <- min(diff(x[,1:p])^2 > 0)
  }
  else{
    T <- tv
  }
  v <- matrix(nrow=k, ncol=p, 0)
  x <- cbind(x, 1:n)
  xrows <- c(1)
  i <- 2
  while(i <= k){
    xs <- x[-xrows,]
    j <- 1
    while(j <= nrow(xs)){
      z <- which(xs[, p+1] == sample(xs[,p+1], 1))
      ridx <- 0
      for(l in xrows){
        if(.sqeucdist(xs[z, 1:p], x[l, 1:p]) > T)
          ridx <- ridx + 1
      }
      if(ridx == length(xrows)){
        xrows <- c(xrows, xs[z, p+1])
        break
      }
      j <- j + 1
    }
    mindist <- T
    for(z1 in xrows[1:length(xrows)-1]){
      for(z2 in xrows[2:length(xrows)]){
        distx <- .sqeucdist(x[z1,1:p], x[z2,1:p])
        if(distx < mindist && mindist > 0)
          mindist <- distx
      }
    }
    T <- mindist
    i <- i + 1
  }
  v <- x[xrows, 1:p] 
  if(p==1){
    v <- as.matrix(v)
    rownames(v) <- paste0("Cl.", 1:k)
  } 
  else if(k > 1){
    colnames(v) <- colnames(x[,1:p])
    rownames(v) <- paste0("Cl.", 1:k)
  }
  result <- list()
    result$v <- v
    result$ctype <- "obj"
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}

rsamp <- function(x, k){
  if(missing(x)) 
    stop("Missing data set")
  if(is.vector(x) || is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x)) 
    stop("Data set must be a numeric vector, data frame or matrix")
  n <- nrow(x); p <- ncol(x)
  if(missing(k))
    stop("Missing input argument k")
  if(!is.numeric(k))
    stop("Argument k must be a positive integer")
  if(k < 1 || k > nrow(x))
    stop(paste0("Argument k should be between 1 and ", n, " for the data set being processed."))
  v <- x[sample(1:n, k, replace=FALSE),]  
  if(p==1){
    v <- as.matrix(v)
    rownames(v) <- paste0("Cl.", 1:k)
  } 
  else if(k > 1){
    colnames(v) <- colnames(x[,1:p])
    rownames(v) <- paste0("Cl.", 1:k)
  }
  result <- list()
    result$v <- v
    result$ctype <- "obj"
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}

rsegment <- function(x, k){ 
  if(missing(x)) 
    stop("Missing data set")
  if(is.vector(x) || is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x)) 
    stop("Data set must be a numeric vector, data frame or matrix")
  n <- nrow(x); p <- ncol(x)
  if(missing(k))
    stop("Missing input argument k")
  if(!is.numeric(k))
    stop("Argument k must be a positive integer")
  if(k < 1 || k > nrow(x))
    stop(paste0("Argument k should be between 1 and ", n, " for the data set being processed."))
  v <- matrix(nrow=k, ncol=p, 0)
  if((n-k) == 0)
    sridx <- 1
  else
    sridx <- sample(1:(n-k), 1)
  eridx <- sridx+k-1
  v <- x[sridx:eridx,]
  if(p==1){
    v <- as.matrix(v)
    rownames(v) <- paste0("Cl.", 1:k)
  } 
  else if(k > 1){
    colnames(v) <- colnames(x[,1:p])
    rownames(v) <- paste0("Cl.", 1:k)
  }
  result <- list()
    result$v <- v
    result$ctype <- "obj"
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}

scseek <- function(x, k, tv){
  if(missing(x)) 
    stop("Missing data set")
  if(is.vector(x) || is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x)) 
    stop("Data set must be a numeric vector, data frame or matrix")
  n <- nrow(x); p <- ncol(x)
  if(missing(k))
    stop("Missing input argument k")
  if(!is.numeric(k))
    stop("Argument k must be a positive integer")
  if(k < 1 || k > n)
    stop(paste0("Argument k should be between 1 and ", n, " for the data set being processed."))
  if(missing(tv))
    tv = "md"
  if(!is.numeric(tv)){
    if(!is.element(tv, c("cd1", "cd2", "md", "mm")))
      stop("Argument tv must be a user specified number or 'cd1', 'cd2', 'md' or 'mm' for internal computation.")
    if(tv=="cd1"){
      distx <- numeric(n-1)
      for(i in 1:n-1)
        distx[i] <- .sqeucdist(x[i, 1:p], x[i+1, 1:p])
      T <- 0.5*mean(distx)/k
    }
    if(tv=="mm"){
      distx <- numeric(n-1)
      for(i in 1:n-1)
        distx[i] <- .sqeucdist(x[i, 1:p], x[i+1, 1:p])
      T <- 0.5*(max(distx)-min(distx))/k
    }
    if(tv=="cd1")
      T <- mean(diff(x[,1:p])^2)
    if(tv=="cd2")
      T <- min(diff(x[,1:p])^2 > 0)
  }
  else{
    T <- tv
  }
  v <- matrix(nrow=k, ncol=p, 0)
  x <- cbind(x, 1:n)
  xrows <- c(1)
  i <- 2
  while(i <= k){
    xs <- x[-xrows,]
    j <- 1
    while(j <= nrow(xs)){
      z <- which(xs[, p+1] == sample(xs[,p+1], 1))
      ridx <- 0
      for(l in xrows){
        if(.sqeucdist(xs[z, 1:p], x[l, 1:p]) > T)
          ridx <- ridx + 1
      }
      if(ridx == length(xrows)){
        xrows <- c(xrows, xs[z, p+1])
        break
      }
      j <- j + 1
    }
    i <- i + 1
  }
  v <- x[xrows, 1:p] 
  if(p==1){
    v <- as.matrix(v)
    rownames(v) <- paste0("Cl.", 1:nrow(v))
  } 
  else if(k > 1){
    colnames(v) <- colnames(x[,1:p])
    rownames(v) <- paste0("Cl.", 1:nrow(v))
  }
  result <- list()
    result$v <- v
    result$ctype <- "obj"
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}

scseek2 <- function(x, k, sfidx, tv){
  if(missing(x)) 
    stop("Missing data set")
  if(is.vector(x) || is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x)) 
    stop("Data set must be a numeric vector, data frame or matrix")
  n <- nrow(x); p <- ncol(x)
  if(missing(k))
    stop("Missing input argument k")
  if(!is.numeric(k))
    stop("Argument k must be a positive integer")
  if(k < 1 || k > n)
    stop(paste0("Argument k should be between 1 and ", n, " for the data set being processed."))
  if(!missing(sfidx)){
    if(!is.numeric(sfidx))
      stop(paste0("Selected feature index must be an integer between 1 and ", p, " for the data set being processed."))
    else 
      sfidx <- round(sfidx)
    if((sfidx < 1) || (sfidx > p))
      stop(paste0("Selected feature index must be an integer between 1 and ", p, " for the data set being processed."))
  }else{
    sfidx <- 1
    if(p > 1){
      cvx <- sqrt(var(x[,1]))/mean(x[,1])
      for(j in 2:p){
        if(sqrt(var(x[,j]))/mean(x[,j]) > cvx){
          cvx <- sqrt(var(x[,j]))/mean(x[,j])
          sfidx <- j
        } 
      }
    } 
  }
  if(missing(tv))
    tv = "cd1"
  if(!is.numeric(tv)){
    if(!is.element(tv, c("cd1", "cd2", "md", "mm")))
      stop("Argument tv must be a user specified number or 'cd1', 'cd2', 'md' or 'mm' for internal computation.")
    if(tv=="md"){
      distx <- numeric(n-1)
      for(i in 1:n-1)
        distx[i] <- .sqeucdist(x[i, 1:p], x[i+1, 1:p])
      T <- 0.5*mean(distx)/k
    }
    if(tv=="mm"){
      distx <- numeric(n-1)
      for(i in 1:n-1)
        distx[i] <- .sqeucdist(x[i, 1:p], x[i+1, 1:p])
      T <- 0.5*(max(distx)-min(distx))/k
    }
    if(tv=="cd1"){
      T <- mean(diff(x[,1:p])^2)
    }
    if(tv=="cd2"){
      T <- min(diff(x[,1:p])^2 > 0)
    }
  }
  else {
    T <- tv
  }
  v <- matrix(nrow=k, ncol=p, 0)
  x <- cbind(x, 1:n)
  xrows <- c(1)
  i <- 2
  while(i <= k){
    xs <- x[-xrows,]
    j <- 1
    while(j <= nrow(xs)){
      z <- which(xs[, p+1] == sample(xs[,p+1], 1))
      ridx <- 0
      for(l in xrows){
        if(.sqeucdist(xs[z, sfidx], x[l, sfidx]) > T)
          ridx <- ridx + 1
      }
      if(ridx == length(xrows)){
        xrows <- c(xrows, xs[z, p+1])
        break
      }
      j <- j + 1
    }
    i <- i + 1
  }
  v <- x[xrows, 1:p] 
  if(p==1){
    v <- as.matrix(x[xrows, 1:p])
    rownames(v) <- paste0("Cl.", 1:nrow(v))
  } 
  else if(k > 1){
    colnames(v) <- colnames(x[,1:p])
    rownames(v) <- paste0("Cl.", 1:nrow(v))
  }
  result <- list()
    result$v <- v
    result$sfidx <- sfidx
    result$ctype <- "obj"
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}

spaeth <- function(x, k){
  if(missing(x)) 
    stop("Missing data set")
  if(is.vector(x) || is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x)) 
    stop("Data set must be a numeric vector, data frame or matrix")
  n <- nrow(x); p <- ncol(x)
  if(missing(k))
    stop("Missing input argument k")
  if(!is.numeric(k))
    stop("Argument k must be a positive integer")
  if(k < 1 || k > n){
    errmsg <- paste0("Argument k should be between 1 and ", n, " for the data set being processed.")
    stop(errmsg)
  }
  v <- matrix(nrow=k, ncol=p, 0)
  xrows <- numeric(n)
  for(j in 1:n)
    xrows[j] <- (j-1) %% k + 1
  for(i in 1:k)
    for(j in 1:p)
      v[i,j] <- mean(x[xrows == i,j])
  colnames(v) <- colnames(x)
  rownames(v) <- paste0("Cl.", 1:k)
  result <- list()
    result$v <- v
    result$ctype <- "avg"
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}

ssamp <- function(x, k){
  if(missing(x)) 
    stop("Missing data set")
  if(is.vector(x) || is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x)) 
    stop("Data set must be a numeric vector, data frame or matrix")
  n <- nrow(x); p <- ncol(x)
  if(missing(k))
    stop("Missing input argument k")
  if(!is.numeric(k))
    stop("Argument k must be a positive integer")
  if(k < 1 || k > n){
    errmsg <- paste0("Argument k should be between 1 and ", n, " for the data set being processed.")
    stop(errmsg)
  }
  v <- matrix(nrow=k, ncol=p, 0)
  colnames(v) <- colnames(x)
  xinc <- floor(n/k)
  rid1 <- runif(1, 1:(xinc))
  xrows <- c()
  for(i in 1:k){
    ridx <- rid1 + xinc * (i-1)
    xrows <- c(xrows, ridx)
  }
  v <- x[xrows,]
  if(k > 1){
    v <- as.matrix(x[xrows,])
    rownames(v) <- paste0("Cl.", 1:k)
  }
  result <- list()
    result$v <- v
    result$ctype <- "obj"
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}

topbottom <- function(x, k){ 
  if(missing(x)) 
    stop("Missing data set")
  if(is.vector(x) || is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x)) 
    stop("Data set must be a numeric vector, data frame or matrix")
  n <- nrow(x); p <- ncol(x)
  if(missing(k))
    stop("Missing input argument k")
  if(!is.numeric(k))
    stop("Argument k must be a positive integer")
  if(k < 1 || k > n){
    errmsg <- paste0("Argument k should be between 1 and ", n, " for the data set being processed.")
    stop(errmsg)
  }
  v <- matrix(nrow=k, ncol=p, 0)
  i1 <- 1; i2 <- nrow(x)
  for(i in 1:k){
     if(i%%2==1){
       v[i,] <- x[i1,]
       i1 <- i1 + 1
     } 
     else{
       v[i,] <- x[i2,]
       i2 <- i2 - 1
     }
  }
  colnames(v) <- colnames(x)
  rownames(v) <- paste0("Cl.", 1:k)
  result <- list()
    result$v <- v
    result$ctype <- "obj"   
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}

uniquek <- function(x, k, sfidx){ 
  if(missing(x)) 
    stop("Missing data set")
  if(is.vector(x) || is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x)) 
    stop("Data set must be a numeric vector, data frame or matrix")
  n <- nrow(x); p <- ncol(x)
  if(missing(k))
    stop("Missing input argument k")
  if(!is.numeric(k))
    stop("Argument k must be a positive integer")
  if(k < 1 || k > n)
    stop(paste0("Argument k should be between 1 and ", n, " for the data set being processed."))
  if(!missing(sfidx)){
    if(!is.numeric(sfidx))
      stop(paste0("Selected feature index must be an integer between 1 and ", p, " for the data set being processed."))
    else
      sfidx <- round(sfidx)
    if(sfidx < 1 ||sfidx > p)
      stop(paste0("Selected feature index must be an integer between 1 and ", p, " for the data set being processed."))
  }else{
    sfidx <- 1
    nuniqx <- length(unique(x[,1]))
    if(p > 1){
      for(j in 2:p){
        if(length(unique(x[,j])) > nuniqx){
          nuniqx <- length(unique(x[,j]))
          sfidx <- j
        }
      }
    }
  }  
  v <- matrix(nrow=k, ncol=p, 0)
  colnames(v) <- colnames(x)
  rownames(v) <- paste0("Cl.", 1:k)
  xrows <- c()
  uniqx <- unique(x[,sfidx])
  if(k < 1 || k > length(uniqx))
    stop(paste0("k should be an integer between 1 and ", length(uniqx), " for the selected feature."))
  uniqx <- sample(uniqx, k, replace=FALSE)
  for(i in 1:k){
    ridx <- which(x[,sfidx] == uniqx[i])[1]
    xrows <- c(xrows, ridx)
  }
  v <- x[xrows,]
  if(k > 1)
    v <- as.matrix(x[xrows,])
  result <- list()
    result$v <- v
    result$ctype <- "obj"
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}

ursamp <- function(x, k){
  if(missing(x)) 
    stop("Missing data set")
  if(is.vector(x) || is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x)) 
    stop("Data set must be a numeric vector, data frame or matrix")
  n <- nrow(x); p <- ncol(x)
  if(missing(k))
    stop("Missing input argument k")
  if(!is.numeric(k))
    stop("Argument k must be a positive integer")
  if(k < 1 || k > n){
    errmsg <- paste0("Argument k should be between 1 and ", n, " for the data set being processed.")
    stop(errmsg)
  }
  v <- matrix(nrow=k, ncol=p, 0)
  for(i in 1:k) 
    for(j in 1:p)
      v[i,j] <- x[sample(n, 1),j]
  colnames(v) <- colnames(x)
  rownames(v) <- paste0("Cl.", 1:k)
  result <- list()
    result$v <- v
    result$ctype <- "obj"
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}

imembrand <- function(n, k, mtype, numseed){
  if(missing(n))
    stop("Missing number of rows argument for generating the membership matrix")
  if(!is.numeric(n))
    stop("Number of rows for the membership matrix must be a positive integer")
  else
    n <- as.integer(n)
  if(missing(k))
    stop("Missing input argument k")
  if(!is.numeric(k))
    stop("Argument k must be a positive integer")
  if(missing(mtype))
    mtype <- "f1"
  if(!is.element(mtype, c("f1","f2","f3","h")))
    stop(paste0("Invalid type ", mtype, ". The valid options are 'f1', 'f2' 'f3' and 'h'."))
  if(!missing(numseed)){
    if(!is.numeric(numseed))
      stop("Argument seed number should be a number")
    set.seed(numseed)
  }
  u <- matrix(nrow=n, ncol=k, 0)
  colnames(u) <- paste0("Cl.", 1:k)
  rownames(u) <- 1:n
  if(mtype=="f1"){
    for(i in 1:n){
      for(j in 1:k)
        u[i,j] <- sample(1:100, 1)
    }
    u <- u/apply(u, 1, sum)
    if(any(rowSums(u) != 1))
      u[,1] <- u[,1] + (1-rowSums(u)) 
  }
  if(mtype=="f2"){
    for(i in 1:n){
      remain = 100.0
      idx <- c(1:k)
      while(length(idx) >= 2){
        j <- sample(idx, 1)
        samp <- sample(0:remain, 1)
        u[i, j] <- samp
        remain <- remain - samp
        idx <- idx[!idx == j] 
      }
      u[i, idx[1]] <- remain
    }
    u <- u/apply(u, 1, sum)
    if(any(rowSums(u) != 1))
      u[,1] <- u[,1] + (1-rowSums(u)) 
  }
  if(mtype=="f3"){
    for(i in 1:n){
      remain <- 100.0 
      for(j in 1:k) {
        if(remain <= 100.0){
          samp <- sample(n, 1)
          if(samp <= remain){
            if(j==k)
              samp <- remain
          }else{
            if(j==k) 
              samp <- remain
            else 
              samp <- remain/j
          }
        }else{
           samp <- 0.0
        }
        u[i,j] <- samp 
        remain <- remain - samp
      }
    }
    u <- u/apply(u, 1, sum)
    if(any(rowSums(u) != 1))
      u[,1] <- u[,1] + (1-rowSums(u)) 
  }
  if(mtype=="h"){
    for(i in 1:n){
      if(i <= k)
        j <- i
      else
        j <- sample(k, 1)
      u[i, j] <- 1
    }
  }
  result <- list()
    result$u <- u
    if(mtype=="h")
      result$mtype <- "hard"
    else 
      result$mtype <- "fuzzy"
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}

imembones <- function(n, k, mtype="hrc", numseed){
  if(missing(k))
    stop("Missing input argument k")
  if(!is.numeric(k))
    stop("Argument k must be a positive integer")
  if(k < 1 || k > n)
    stop(paste0("k should be an integer between 1 and ", n, "."))
  if(!missing(numseed)){
    if(!is.numeric(numseed))
      stop("Argument seed number should be a number")
    set.seed(numseed)
  }
  if(!is.element(mtype, c("hfc","hlc", "hrc")))
    stop(paste0("Invalid initialization type ", mtype, ". Available options are 'hfc', 'hlc' and 'hrc'"))
  u <- matrix(nrow=n, ncol=k, 0)
  colnames(u) <- paste0("Cl.", 1:k)
  rownames(u) <- 1:n
  if(mtype=="hfc")
    u[,1] <- 1.0
  if(mtype=="hlc")
    u[,k] <- 1.0
  if(mtype=="hrc"){
    cidx <- runif(1,1,k)
    u[,cidx] <- 1.0
  }
  result <- list()
    result$u <- u
    result$mtype <- "hard"
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}
figen <- function(x, k, mtype="f", sfidx){
  if(missing(x)){
    stop("Missing data set")
  } 
  else{
    if(is(x,"data.frame"))
      x <- as.matrix(x)
    if(!is(x, "matrix"))
      stop("Data set must be a numeric data frame or matrix")
  }
  n <- nrow(x); p <- ncol(x)
  if(missing(k))
    stop("Missing input argument k")
  if(!is.numeric(k))
    stop("Argument k must be a positive integer")
  if(k < 1 || k > n)
    stop(paste0("k should be an integer between 1 and ", n, " for the data set being processed."))
  if(!mtype %in% c("f", "h")) 
    stop("Invalid membership type. The valid membership types are 'f' anf 'h'.")
  if(!missing(sfidx)){
    if(!is.numeric(sfidx))
      stop(paste0("Selected feature index must be an integer between 1 and ", p, " for the data set being processed."))
    else 
      sfidx <- round(sfidx)
    if(sfidx < 1 |sfidx > p)
      stop(paste0("Selected feature index must be an integer between 1 and ", p, " for the data set being processed."))
  }else{
    if(mean(x[, 1])==0) 
      mx <- 1
    else
      mx <- mean(x[, 1])
    cvx <- sqrt(var(x[, 1]))/mx
    for(j in 2:p){
      if(mean(x[, j])==0) 
        mx <- 1
      else
        mx <- mean(x[, j])
      if(sqrt(var(x[,j]))/mx > cvx){
        cvx <- sqrt(var(x[,j]))/mx
        sfidx <- j
      }
    }
  }
  u <- matrix(nrow=n, ncol=k, 0)
  colnames(u) <- paste0("Cl.", 1:k)
  rownames(u) <- 1:n
  xrange <- diff(range(x[,sfidx]))
  clbreaks <- seq(min(x[,sfidx]),max(x[,sfidx]), xrange/k)
  clmids <- c()
  for(i in 1:length(clbreaks)-1)
    clmids <- c(clmids, (clbreaks[i]+clbreaks[i+1])/2)
  if(mtype=="f"){
    for(i in 1:n){
      ridx <- which.min(abs(x[i,sfidx] - clmids))
      u[i,ridx] <-  xrange
      xs <- 1:k; xs <- xs[-ridx]
      for(j in xs)
        u[i,j] <- xrange-abs(j-ridx)*xrange/k
    }
    u <- u/rowSums(u)
  }else{
    for(i in 1:n){
      ridx <- which.min(abs(x[i, sfidx] - clmids))
      u[i, ridx] <- 1
    }
  }
  u <- u/apply(u, 1, sum)
  if(any(rowSums(u) != 1))
    u[,1] <- u[,1] + (1-rowSums(u)) 
  result <- list()
    result$u <- u
    if(mtype=="h")
      result$mtype <- "hard"
    else 
      result$mtype <- "fuzzy"
    result$call <- match.call()
  class(result) <- c("inaparc")
  return(result)
}

is.inaparc <- function(x){
  return(class(x)=="inaparc")
}

get.algorithms <- function(atype="prototype"){
  atype <- match.arg(atype, c("prototype", "membership"))
  if(atype=="prototype"){
    algonames <- vector(mode = "character", length = 27)
    algonames <- c("aldaoud", "ballhall", "crsamp", "firstk", "forgy", "hartiganwong", 
                   "inofrep", "inscsf", "insdev", "kkz", "kmpp", "ksegments", "ksteps",
                   "lastk", "lhsmaximin", "lhsrandom", "maximin", "mscseek", "rsamp",
                   "rsegment", "scseek", "scseek2", "spaeth", "ssamp", "topbottom", 
                   "uniquek", "ursamp")
  }
  if(atype=="membership"){
    algonames <- vector(mode = "character", length = 3)
    algonames <- c("imembones", "imembrand", "figen")
  }
  return(algonames)
}

.sqeucdist <- function(a, b){
  if(missing(a) || missing(b))
    stop("Missing argument to compute the distance")
  if(!is.numeric(a) || !is.numeric(b))
    stop("Input arguments must be numeric to compute the distance")
  return(t(a-b)%*%(a-b))
}

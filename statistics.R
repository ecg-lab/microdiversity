# Statistical framework for comparing empirical prokaryotic gene families to neutral simulations
# Tim Straub
# Zhaxybayeva Lab
# Department of Biological Sciences
# Dartmouth College

# requires zoo package
library(zoo)

# define user functions

# read in cluster pattern file from mothur
read.cluster <- function(file){
  read.delim(file, row.names=1)
}

# get all possible distances from cd plot
.get.dists <- function(cdplot){
  sapply(row.names(cdplot), as.numeric)
}

# get bins for cd plot
.get.bins <- function(dists, precision=.01){
  bins <- hist(dists, breaks=seq(0,ceiling(max(dists)*100)/100,precision), plot=F)$breaks
}

# find rows for bins
.find.row <- function(names, bin){
  for(name in names){
    if(as.numeric(name) >= bin){
      return(name)
    }
  }
  return(name)
}

# bin cdplot values
bin.cdplot.values <- function(data, dists, precision=.01){
  bins <- .get.bins(dists, precision)
  data <- na.locf(data)
  names <- row.names(data)
  keep.rows <- NULL
  i <- 1
  for(x in bins){
    keep.rows[i] <- .find.row(names, x)
    i <- i + 1
  }
  results <- t(data[keep.rows,])
  colnames(results) <- bins
  results
}

# calculate distances for each gene family from null distribution for multiple ranges
get.max.distance.vector <- function(data, null, bins){
  d <- matrix(nrow=nrow(data), ncol=length(bins))
  if(missing(null)){
    null <- data
  }
  med <- apply(null, 2, median)
  medabsdev <- apply(null, 2, mad)
  dists <- sapply(colnames(data), as.numeric)
  for(i in 1:nrow(data)){
    temp <- NULL
    for(j in 1:length(data[i,])){
      if(medabsdev[j] > 1){
        temp[j] <- (data[i,j] - med[j]) / medabsdev[j]
      } else {
        temp[j] <- data[i,j] - med[j]
      }
    }
    temp2 <- NULL
    for(k in 1:length(bins)){
      if(k == 1){
        temp3 <- temp[dists < bins[k] & dists > 0]
      } else {
        temp3 <- temp[dists < bins[k] & dists >= bins[k-1]]
      }
      maxtemp <- max(temp3, na.rm=T)
      mintemp <- min(temp3, na.rm=T)
      if(abs(mintemp) > maxtemp){
        temp2[k] <- mintemp
      } else {
        temp2[k] <- maxtemp
      }
    }
    d[i,] <- temp2
  }
  rownames(d) <- rownames(data)
  colnames(d) <- sapply(bins, as.character)
  d
}

# plot cd plot
plot.cdplot <- function(data, title='', xmax=-1, ymax=-1, log='y', col='black', axes=T, xlab='Divergence', ylab='# of Clusters', ...){
  
  data_locf <- na.locf(data)
  
  if(typeof(data)== 'list'){
    if(xmax < 0){
      xmax <- max(as.numeric(rownames(data)))
    }
  } else if(typeof(data) == 'integer'){
    if(xmax < 0){
      xmax <- max(as.numeric(names(data)))
    }
  }
  if(ymax < 0){
    ymax <- max(data, na.rm=T)
  }
  
  if(!axes){
    xlab=''
    ylab=''
    title=''
  } else if(missing(xlab) ) {
    xlab='Divergence'
  } else if(missing(ylab)) {
    ylab='# of Clusters'
  }
  
  if(length(col) == 1){
    coltemp = col
  }
  
  plot(x=1, y=1, type="n", main=title, xlim=c(xmax, 0), ylim=c(1, ymax), xlab=xlab, ylab=ylab, log=log, axes=axes)
  if(typeof(data)== 'list'){
    for(i in 1:ncol(data)){
      if(length(col) >= ncol(data)){
        coltemp = col[i]
      } else if(length(col) > 1){
        coltemp = col[ i %% length(col)]
      }
      temp <- data[,i]
      names(temp) <- rownames(data)
      temp <- na.omit(temp)
      par(new=T)
      plot(names(temp), temp, type="s", xlim=c(xmax, 0), ylim=c(1, ymax), main="", xlab="", ylab="", axes=F, log=log, col=coltemp, ...)
    }
  } else if(typeof(data) == 'integer') {
    temp <- na.omit(data)
    par(new=T)
    plot(names(temp), temp, type="s", xlim=c(xmax, 0), ylim=c(1, ymax), main="", xlab="", ylab="", axes=F, log=log, col=coltemp, ...)
  }
}

# examples

# read in ltt
example.dat <- na.locf(read.cluster('./example/dat.txt'))
example.sim <- na.locf(read.cluster('./example/sim.txt'))

# get divergences
example.dist <- c(.get.dists(example.dat), .get.dists(example.sim))

# bin data for use with get.max.distance.vector
example.dat.cd <- bin.cdplot.values(example.dat, example.dist)
example.sim.cd <- bin.cdplot.values(example.sim, example.dist)

# get distances from null from 0 to 0.05 divergence
example.d <- get.max.distance.vector(example.dat.cd, example.sim.cd, seq(0.05, 0.1, 0.05))
# get null distances only
example.d.null <- get.max.distance.vector(example.sim.cd, bins = seq(0.05, 0.1, 0.05))

# above simulations
example.above <- example.dat[ , example.d[,1] > quantile(example.d.null[,1], 0.975) ]

# below simulations
example.below <- example.dat[ , example.d[,1] < quantile(example.d.null[,1], 0.025) ]

# indistinguishable from simulations
example.indistinguishable <- example.dat[ , example.d[,1] <= quantile(example.d.null[,1], 0.975) & example.d[,1] >= quantile(example.d.null[,1], 0.025) ]


# plot some cd plots
pdf('./example/example.pdf', width=11, height=8.5)
par(mfrow=c(2,2))
plot.cdplot(example.sim, 'Simulations', xmax = 0.3, ymax = 200)

plot.cdplot(example.sim, 'Indistinguishable', xmax = 0.3, ymax = 200)
par(new=T)
plot.cdplot(example.indistinguishable, title = '', axes = F, xlab = "", ylab = "", col = 'purple', xmax = 0.3, ymax = 200)

plot.cdplot(example.sim, 'Below', xmax = 0.3, ymax = 200)
par(new=T)
plot.cdplot(example.below, title = '', axes = F, xlab = '', ylab = '', col = 'steelblue1', xmax = 0.3, ymax = 200)

plot.cdplot(example.sim, 'Above', xmax = 0.3, ymax = 200)
par(new=T)
plot.cdplot(example.above, title = '', axes = F, xlab = '', ylab = '', col = 'orange', xmax = 0.3, ymax = 200)
dev.off()
##############################################################################################
getDistanceFromLatLonInKm = function(lat1,lon1,lat2,lon2) {
  ##############################################################################################
  R = 6371;                      # Radius of the earth in km
  dLat = (lat2-lat1)* (pi/180);  # deg2rad below
  dLon = (lon2-lon1)* (pi/180);
  a =
    sin(dLat/2) * sin(dLat/2) +
    cos((lat1)* (pi/180)) * cos((lat2)* (pi/180)) *
    sin(dLon/2) * sin(dLon/2)

  c = 2 * atan2(sqrt(a), sqrt(1-a));
  d = R * c; # Distance in km
  return(d)
}

##############################################################################################
conversion <- function(albers,coord){
  ##############################################################################################
  dmat1 <- sqrt(outer(albers[,1], albers[,1], "-")^2 + outer(albers[,2], albers[,2], "-")^2)
  dmat2 <- matrix(0,dim(coord)[1],dim(coord)[1])
  for(i in 1:dim(coord)[1]){
    for(j in 1:i){
      dmat2[i,j] <- getDistanceFromLatLonInKm(coord[i,2],coord[i,1],coord[j,2],coord[j,1])
    }
  }
  dmat2 <- t(dmat2) + dmat2
  result <- mean(dmat2/dmat1,na.rm=T)
  return(result)
}


#' Estimate local gradients.
#'
#' @param dates vector of length \eqn{N} representing the observed time when the process first appears.
#' @param covar matrix of dimension \eqn{N} x \eqn{p} representing explanatory variables except for coordinates. If not specified, only \code{coord.longlat} is used as an explanatroy variable.
#' @param coord.longlat matrix of dimension \eqn{N} x 2 representing the longitude and lattitude of each observation.
#' @param newcoord.longlat matrix of dimension \eqn{n} x 2 representing the longitude and lattitude locations to estimate the gradient. The default is \code{coord.longlat}.
#' @param Albers vector of length \eqn{2} representing the parameters of Albers equal area conic projection. If process of interest is on the longitude and lattitude domain, user must specify \code{Albers}. If not specified, projection is not conducted and mean speed of spread is not evaluated in the \code{\link{summary.localgrad}}. See \code{\link{mapproject}} for an explanation.
#' @param n.samp total number of posterior samples. Defaults to 10,000.
#' @param thin the thinning interval between consecutive observations. Defaults to 10.
#' @param knots if specified, use kernel convolution method over knot process to speed up computation. See \code{\link{spLM}} for an explanation.
#' @param amcmc if specified, use adaptive mcmc algorithm. See \code{\link{spLM}} for an explanation.
#' @param verbose if TRUE, progress of the sampler is printed to the screen. Otherwise, nothing is printed.
#' @return \code{localgrad} returns an object of class \dQuote{\code{localgrad}}, which is a list containing the following components.
#' \item{dates}{ vector of length \eqn{N} representing the observed time when the process first appears. }
#' \item{coord.longlat}{ matrix of dimension \eqn{N} x 2 representing the longitude and lattitude of each observation. }
#' \item{newcoord.longlat}{ matrix of dimension \eqn{n} x 2 representing locations to estimate the gradient. }
#' \item{coord }{ matrix of dimension \eqn{N} x 2 representing the Albers conic projected observation. If projection is not used, \code{coord} is identical to \code{coord.longlat}. }
#' \item{newcoord }{ matrix of dimension \eqn{n} x 2 representing the Albers conic projected locations to be estimated. If projection is not used, \code{newcoord} is identical to \code{newcoord.longlat}. }
#' \item{conv.factor}{ conversion factor to minimize distortion between Albers conic projected coordinates and longitude latitude coordinates. If projection is not used, \code{conv.factor} is 1.}
#' \item{samples }{ list of posterior samples for the gradient at each location in newcoord. Each posterior sample consists of an \eqn{x} and \eqn{y} component of the gradient.}
#' \item{means}{ mean speed of spread at each location in newcoord. }
#' \item{sig}{ vector of representing whether the gradients are significant or not in the \code{newcoord}. }
#' \item{post.samp}{ thined posterior samples for \eqn{\Theta} from \code{\link{spLM}} function. }
#' @references
#' Banerjee, S., and Gelfand, A. E. (2006). Bayesian wombling: Curvilinear gradient assessment under spatial process models. \emph{Journal of the American Statistical Association}, \bold{101}(476), 1487-1501.
#' @references
#' Banerjee, S., Gelfand, A. E., and Sirmans, C. F. (2003). Directional rates of change under spatial process models. \emph{Journal of the American Statistical Association}, \bold{98}(464), 946--954.
#' @references
#' Goldstein, J., Jaewoo Park, Haran, M., Liebhold, A., and Bjornstad, O. N. (2018). Quantifying Spatio-Temporal Variation of Invasion Spread.\emph{arXiv preprint arXiv:1506.02685v3}
#' @details
#' This function estimates gradients of surface by specifying a Gaussian process model as in Goldstein et al. (2018) using theoretical results in Banerjee et al. (2003); Banerjee and Gelfand (2006).
#' If the surface is the waiting time to first appearance of the process, then the reciprocal of the gradient length is a measure of the invasion speed: Fast spread leads to shallow waiting time surfaces, while slow spread results in steep surfaces.
#'
#' Consider \eqn{N} observations of the year of first appearance \eqn{Y(s)} at location \eqn{s}. \eqn{Y(s)} can be modeled using an isotropic
#' Gaussian process with mean \eqn{\mu(s)} and covariance \eqn{K(D)}, where \eqn{D} is the \eqn{N} x \eqn{N} matrix of pairwise Euclidian distances of \eqn{s}.
#' Then the gradient of \eqn{Y(s)} at \eqn{s}, the vector of directional derivatives can be defined as \eqn{dY(s)}.
#' The conditional distribution of the gradient can be derived as:
#'
#' \deqn{dY(s)\mid Y,\Theta \sim N(d\mu(s)- dK(D)^\prime K(D)^{-1}(Y-\mu), -H(D)-dK(D)^\prime K(D)^{-1}dK(D) ),}{dY(s)|Y,\Theta ~ N(d\mu(s)- dK(D)'K(D)^{-1}(Y-\mu), -H(D)-dK(D)'K(D)^{-1}dK(D) ),}
#'
#' where \eqn{\Theta = (\beta^\prime, \eta)^\prime}{\Theta = (\beta', \eta)'} is the parameters for the Gaussian process;
#' \eqn{\beta} is the regression parameter for \eqn{\mu(s)} and \eqn{\eta} is the covariance parameter in \eqn{K(D)} and \eqn{H(D)} is the Hessian matrix of \eqn{K(D)}.
#'
#' \code{localgrad} function conducts Bayesian inference for \eqn{\Theta} through the \code{spLM} function.
#' Simulation from the predictive distribution of the gradient is done by composition; given each posterior sample of
#' \eqn{\Theta}, a sample for \eqn{dY(s)} is drawn from the conditional distribution of the gradient.
#' Then speed of spread and the dominant direction of the spread are calculated and tested.
#'
#'
#'@seealso \code{\link{spLM}}
#'@export
#'@importFrom mapproj mapproject
#'@importFrom geoR variofit
#'@importFrom geoR variog
#'@importFrom spBayes spLM
#'@importFrom spBayes spRecover
#'@import stats
#'@examples \dontrun{
#'
#' ###########    hemlock data example    ##########
#' data(hemlock)
#' coord = cbind(hemlock$long,hemlock$lat)    # longitude latitude coordinates
#' dates = hemlock$first.year                 # quaratine data by county
#'
#' # adaptive mcmc algorithm
#' # specify batches, batch length, target acceptance rate
#' amcmc=list(n.batch=10,batch.length=100,accept.rate=0.3)
#' out.grad = localgrad(dates=dates,coord.longlat=coord, Albers=c(29.5,45.5),
#' n.samp=1000,thin=10, amcmc=amcmc)
#' }
##############################################################################################
localgrad <- function(dates, covar = NULL, coord.longlat, newcoord.longlat = NULL, Albers = NULL, n.samp = 10000, thin = 10, knots = NA, amcmc = NA, verbose = FALSE){
  ##############################################################################################

  # without Albers equal area conic projection
  if(is.null(Albers)){
    coord = coord.longlat
    conv.factor = 1
    # Default to predicting gradient at observed locations
    if(is.null(newcoord.longlat)){
      newcoord.longlat <- coord.longlat
    }
    newcoord = newcoord.longlat

    # Albers equal area conic projection
  }else{
    proj = mapproject(coord.longlat[,1], coord.longlat[,2], projection="albers", parameters=Albers, orientation=NULL)
    coord = cbind(proj$x,proj$y)
    # scale to km (minimize distortion)
    conv.factor = conversion(coord,coord.longlat)

    # Default to predicting gradient at observed locations
    if(is.null(newcoord.longlat)){
      newcoord.longlat <- coord.longlat
    }

    # Albers equal area conic projection for newcoordinate
    proj.n = mapproject(newcoord.longlat[,1], newcoord.longlat[,2], projection="albers", parameters=Albers, orientation=NULL)
    newcoord = cbind(proj.n$x,proj.n$y)
  }

  # Construct covariates
  n.obs <- length(dates)
  n.pred <- dim(newcoord)[1]
  covariates = cbind(coord,covar)
  n.beta = dim(covariates)[2]+1
  formula = dates~covariates

  # matrix of pairwise distances
  dmat <- sqrt(outer(coord[,1], coord[,1], "-")^2 + outer(coord[,2], coord[,2], "-")^2)

  # get prior estimates of covariance parameters from variogram (trend is only using coordinates)
  vario.out <- variofit(variog(coords = coord, data=dates, trend = ~ coord[,1] + coord[,2], messages=verbose), fix.kappa = TRUE, kappa = 1.5, weights="cressie",max.dist=as.numeric(quantile(dmat[dmat!=0],.5)), messages=verbose)
  tausq.vario <- vario.out$nugget
  sigmasq.vario <- vario.out$cov.pars[1]
  phi.vario <- 1/vario.out$cov.pars[2] # transform to agree w/ spBayes parameterization

  # provision for bad variogram fit
  if(phi.vario == 0) { phi.vario = min(dmat[dmat!=0]) }
  if(sigmasq.vario == 0) { sigmasq.vario = min(dmat[dmat!=0]) }
  if(tausq.vario == 0) { tausq.vario = sigmasq.vario/100 }

  # fit a Gaussian process, Matern nu=3/2 using spBayes package
  # starting value for the regression parameter
  start.beta <- c(mean(dates),rep(1,n.beta-1))
  priors <- list("phi.Unif"=c(1/max(dmat),1/min(dmat[dmat!=0])), "sigma.sq.IG"=c(2, sigmasq.vario),
                 "tau.sq.IG"=c(2, tausq.vario), "nu.unif"=c(1,2))
  starting <- c(start.beta,"phi"=phi.vario, "sigma.sq"=sigmasq.vario, "tau.sq"=tausq.vario, "nu"=1.5)
  tuning <- lapply(starting,function(x){x/100})
  tuning$nu <- 0 # fix smoothness parameter at 1.5

  if(is.na(knots) && !exists('amcmc')) {
    krige.out <- spLM(formula, coords=coord, starting=starting, tuning=tuning, priors=priors,
                      cov.model = "matern",  n.samples = n.samp, n.report = 500, verbose = verbose)
  } else if(!is.na(knots) && !exists('amcmc')) {
    krige.out <- spLM(formula, coords=coord, starting=starting, tuning=tuning, priors=priors,
                      cov.model = "matern",  n.samples = n.samp, n.report = 500, verbose = verbose, knots=knots)
  } else if(is.na(knots) && exists('amcmc')) {
    krige.out <- spLM(formula, coords=coord, starting=starting, tuning=tuning, priors=priors,
                      cov.model = "matern",  n.samples = n.samp, n.report = 500, verbose = verbose, amcmc=amcmc)
  } else {
    krige.out <- spLM(formula, coords=coord, starting=starting, tuning=tuning, priors=priors,
                      cov.model = "matern",  n.samples = n.samp, n.report = 500, verbose = verbose, knots=knots, amcmc=amcmc)
  }
  burn.in <- 0.5*n.samp
  # recover beta and spatial random effects
  krige.samp <- spRecover(krige.out, start=burn.in, verbose=verbose)
  posterior.samples <- cbind(krige.samp$p.beta.recover.samples,krige.samp$p.theta.recover.samples)


  # thin posterior samples to numsamp
  numsamp <- n.samp/thin
  thin.ind <- seq(1,dim(posterior.samples)[1],length=floor(numsamp))
  posterior.samples <- posterior.samples[thin.ind,]

  # draw posterior gradients from each location in newcoord following (Banerjee 2003)
  posterior.out <- list()
  for(i in 1:n.pred) {
    posterior.out[[i]] <- matrix(nrow=0,ncol=2)
  }
  if(verbose==TRUE){cat("Computing gradients for",numsamp,"samples:\n")}

  if(dim(covariates)[1]!=0){ covar.mat = as.matrix(cbind(1,covariates)) }else { covar.mat = 1 }

  for(j in 1:dim(posterior.samples)[1]) {
    beta <- posterior.samples[j,1:(n.beta)]
    sigmasq <- posterior.samples[j,n.beta+1]
    tausq <- posterior.samples[j,n.beta+2]
    phi <- posterior.samples[j,n.beta+3]

    Kinv <- solve( sigmasq * (1+phi*dmat) * exp(-phi*dmat) + tausq*diag(n.obs) )
    mu <- covar.mat %*% beta
    grad.mu <- beta[2:3]
    grad.mu = matrix(rep(grad.mu,n.pred),nrow=n.pred,byrow=TRUE)

    # now get 10 gradient samples at each location (total samples 10*numsamp)
    for(i in 1:n.pred) {
      s0 <- newcoord[i,]
      delta <- cbind(s0[1] - coord[,1], s0[2] - coord[,2])
      delta.mag <- sqrt( (s0[2]-coord[,2])^2 + (s0[1]-coord[,1])^2 )

      gam <- -sigmasq*phi^2*exp(-phi*delta.mag)*delta
      mean.grad <- grad.mu[i,] - as.vector( t(gam)%*%Kinv%*%(dates-mu) )
      var.grad <- sigmasq*phi^2*diag(2) - t(gam)%*%Kinv%*%gam

      # multivariate normal sample of size 10
      z <- matrix(rnorm(10*2),10,2) %*% chol(var.grad)
      y <- t(mean.grad + t(z))
      posterior.out[[i]] <- rbind(posterior.out[[i]], y)
    }
    if(verbose==TRUE){
      cat(10*j,",",sep="")
      if(j%%20==0) cat("\n")
    }
  }

  # test for significant gradients
  normal.density <- rep(NA,n.pred)
  for(i in 1:n.pred) {
    pts <- posterior.out[[i]]
    rho <- cor(pts[,1],pts[,2])
    mu.x <- mean(pts[,1]); mu.y <- mean(pts[,2])
    sd.x <- sd(pts[,1]); sd.y <- sd(pts[,2])
    normal.density[i] <- 1/(1-rho^2) * ( mu.x^2/sd.x^2 + mu.x^2/sd.x^2 - 2*rho*mu.x*mu.y/(sd.x*sd.y) )
  }
  is.sig <- normal.density > pchisq(.95,2)

  # transform gradients into speeds of spread (polar transformation)
  spread.means <- t(sapply( posterior.out, function(x){apply(x,2,mean)} ))
  x <- spread.means[,1]
  y <- spread.means[,2]
  r <- sqrt(x^2+y^2)
  theta <- atan(y/x) + (x<=0)*pi
  r <- 1/r
  theta <- theta + pi
  speed.means <- cbind(r*cos(theta),r*sin(theta))

  res <- list(call=deparse(match.call(),width.cutoff = 500), dates = dates, coord.longlat = coord.longlat, newcoord.longlat = newcoord.longlat,
              coord = coord, newcoord = newcoord, conv.factor=conv.factor, samples = posterior.out, means = speed.means, sig = is.sig, post.samp = posterior.samples)
  class(res) <- "localgrad"
  res
}





#' Draw the vector field plot for the significant gradient from the class \dQuote{\code{localgrad}}.
#'
#' @param object an object of class \code{localgrad}, typically the result of a call to \code{\link{localgrad}}.
#' @param main an overall title for the plot: see \code{\link{plot}}. The default title is "Speed of spread".
#' @param scale scale of arrow vector. The default is 1.
#' @param database character string naming a geographical database in \code{\link{map}} function. If specified, draw geographical map on the gradient plot. The default is not using \code{\link{map}} function.
#' @param xlab a title for the \eqn{x} axis: see \code{\link{plot}}. The default title is "Longitude".
#' @param ylab a title for the \eqn{y} axis: see \code{\link{plot}}. The default title is "Latitude".
#' @param xlim numeric vectors of length 2, giving the ranges for \eqn{x} axis: see \code{\link{plot}}. The default is range of \eqn{x} axis in the predicted coordinates.
#' @param ylim numeric vectors of length 2, giving the ranges for \eqn{y} axis: see \code{\link{plot}}. The default is range of \eqn{y} axis in the predicted coordinates.
#' @param \dots additional arguments.
#' @details
#' This function draws the vector field plot for class \dQuote{\code{localgrad}}. Arrows in the plot indicate the direction and speed of spread.
#' As the speed of spread increases, the length of the arrows become longer.
#' The color of each arrow represents the first appearance time. Red indicates the earliest appearance, and blue represents the latest appearance.
#'
#' @seealso \code{\link{localgrad}}, \code{\link{map}}
#' @export
#' @import maps
#' @import grDevices
#' @import graphics
#' @examples \dontrun{
#'
#' ###########    hemlock data example    ##########
#' data(hemlock)
#' coord = cbind(hemlock$long,hemlock$lat)    # longitude latitude coordinates
#' dates = hemlock$first.year                 # quaratine data by county
#'
#' # adaptive mcmc algorithm
#' # specify batches, batch length, target acceptance rate
#' amcmc=list(n.batch=10,batch.length=100,accept.rate=0.3)
#' out.grad = localgrad(dates=dates,coord.longlat=coord, Albers=c(29.5,45.5),
#' n.samp=1000,thin=10, amcmc=amcmc)
#'
#' # gradient plot
#' plotgrad(out.grad,cex=1,pch=".",database="state")
#' }
##############################################################################################
plotgrad = function(object, main="Speed of spread",scale=1,database=NULL,
                    xlab=expression(paste(degree,"Longitude")),ylab=expression(paste(degree,"Latitude")),xlim=NA, ylim=NA,...){
  ##############################################################################################
  dates <- object$dates
  means <- object$means
  sig <- object$sig
  plotcoord <- object$newcoord.longlat
  plotcoord <- plotcoord[sig,]
  plotmeans <- means[sig,]

  if(is.na(xlim[1])) { xlim=range(plotcoord[,1]) }
  if(is.na(ylim[1])) { ylim=range(plotcoord[,2]) }

  colramp = colorRampPalette(c("red", "blue"))(10)
  # divide into 10 equal intervals, assign a color to each
  vect = c(seq(min(dates),max(dates),length=12)[2:11],max(dates)+1)
  plotcol = colramp[findInterval(dates,vect,all.inside=TRUE)]

  plot(plotcoord,xlab=xlab,ylab=ylab,main=main,xlim=xlim,ylim=ylim,col=plotcol,...)
  # goal: plot each arrow by hand
  # goal: color code object by date (optionally)
  #maxlength=4 # length of longest arrow
  #scale=0.1 # scaling factor for arrow
  r = sqrt(plotmeans[,1]^2 + plotmeans[,2]^2)
  theta = atan2(plotmeans[,2],plotmeans[,1])
  #r.scale = r/max(r)*maxlength
  r.scale = r/max(r)*scale
  arrows( x0 = plotcoord[,1], y0 = plotcoord[,2], x1 = plotcoord[,1] + r.scale*cos(theta), y1 = plotcoord[,2] + r.scale*sin(theta), length = 0.1, col=plotcol[sig] )

  if(!is.null(database)){ map(database=database, fill = FALSE, col="darkgray", add=TRUE,...) }
}



#' Print a summary of a class \dQuote{\code{localgrad}}.
#' @param object an object of class \code{localgrad}, typically the result of a call to \code{\link{localgrad}}.
#' @param \dots additional arguments.
#' @details
#' This is the summary function for class \dQuote{\code{localgrad}}. Among the predicted coordinates this function displays the number of significant locations.
#' If the process of interest is on the longitude and lattitude domian, this also returns the estimated mean and median speed.
#'@seealso \code{\link{localgrad}}
#'@method summary localgrad
#'@export
#' @examples \dontrun{
#'
#' ###########    hemlock data example    ##########
#' data(hemlock)
#' coord = cbind(hemlock$long,hemlock$lat)    # longitude latitude coordinates
#' dates = hemlock$first.year                 # quaratine data by county
#'
#' # adaptive mcmc algorithm
#' # specify batches, batch length, target acceptance rate
#' amcmc=list(n.batch=10,batch.length=100,accept.rate=0.3)
#' out.grad = localgrad(dates=dates,coord.longlat=coord, Albers=c(29.5,45.5),
#' n.samp=1000,thin=10, amcmc=amcmc)
#'
#' # summary of inference results
#' summary(out.grad)
#' }
##############################################################################################
summary.localgrad<-function(object, ...){
  ##############################################################################################
  means.km <- object$conv.factor*object$means # scale to km
  mag.speed <- sqrt( means.km[,1]^2 + means.km[,2]^2 )

  mean.speed <- mean(mag.speed)
  median.speed <- median(mag.speed)

  n.pred <- length(object$samples)
  n.sig <- sum(object$sig)

  cat("\nCall:\n",object$call)
  cat("\n\n")
  cat("Speed of spread estimated at",n.pred,"locations.\n","Of these,",n.sig,"are significantly nonzero.\n")
  cat("\n\n")
  if(!object$conv.factor==1){
  cat("Mean speed of spread is estimated as",mean.speed,"km.\n","Median speed of spread is estimated as",median.speed,"km.\n")
  }
}


#' Test class \dQuote{\code{localgrad}} for potential regions of long-range jumps
#' using a Rayleigh test of a uniform circular distribution (Jammalamadaka and Sengupta, 2001).
#' @param object an object of class \code{localgrad}, typically the result of a call to \code{\link{localgrad}}.
#' @param r radius of jump to perform the Rayleigh test. The default is the 2.5 percent quantile of the distance matrix of Albers conic projected coordinates. If projection is not used in \code{\link{localgrad}}, the default is the 2.5 percent quantile of the distance matrix of original coordinates.
#' @param plot TRUE: plot points around which the spread is radially outward.
#' @param add whether we add to the original plot or not.
#' @param main an overall title for the plot: see \code{\link{plot}}. The default title is "Test for long-range jumps".
#' @param xlab a title for the \eqn{x} axis: see \code{\link{plot}}. The default title is "Longitude".
#' @param ylab a title for the \eqn{y} axis: see \code{\link{plot}}. The default title is "Latitude".
#' @param col color code or name, see \code{\link{par}}.
#' @param cex character (or symbol) expansion: a numerical vector. This works as a multiple of \code{\link{par}}("cex").
#' @param \dots additional arguments.
#' @return \code{longrange.localgrad} returns an object of class \dQuote{\code{longrange.localgrad}}, which is a list containing the following components.
#' \item{rayleigh.n}{ number of points in the neighborhood used for the Rayleigh test at each location. }
#' \item{rayleigh.z}{ summary statistic for the Rayleigh test at each location. }
#' \item{sig}{ binary vector, 1 for locations where the Rayleigh test is significant. }
#' @references
#' Goldstein, J., Jaewoo Park, Haran, M., Liebhold, A., and Bjornstad, O. N. (2018). Quantifying Spatio-Temporal Variation of Invasion Spread.\emph{arXiv preprint arXiv:1506.02685v3}
#' @references
#' Jammalamadaka, S. R. and Sengupta, A. (2001). Topics in circular statistics. \emph{(volume 5). World Scientific.}
#' @details
#' This function conducts the Rayleigh test to the gradients from the class \dQuote{\code{localgrad}}.
#' Then flag at the potential sites of long range jumps.
#' @seealso \code{\link{localgrad}}
#' @export
#' @import graphics
#' @examples \dontrun{
#'
#' ###########    hemlock data example    ##########
#' data(hemlock)
#' coord = cbind(hemlock$long,hemlock$lat)    # longitude latitude coordinates
#' dates = hemlock$first.year                 # quaratine data by county
#'
#' # adaptive mcmc algorithm
#' # specify batches, batch length, target acceptance rate
#' amcmc=list(n.batch=10,batch.length=100,accept.rate=0.3)
#' out.grad = localgrad(dates=dates,coord.longlat=coord, Albers=c(29.5,45.5),
#' n.samp=1000,thin=10, amcmc=amcmc)
#'
#' # gradient plot
#' plotgrad(out.grad,cex=1,pch=".",database="state")
#'
#' # rayleigh test
#' out.longrange = longrange.localgrad(out.grad,add=TRUE)
#' }
##############################################################################################
longrange.localgrad <- function(object, r=NULL, plot=TRUE, add = FALSE, main = "Test for long range jumps",
                                xlab ="Longitude", ylab = "Latitude", col=1, cex=5, ...){
  ##############################################################################################
  conv.factor <- object$conv.factor
  plotcoord <- object$newcoord.longlat
  coord <- object$coord
  newcoord <- object$newcoord
  means <- object$means
  sig <- object$sig

  # consider only signficant gradients
  newcoord.test = newcoord[sig,]
  means.test = means[sig,]

  # normalized directions of spread
  direction.spread <- means.test / sqrt( means.test[,1]^2 + means.test[,2]^2)

  if(is.null(r)){
    dmat.coord <- sqrt(outer(coord[,1], coord[,1], "-")^2 + outer(coord[,2], coord[,2], "-")^2)
    r <- quantile(dmat.coord,.025)
  }

  dmat <- sqrt(outer(newcoord[,1], newcoord.test[,1], "-")^2 + outer(newcoord[,2], newcoord.test[,2], "-")^2)
  dkl <- ifelse(dmat > 0 & dmat < r, 1, NA)

  rayleigh.r <- rep(NA,dim(newcoord)[1])
  for(i in 1:dim(newcoord)[1]) {	# note: vectorize this
    coord.neigh <- direction.spread[!is.na(dkl[i,]),]
    if(length(coord.neigh) <= 2) rayleigh.r[i] = NA
    else {
      sums.neigh <- colSums(coord.neigh)
      rayleigh.r[i] <- sqrt( sums.neigh[1]^2 + sums.neigh[2]^2 ) / dim(coord.neigh)[1]
    }
  }
  rayleigh.n <- rowSums(dkl,na.rm=TRUE)
  rayleigh.z <- rayleigh.n * rayleigh.r^2
  sig2 <- rayleigh.z < 3 & rayleigh.n > 5

  if(plot) {
    if(!add) plot(plotcoord[sig2,],col=col,pch=".",cex=cex,xlim=range(plotcoord[,1]),ylim=range(plotcoord[,2]), xlab=xlab, ylab=ylab, main=main)
    else points(plotcoord[sig2,],col=col,pch=".",cex=cex)
  }
  res = list(rayleigh.n=rayleigh.n,rayleigh.z=rayleigh.z,coord.sig=plotcoord[sig2,])
  class(res) <- "longrange.localgrad"


  # print
  if(conv.factor==1){
    cat("Long range jumps tested for" ,r,".\n")
  }else{
    cat("Long range jumps tested for ",conv.factor*r,"km.\n") # scale to km #
  }
  return(res)
}


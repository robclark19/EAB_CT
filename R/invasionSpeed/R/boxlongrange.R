
##############################################################################################
gam.integrand = function(t,j,s0,u,ut,sigmasq,phi,coord) {
  ##############################################################################################
  delta = unlist(s0) - coord[j,] + u%*%t
  delta.mag = sqrt(delta[1,]^2 + delta[2,]^2)
  gradK = -sigmasq*phi^2*exp(-phi*delta.mag) * t(delta)
  gradK %*% ut
}

##############################################################################################
K.integrand = function(t,s0,u,ut,sigmasq,phi,coord) {
  ##############################################################################################

  if(diff(t) == 0){return(0)}
  delta1 = u[1] * (t[2] - t[1])
  delta2 = u[2] * (t[2] - t[1])
  delta.mag = sqrt(delta1^2 + delta2^2)
  #if(delta.mag == 0){ return(phi*ut[1]^2 + phi*ut[2]^2) } # limit of H as distance -> 0
  H11 = -sigmasq*phi^2*exp(-phi*delta.mag)*(1-phi*delta1^2/delta.mag)
  H22 = -sigmasq*phi^2*exp(-phi*delta.mag)*(1-phi*delta2^2/delta.mag)
  -(H11*ut[1]^2 + H22*ut[2]^2)
}



##############################################################################################
flux.line = function(s0,u,ut,tstar,jump.data) {
  ##############################################################################################
  with(jump.data, {
    beta = params[1:n.beta]
    sigmasq = params[n.beta+1]
    tausq = params[n.beta+2]
    phi = params[n.beta+3] # make sure parameterization agrees w/ Banerjee

    u = matrix(u); ut = matrix(ut);
    mu.flux = tstar*t(ut)%*%beta[2:3] 	# correct for beta0+beta1*x+beta2*y

    if(dim(covariates)[1]!=0){ covar.mat = as.matrix(cbind(1,covariates)) }
    else{ covar.mat = 1 }
    mu <- covar.mat %*% beta

    # integrate gamma for each j
    n = length(dates)
    gamma.flux = rep(NA,n)
    for(j in 1:n) {
      gamma.flux[j] = integrate(gam.integrand,lower=0,upper=tstar,j=j,s0=s0,u=u,ut=ut,
                                sigmasq=sigmasq,phi=phi,coord=coord)$value
    }
    # integrate to find K only once
    K.flux = adaptIntegrate(K.integrand, lowerLimit=c(0,0), upperLimit=c(tstar,tstar),
                            s0=s0,u=u,ut=ut,sigmasq=sigmasq,phi=phi,coord=coord)$integral
    # get mean and variance of the total flux across the line segment
    #Kinv = solve( sigmasq * (1+phi*dmat) * exp(-phi*dmat ) + tausq*sigmasq*diag(n) )
    Kinv = solve( sigmasq * (1+phi*dmat) * exp(-phi*dmat ) + tausq*diag(n) )

    mean.flux = mu.flux + as.vector( -t(gamma.flux)%*%Kinv%*%(dates - mu) )
    var.flux = K.flux - t(gamma.flux)%*%Kinv%*%gamma.flux
    return(c(mean.flux,var.flux))
  })
}


##############################################################################################
flux.box = function(center,jump.data) {
  ##############################################################################################

  with(jump.data, {
    flux.up = flux.line(s0 = c(center[1] - r, center[2] + r), u = c(1,0), ut = c(0,1), tstar = 2*r, jump.data)
    flux.right = flux.line(c(center[1] + r, center[2] + r), c(0,-1), c(1,0), 2*r, jump.data)
    flux.down = flux.line(c(center[1] + r, center[2] - r), c(-1,0), c(0,-1), 2*r, jump.data)
    flux.left = flux.line(c(center[1] - r, center[2] - r), c(0,1), c(-1,0), 2*r, jump.data)

    # get 10000 posterior samples of the average flux out of the box
    #samp.up = rnorm(10000,mean=flux.up[1],sd=sqrt(flux.up[2]))
    #samp.right = rnorm(10000,mean=flux.right[1],sd=sqrt(flux.right[2]))
    #samp.down = rnorm(10000,mean=flux.down[1],sd=sqrt(flux.down[2]))
    #samp.left = rnorm(10000,mean=flux.left[1],sd=sqrt(flux.left[2]))
    #samps = list(samp.up/(2*r),samp.down/(2*r),samp.left/(2*r),samp.right/(2*r))

    #samp.flux = (samp.up+samp.right+samp.down+samp.left)/(8*r)
    # use xtable to print a nice table showing the mean + stdev of the average flux along each side,
    # in addition to the overall average

    #meangrad = sapply(samps,mean)
    #sdgrad = sapply(samps,sd)

    meangrad = c(flux.up[1]/(2*r),flux.down[1]/(2*r),flux.left[1]/(2*r),flux.right[1]/(2*r))
    sdgrad = c(sqrt(flux.up[2])/(2*r),sqrt(flux.down[2])/(2*r),sqrt(flux.left[2])/(2*r),sqrt(flux.right[2])/(2*r))
    #sides = c("Top","Bottom","Left","Right","Total")
    #out = data.frame(meangrad,sdgrad)
    #rownames(out) = sides
    #colnames(out) = c("Avg Gradient","StDev")
    #print(xtable(out))

    # return whether each side has a significantly negative gradient (-1), significantly positive gradient (+1) or non-significant gradient (0)
    sig.sides = c(0,0,0,0)
    sig.sides = sig.sides - ((meangrad + 1.96*sdgrad) < 0) # significantly small gradient (long range jump)
    sig.sides = sig.sides + ((meangrad - 1.96*sdgrad) > 0) # significantly large gradient
    return( sig.sides )
  })
}



#' Scan over predicted coordinates for areas of long-range jumps.
#' @param dates vector of length \eqn{N} representing the observed time when the process first appears.
#' @param coord.longlat matrix of dimension \eqn{N} x 2 representing the longitude and lattitude of each observation.
#' @param newcoord.longlat matrix of dimension \eqn{n} x 2 representing the longitude and lattitude locations to estimate the jumps. The default is \code{coord.longlat}.
#' @param r radius about the jumps. The default is the 2.5 percent quantile of the distance matrix of coordinates.
#' @param params estiamted parameter \eqn{\Theta = (\beta^\prime, \eta)^\prime}{\Theta = (\beta', \eta)'}  for the Gaussian process model using \code{\link{localgrad}}. \eqn{\beta} represents the regression coefficients and \eqn{\eta} represents the covariance parameters.
#' @param Albers vector of length \eqn{2} representing the parameters of Albers equal area conic projection. If process of interest is on the longitude and lattitude domain, user must specify \code{Albers}. If specified, use the Albers equal area conic projection. The default setting does not use projections. See \code{\link{mapproject}} for an explanation.
#' @param covar matrix of dimension \eqn{N} x \eqn{p} of explanatory variables except for coordinates.
#' @return \code{jump.scan} returns an object of class \dQuote{\code{jump.scan}}, which is a list containing the following components.
#' \item{coord.longlat}{ matrix of dimension \eqn{N} x 2 representing the longitude and lattitude of each observation. }
#' \item{newcoord.longlat}{ matrix of dimension \eqn{n} x 2 representing longitude and lattitude locations to estimate the jumps.}
#' \item{conv.factor}{ conversion factor to minimize distortion between Albers conic projected coordinates and longitude latitude coordinates. If projection is not used, the \code{conv.factor} is 1.}
#' \item{r}{ radius about jumps.}
#' \item{all.sides}{ Matrix of dimension \eqn{n} x 4 representing significant sides. Return whether each side has a significantly negative gradient (-1), significantly positive gradient (+1) or non-significant gradient (0). }
#' @references
#' Banerjee, S., and Gelfand, A. E. (2006). Bayesian wombling: Curvilinear gradient assessment under spatial process models. \emph{Journal of the American Statistical Association}, \bold{101}(476), 1487-1501.
#' @references
#' Goldstein, J., Jaewoo Park, Haran, M., Liebhold, A., and Bjornstad, O. N. (2018). Quantifying Spatio-Temporal Variation of Invasion Spread.\emph{arXiv preprint arXiv:1506.02685v3}
#' @details
#' This function directly tests if there is a significant radial expansion around a predicted point as in Goldstein et al. (2018), using the distribution of the gradient around a point (Banerjee and Gelfand, 2006).
#' When centered on each location, we can test the gradient normal to four sides of a box with sides of radius r.
#' As a heuristic we say if the spread is significantly outside of at least two sides of the box, and does not go significantly into any side of the box,
#' we can flag the location as a potential site of a long-range jump.
#'@seealso \code{\link{plotboxes}}
#' @export
#' @importFrom mapproj mapproject
#' @importFrom stats integrate
#' @importFrom cubature adaptIntegrate
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
#' # posterior mean of parameters
#' par.mean = apply(out.grad$post.samp,2,mean)
#'
#' # test for long range jump
#' out.jump = jump.scan(dates=dates, coord.longlat=coord, params = par.mean, Albers= c(29.5,45.5) )
#' }
##############################################################################################
jump.scan = function(dates, coord.longlat, newcoord.longlat=NULL, r=NULL, params, Albers = NULL, covar = NULL) {
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
  covariates = cbind(coord,covar)
  n.beta = dim(covariates)[2]+1
  all.sides = matrix(NA,ncol=4,nrow=dim(newcoord)[1])

  # compute pairwise distance matrix 'dmat' of coord
  dmat <- sqrt(outer(coord[,1], coord[,1], "-")^2 + outer(coord[,2], coord[,2], "-")^2)
  if(is.null(r)){   r <- quantile(dmat,.025)    }
  jump.data = list(newcoord=newcoord, r=r, params=params, dates=dates, coord=coord, n.beta=n.beta, covariates=covariates, dmat=dmat)

  for(i in 1:dim(newcoord)[1]) {
    all.sides[i,] = flux.box(newcoord[i,], jump.data)
  }
  res <- list(coord.longlat=coord.longlat,newcoord.longlat=newcoord.longlat, conv.factor=conv.factor, r=r, all.sides=all.sides)
  class(res) <- "jump.scan"
  res
}


#' Plot the significant long range jumps from the class \dQuote{\code{jump.scan}}.
#' @param r length of boxes or arrows on the plot. The default is 1.
#' @param object an object of class \code{jump.scan}, typically the result of a call to \code{\link{jump.scan}}.
#' @param num.sides number of significant sides in the box. It should be the integer between 1 and 4.
#' @param method the method to be used. If method is \dQuote{\code{box}}, draw a box for significant sides. If the method is \dQuote{\code{arrow}}, draw arrows for significant sides.
#' @param point if TRUE, plot points where long-range jumps are detected.
#' @param \dots additional arguments.
#' @details
#' This function draws boxes representing long range jumps on the plot using an object of class \code{jump.scan}.
#' @seealso \code{\link{jump.scan}}
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
#' # posterior mean of parameters
#' par.mean = apply(out.grad$post.samp,2,mean)
#'
#' # test for long range jump
#' out.jump = jump.scan(dates=dates, coord.longlat=coord, params = par.mean, Albers= c(29.5,45.5) )
#'
#' # gradient plot
#' plotgrad(out.grad,cex=1,pch=".",database="state")
#'
#' # long-range jump plot
#' plotboxes(object = out.jump, num.sides = 2, method="box")
#' }
##############################################################################################
plotboxes = function(r = NULL, object, num.sides, method = c("box", "arrow"),point=FALSE,...) {
  ##############################################################################################

  coord.longlat <- object[[1]]
  newcoord.longlat <- object[[2]]
  conv.factor <- object[[3]]
  jump <- object[[4]]
  all.sides <- object[[5]]
  if(is.null(r)){   r <- 1    }

  jump.centers = NULL
  jump.sides = NULL

  x0 = NULL; x1 = NULL; y0 = NULL; y1 = NULL;
  color.sides = NULL
  for(i in 1:dim(newcoord.longlat)[1]) {
    if(max(all.sides[i,]) <= 0 && sort(all.sides[i,])[num.sides] == -1) {
      jump.centers = rbind(jump.centers,newcoord.longlat[i,])
      jump.sides = rbind(jump.sides,all.sides[i,])
    }
  }
  if(length(jump.centers)==0) {
    print("No sites of potential long-range jumps.")
    return(0)
  }

  for(i in 1:dim(jump.centers)[1]) {
    center = jump.centers[i,]
    color.sides = c(color.sides,-jump.sides[i,] + 1) # black if nonsignificant, red if significant
    upper_x = center[1] + r
    lower_x = center[1] - r
    upper_y = center[2] + r
    lower_y = center[2] - r
    x0 = c(x0,lower_x,upper_x,lower_x,upper_x)
    x1 = c(x1,upper_x,lower_x,lower_x,upper_x)
    y0 = c(y0,upper_y,lower_y,lower_y,upper_y)
    y1 = c(y1,upper_y,lower_y,upper_y,lower_y)
    # plot in order top, bottom, left, right
  }
  color.sides[color.sides==1] = NA # remove nonsignificant segments  plotgrad(out.grad,cex=1,pch=".",database="state")

  if(method=="box"){               # Draw box
    segments(unlist(x0), unlist(y0), unlist(x1), unlist(y1), col=color.sides,...)}else{
      # Draw arrow
      points(jump.centers,cex=1,col=2,pch=7,lwd=2)
      colorsig <- matrix(color.sides,ncol=4,byrow=T)
      arrows(jump.centers[,1],jump.centers[,2],jump.centers[,1],y1[seq(1,length(color.sides),by=4)],col=colorsig[,1],length=0.1,...)  #1. (x,y+r)
      arrows(jump.centers[,1],jump.centers[,2],jump.centers[,1],y1[seq(2,length(color.sides),by=4)],col=colorsig[,2],length=0.1,...)  #2. (x,y-r)
      arrows(jump.centers[,1],jump.centers[,2],x1[seq(2,length(color.sides),by=4)],jump.centers[,2],col=colorsig[,3],length=0.1,...)  #3. (x-r,y)
      arrows(jump.centers[,1],jump.centers[,2],x1[seq(1,length(color.sides),by=4)],jump.centers[,2],col=colorsig[,4],length=0.1,...)  #4. (x+r,y)
    }

  # print
  if(conv.factor==1){
    cat("Long range jumps tested for" ,jump,".\n")
  }else{
    cat("Long range jumps tested for ",conv.factor*jump,"km.\n") # scale to km #
  }
}


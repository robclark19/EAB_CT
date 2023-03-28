\name{hemlock}
\alias{hemlock}
\docType{data}
\title{Spread of Hemlock wolly adelgid.}
\description{hemlock contains the hemlock data analyzed in Goldstein et al (2015).}
\usage{data("hemlock")}
\format{
  A data frame with 340 observations on the following 8 variables.
  \describe{
    \item{\code{FIPS}}{identification code.}
    \item{\code{first.year}}{the observed time that hemlock first appears.}
    \item{\code{long}}{longitude of each observation.}
    \item{\code{lat}}{lattitude of each observation.}
    \item{\code{ba}}{basal area of hemlock (Morin et al., 2004).}
    \item{\code{plant}}{plant hardiness zone (Cathey, 1990).}
    \item{\code{x}}{x-axis of Albers conic projected locations.}
    \item{\code{y}}{y-axis of Albers conic projected locations.}
  }
}

\source{Historical spread of the HWA was recorded at the county level. Records are from the US Forest Service Forest Health Protection.}

\references{
Cathey, H. M. (1990). USDA plant hardiness zone map. \emph{U.S.D.A. Misc. Publ.}


Goldstein, J., Haran, M., Bjornstad, O. N., and Liebhold, A. (2015). Quantifying Spatio-Temporal Variation of Invasion Spread.\emph{arXiv preprint arXiv:1506.02685}

Morin, R. S., Liebhold, A. M., Luzader, E. R., Lister, A. J., Gottschalk, K. W., and Twardus, D. B. (2004).Mapping host-species abundance of three major exotic forest pests.\emph{USDA Forest Service Northeastern Research Station Research Paper}
}
\examples{data(hemlock)}
\keyword{datasets}

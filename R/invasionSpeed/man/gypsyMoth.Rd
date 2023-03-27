\name{gypsyMoth}
\alias{gypsyMoth}
\docType{data}
\title{Spread of Gypsy moth.}
\description{gypsyMoth contains the Gypsy moth data analyzed in Goldstein et al (2015).}
\usage{data("gypsyMoth")}
\format{
  A data frame with 546 observations on the following 10 variables.
  \describe{
    \item{\code{FIPS}}{identification code.}
    \item{\code{first.year}}{the observed time that hemlock first appears.}
    \item{\code{lat}}{lattitude of each observation.}
    \item{\code{long}}{longitude of each observation.}
    \item{\code{size}}{    }
    \item{\code{population}}{    }
    \item{\code{income}}{    }
    \item{\code{basalperecent}}{percent of the forest basal area comprised of oaks (Liebhold et al., 1997).}
    \item{\code{x}}{x-axis of Albers conic projected locations.}
    \item{\code{y}}{y-axis of Albers conic projected locations.}
  }
}

\source{Historical spread of the Gypsy moth was recorded at the county level. Records are maintained by the United States Department of Agriculture.}

\references{

Goldstein, J., Jaewoo Park, Haran, M., Liebhold, A., and Bjornstad, O. N. (2018). Quantifying Spatio-Temporal Variation of Invasion Spread.\emph{arXiv preprint arXiv:1506.02685v2}


Liebhold, A. M., Gottschalk, K. W., Luzader, E. R., Mason, D. A., Bush, R., and Twardus, D. B. (1997). Gypsy moth in the United States, an atlas. \emph{General Technical Report-Northern Research Station, USDA Forest Service}


}
\examples{data(hemlock)}
\keyword{datasets}

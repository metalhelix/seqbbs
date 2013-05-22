\name{seqbbs}
\alias{seqbbs}
\title{Runs SeqBBS algorithm on pre-calculated ratios}
\usage{
  seqbbs(ratios, window = 12, threshold = 0.7)
}
\arguments{
  \item{ratios}{precomputed ratios between samples}

  \item{window}{window size to use.  Expected range: 12 -
  20}

  \item{threshold}{threshold value for the posterior
  probability.  Expected range: 0.55 - 0.95}
}
\value{
  SeqBBSData object with all data
}
\description{
  Details on this algorithm. found here:
}
\examples{
ratios <- c(0.2, 0.3, 0.4, 0.6)
seqbbs(ratios, window = 2, threshold = 0.8)
}


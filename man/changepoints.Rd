\name{changepoints}
\alias{changepoints}
\title{Pulls out change points from sebbs data given a threshold}
\usage{
  changepoints(seqbbs_data,
    threshold = seqbbs_data@threshold, confidence = 0.95)
}
\arguments{
  \item{seqbbs_data}{SeqBBSData object to work on.  Created
  using seqbbs function.}

  \item{threshold}{threshold value for the posterior
  probability.  Expected range: 0.55 - 0.95 Defaults to
  threshold in seqbbs}

  \item{confidence}{confidence interval value.  Defaults to
  0.95}
}
\value{
  dataframe of changepoint loci along with mean ratio and
  confidence interval
}
\description{
  Pulls out change points from sebbs data given a threshold
}
\examples{
ratios <- c(0.2, 0.3, 0.4, 0.6)
data <- seqbbs(ratios, window = 2, threshold = 0.8)
results <- changepoints(data)
}


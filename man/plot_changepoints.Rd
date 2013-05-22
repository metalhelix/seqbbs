\name{plot_changepoints}
\alias{plot_changepoints}
\title{Plots data along with changepoints at a particular threshold}
\usage{
  plot_changepoints(seqbbs_data,
    threshold = seqbbs_data@threshold, col = "red",
    pch = 1, basecol = "blue", basepch = 18,
    xlab = "Genomic Position",
    ylab = "Log2 Ratio of Reads", means = TRUE,
    meancol = "red", ...)
}
\arguments{
  \item{seqbbs_data}{SeqBBSData object to work on.  Created
  using seqbbs function.}

  \item{threshold}{threshold value for the posterior
  probability.  Expected range: 0.55 - 0.95 Defaults to
  threshold in seqbbs}

  \item{col}{color to use for changepoints}

  \item{pch}{symbol to use for changepoints}

  \item{basecol}{color to use for underlying data}

  \item{basepch}{symbol to use for underlying data}

  \item{xlab}{x axis label}

  \item{ylab}{y axis label}

  \item{means}{display mean lines or not}

  \item{meancol}{color to use for mean lines if displayed}

  \item{...}{other options passed into plot function}
}
\description{
  Plots data along with changepoints at a particular
  threshold
}


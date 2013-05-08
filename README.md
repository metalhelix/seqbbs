# seqbbs

R package to accompany paper:

## Install

To install the latest (potentially unstable) version, try using [devtools](https://github.com/hadley/devtools):

    install.packages("devtools")
    library(devtools)
    install_github('seqbbs')

## Use

`seqbbs` function runs main algorithm, returning a SeqBBSData object:

    ratios <- c(0.2, 0.3, 0.4, 0.6)
    data <- seqbbs(ratios, window = 2, threshold = 0.8)

You can plot that data directly using `plot_changepoints` and/or `plot_posteriors`:

    plot_changepoints(data)
    plot_posteriors(data)

You can also get a dataframe of changepoints, given a particular threshold

    results <- changepoints(data)


## Ratio Input

`seqbbs` requires a list of ratios. This can be produced from a BAM file in the following way:
    

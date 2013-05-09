# seqbbs

R package to accompany paper:

## Install

To install the latest (potentially unstable) version, try using [devtools](https://github.com/hadley/devtools):

    install.packages("devtools")
    library(devtools)
    install_github('seqbbs', username = 'metalhelix')

## Use

`seqbbs` function runs main algorithm, returning a SeqBBSData object:

    library('seqbbs')

    test_filename <- system.file("extdata", "paper.txt", package="seqbbs")
    ratios <- read.table(test_filename, header = FALSE)

    data <- seqbbs(ratios, window = 20, threshold = 0.7)

You can plot that data directly using `plot_changepoints` and/or `plot_posteriors`:

    plot_changepoints(data)
    plot_posteriors(data)

You can also get a dataframe of changepoints, given a particular threshold

    results <- changepoints(data)
    results

       loci log_mean_ratio mean_ratio confidence_lower_bound confidence_upper_bound
    1     0   -1.134965092  0.4553625              0.4508998              0.4598695
    2     2   -1.203170161  0.4381272              0.4322632              0.4440707
    3   258   -1.197257722  0.4382673              0.4312094              0.4454407
    4   358   -1.180306888  0.4429142              0.4258967              0.4606116
    5   371    0.815682310  1.8286183              1.6610245              2.0131221
    6   393   -0.055610043  1.0280860              0.9717306              1.0877097
    7   504    0.009137405  1.0092367              0.9946625              1.0240244
    8   576    0.057639917  1.0449468              1.0292259              1.0609078
    9   669    0.739912756  1.6817510              1.6426376              1.7217957
    10  736    1.052421293  2.0854567              1.9805270              2.1959457
    11  747    0.771789124  1.7145819              1.6640220              1.7666780
    12  772    0.287218176  1.2600919              1.1951205              1.3285954
    13  833    1.010216981  2.0331663              1.9259709              2.1463279
    14  850    1.290718032  2.4991125              2.3891859              2.6140967
    15  906    1.332075525  2.5510044              2.4862401              2.6174557
    16 1012    1.631469509  4.0091719              3.4669349              4.6362160
    17 1077    1.240876767  2.6971268              2.3733244              3.0651069
    18 1120    2.457895760  5.5510519              5.2918697              5.8229282
    19 1144    1.033351104  2.1127953              2.0178664              2.2121900
    20 1224    0.805114490  2.2791780              1.7965850              2.8914036
    21 1249    0.749060517  2.1260144              1.6836537              2.6846002
    22 1272    0.524096900  1.4422055              1.3933667              1.4927561
    23 1285   -0.373881710  0.7757135              0.7641330              0.7874695


## Ratio Input

`seqbbs` requires a list of ratios. This can be produced from a BAM file in the following way:
    

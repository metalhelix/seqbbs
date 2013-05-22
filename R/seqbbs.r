
setClass( "SeqBBSData",
          representation = representation(
            window = "numeric",
            threshold = "numeric",
            log_ratios = "numeric",
            all_posteriors = "matrix",
            max_posteriors = 'numeric',
            change_points = 'numeric'
            )
          )

#' Runs SeqBBS algorithm on pre-calculated ratios
#'
#' Details on this algorithm. found here:
#'
#' @param ratios precomputed ratios between samples
#' @param window window size to use.
#'    Expected range: 12 - 20  
#' @param threshold threshold value for the posterior probability.
#'    Expected range: 0.55 - 0.95
#' @export
#' @return SeqBBSData object with all data
#' @examples
#' ratios <- c(0.2, 0.3, 0.4, 0.6)
#' seqbbs(ratios, window = 2, threshold = 0.8)
seqbbs <- function(ratios, window = 12, threshold = 0.70) {
  
  # TODO: be flexible about input type or
  # convert to vector? initially
  
  # TODO: reduce memory by removing 
  # unneeded temp vars.

  # TODO: double check if this is right
  ratios <- as.vector(as.matrix(ratios))
  
  # original: n
  size <- length(ratios)

  # original: x
  log_ratios <- log2(ratios)

  # The sum of all reads ratios of a chromosome
  # original: x_s
  log_ratios_sum <- cumsum(log_ratios)
  
  # squared value of each observation in x
  # original: x2
  log_ratios_squared <- log_ratios ^ 2
  
  # cusum vector of x2
  # original: xx2_s
  log_ratios_squared_sum <- cumsum(log_ratios_squared)
  
  # The number of windows this sequence will 
  # have based on the given window size
  # original: d
  window_number <- floor(size / window)

  # Vector stores the window size for all the d+1 windows
  # original: wins
  windows <- rep(window,window_number + 1)
  
  # windows[1] is the first element of the window length vector.
  # we set it to an intial value of 0;
  windows[1] <- 0

  # This is the window size of the last window
  last_window_size <- size - window * (window_number - 1)  
  windows[window_number + 1] <- last_window_size
  
  # Vector with each element representing the cumulative length of the windows
  # original: swins
  windows_sum <- cumsum(windows)
  
  # Vector that contains the change point locations
  # original: cp
  change_points <- rep(0,window_number)
  
  # The maximum posterior probability vector, 
  # with each element being the max posterior of each segement, 
  # where all segements are divided by the sliding window
  # original: mp
  max_posteriors <- rep(0, window_number)
  
  # To remove the value of gamma function from infinity.
  # original: e
  e <- 1
  
  # The vector that contains the posterior probabilities 
  # of every position in a segement.
  # original: all_p
  all_posteriors <- t(rep(0,size))
  
  # TODO: vectorize
  
  for(k in 1:window_number) {
    
    # The number of observations contained in each window (segement)
    # original: n1
    window_obervations <- windows[k + 1]
    
    # The starting point of a window
    # original: k1
    window_start <- windows_sum[k] + 1
    
    # The ending point of the same window
    # original: k2
    window_end <- windows_sum[k + 1]
  
    # TODO: what do these variables mean?
    prob <- t(rep(0,size))
    var1 <- 0
    var2 <- 0
    A <- 0
    
    if(k == 1) {
      for(j in window_start:(window_end - 1)) {
        n1 <- window_end
        tau <- j
        # original: 
        #var1=(tau*xx2_s(tau)-(x_s(tau)^2));
        #var2=((n1-tau)*(xx2_s(n1)-xx2_s(tau))-(x_s(n1)-x_s(tau))^2);
        var1 <- (tau * log_ratios_squared_sum[tau] - (log_ratios_sum[tau] ^ 2))
        var2 <- ((n1 - tau) * (log_ratios_squared_sum[n1] - log_ratios_squared_sum[tau]) - (log_ratios_sum[n1] - log_ratios_sum[tau])^2)
        
        # original:
        #A=(var1^((tau-1)/2))*(var2^((n1-tau-1)/2));
        A <- (var1 ^ ((tau - 1) / 2)) * (var2 ^ ((n1 - tau - 1) / 2))
        
        # original:
        # prob(j)=(1/(n1-1))*gamma(e+(tau-1)/2)*gamma(e+(n1-tau-1)/2)*(tau^((tau-2)/2))*((n1-tau)^((n1-tau-2)/2))*(1/A);
        prob[j] <- (1 / (n1 - 1)) * gamma(e + (tau - 1) / 2) * gamma(e + (n1 - tau - 1) / 2) * (tau ^ ((tau - 2) / 2)) * ((n1 - tau) ^ ((n1 - tau - 2) / 2)) * (1 / A)
      }
      
    } else {
      
      qx_s <- log_ratios_sum[windows_sum[k] - 1]
      qxx2_s <- log_ratios_squared_sum[windows_sum[k] - 1]
      for(j in (window_start - 1):(window_end - 1)) {
        n1 <- window_end - window_start + 2
        tau <- j - window_start + 2
        # original:         
        #var1=(tau*(xx2_s(j)-qxx2_s)-(x_s(j)-qx_s)^2);
        #var2=((n1-tau)*(xx2_s(k2)-xx2_s(j))-(x_s(k2)-x_s(j))^2);
        var1 <- (tau * (log_ratios_squared_sum[j] - qxx2_s) - (log_ratios_sum[j] - qx_s) ^ 2)
        var2 <- ((n1 - tau) * (log_ratios_squared_sum[window_end] - log_ratios_squared_sum[j]) - (log_ratios_sum[window_end] - log_ratios_sum[j]) ^ 2)
        # original:
        #A=(var1^((tau-1)/2))*(var2^((n1-tau-1)/2));
        A <- (var1 ^ ((tau - 1) / 2)) * (var2 ^ ((n1 - tau - 1) / 2))
        # origianl:
        #prob(j)=(1/(n1-1))*gamma(e+(tau-1)/2)*gamma(e+(n1-tau-1)/2)*((tau)^((tau-2)/2))*((n1-tau)^((n1-tau-2)/2))*(1/A);
        prob[j] <- (1 / (n1 - 1)) * gamma(e + (tau - 1) / 2) * gamma(e + (n1 - tau - 1) / 2) * ((tau) ^ ((tau - 2) / 2)) * ((n1 - tau) ^ ((n1 - tau - 2) / 2)) * (1 / A)
      } #end for j
    } #end if k

    # original: p1
    p1 <- sum(prob)
  
    # the posterior probability vector
    # original: p
    p <- prob / p1
    
    # original: c1
    c1 <- which.max(p)
    # original: mp_1
    mp_1 <- p[c1]
    
    change_points[[k]] <- c1
    max_posteriors[[k]] <- mp_1
    
    if(k == 1) {
      for(j in window_start:(window_end - 1)) {
        all_posteriors[j] <- p[j]
      }
    } else {
      for(j in (window_start - 1):(window_end - 1)) {
        all_posteriors[j] <- p[j]
      }
    } #end if k
      
  } #end for k

  seqbbs_out = new("SeqBBSData", window = window,
                  threshold = threshold,
                  log_ratios = log_ratios, 
                  all_posteriors = all_posteriors,
                  max_posteriors = max_posteriors,
                  change_points = change_points)
  
  seqbbs_out
} #end function

#' Pulls out change points from sebbs data given a threshold
#'
#'
#' @param seqbbs_data SeqBBSData object to work on.
#'    Created using seqbbs function.
#' @param threshold threshold value for the posterior probability.
#'    Expected range: 0.55 - 0.95
#'    Defaults to threshold in seqbbs
#' @param confidence confidence interval value. 
#'    Defaults to 0.95
#' @export
#' @return dataframe of changepoint loci along with mean ratio and confidence interval
#' @examples
#' ratios <- c(0.2, 0.3, 0.4, 0.6)
#' data <- seqbbs(ratios, window = 2, threshold = 0.8)
#' results <- changepoints(data)
changepoints <- function(seqbbs_data, threshold = seqbbs_data@threshold, confidence = 0.95) {
  thresholded_iv <- seqbbs_data@max_posteriors >= threshold
  change_points_thresholded <- seqbbs_data@change_points[thresholded_iv]
  log_ratios_thresholded <- seqbbs_data@log_ratios[change_points_thresholded]
  cp_extended <- c(0, change_points_thresholded, length(seqbbs_data@log_ratios))
  
  log_ratios_thresholded_means = c()
  mean_ratios = c()
  lower_cis = c()
  upper_cis = c()
  for(k in 1:(length(change_points_thresholded) + 1)) {
    diff <- (cp_extended[k + 1] - cp_extended[k])
    sqrt_diff <- sqrt(diff)
    log_ratios_thresholded_means[k] <- sum(seqbbs_data@log_ratios[(cp_extended[k] + 1):cp_extended[k + 1]] / diff)
    threshold_sd <- sd(seqbbs_data@log_ratios[((cp_extended[k] + 1):cp_extended[k + 1])])
    mean_ratios[k] <- 2 ^ (log_ratios_thresholded_means[k] + ((threshold_sd)^2) / 2)
    
    lcl <- log_ratios_thresholded_means[k] - qnorm(0.5 + confidence / 2) * threshold_sd / sqrt_diff
    lower_cis[k] <- 2 ^ (lcl) * 2 ^ ((threshold_sd ^ 2) / 2)
    
    ucl <- log_ratios_thresholded_means[k] + qnorm(0.5 + confidence / 2) * threshold_sd / sqrt_diff
    upper_cis[k] <- 2 ^ (ucl) * 2 ^ ((threshold_sd ^ 2) / 2)
  }
  
  out <- data.frame('loci' = c(0,change_points_thresholded), 
                    'log_mean_ratio' = log_ratios_thresholded_means, 
                    'mean_ratio' =  mean_ratios,
                    'confidence_lower_bound' = lower_cis,
                    'confidence_upper_bound' = upper_cis)
  out
}

#' Plots data along with changepoints at a particular threshold
#'
#' @param seqbbs_data SeqBBSData object to work on.
#'    Created using seqbbs function.
#' @param threshold threshold value for the posterior probability.
#'    Expected range: 0.55 - 0.95
#'    Defaults to threshold in seqbbs
#' @param col color to use for changepoints
#' @param pch symbol to use for changepoints
#' @param basecol color to use for underlying data
#' @param basepch symbol to use for underlying data
#' @param xlab x axis label
#' @param ylab y axis label
#' @param means display mean lines or not
#' @param meancol color to use for mean lines if displayed
#' @param ... other options passed into plot function
#' @export
plot_changepoints <- function(seqbbs_data, 
                              threshold = seqbbs_data@threshold,
                              col = 'red',
                              pch = 1,
                              basecol = 'blue',
                              basepch = 18,
                              xlab = 'Genomic Position',
                              ylab = 'Log2 Ratio of Reads',
                              means = TRUE,
                              meancol = 'red', ...) {
    
  thresholded_changepoints <- changepoints(seqbbs_data, threshold = threshold)
  
  plot(1:length(seqbbs_data@log_ratios), seqbbs_data@log_ratios,  pch = basepch, col = basecol, ...)  
  points(thresholded_changepoints$loci, thresholded_changepoints$log_mean_ratio, pch = pch, col = col)
  
  if(means) {
    bars <- c(thresholded_changepoints$loci, length(seqbbs_data@log_ratios))
    for(k in 1:(length(thresholded_changepoints$loci))) {
      # xbar = sum(log_ratios[(bars[k] + 1):bars[k + 1]] / (bars[k + 1] - bars[k]))      
      # segments(bars[k] + 1, xbar, bars[k + 1], xbar, col = col) 
      xbar <- sum(seqbbs_data@log_ratios[(bars[k] + 1):bars[k + 1]] / (bars[k + 1] - bars[k]))
      segments(bars[k] + 1, xbar, bars[k + 1], xbar, col = meancol)
    }
  }  
}

#' Plots posteriors for SeqBBSData
#'
#' @param seqbbs_data SeqBBSData object to work on.
#'    Created using seqbbs function.
#' @param xlab x axis label
#' @param ylab y axis label
#' @param ... other options passed into plot function
#' @export
plot_posteriors <- function(seqbbs_data,
                            xlab = 'Genomic Position, 100kb',
                            ylab = 'Posterior Probs, all windows', ...) {
  barplot(seqbbs_data@all_posteriors, xlab = xlab, ylab = ylab, ...)
}

# ----

# threshold <- 0.70
# window <- 20
# 
# test_filename <- paste("inst","extdata", "paper.txt", sep="/")
# ratios <- read.table(test_filename, header = FALSE)
# 
# seqbbs_data <- seqbbs(ratios, window = window, threshold = threshold)
# 
# #class(seq_out)
# 
# opar <- par(mfcol=c(2,1), mar = c(0, 0, 0, 0) + 2)
# plot_changepoints(seqbbs_data)
# plot_posteriors(seqbbs_data)
# par(opar)
# 
# thresholded_changepoints <- changepoints(seqbbs_data, threshold = threshold)
# 
# signif(thresholded_changepoints$mean_ratio, digits = 4)


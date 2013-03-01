#'
#'
#'
setwd("~/code/R/seqbbs/")
#ratios <- read.table("test/test.txt", header = FALSE)

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

  # is the total number of breakpoints found for the given threshold
  # original: counter
  thresholded_iv <- max_posteriors >= threshold
  max_posteriors_thresholded <- max_posteriors[thresholded_iv]
  change_points_thresholded <- change_points[thresholded_iv]
  counter <- sum(thresholded_iv)

  opar <- par(mfcol=c(2,1), mar = c(0, 0, 0, 0) + 2)
  plot(1:size, log_ratios,  pch = 18, col = "blue")
  #points(change_points, log_ratios[change_points], pch = 1, col = "red")
  points(change_points_thresholded, log_ratios[change_points_thresholded], pch = 1, col = 'red')

  bars <- c(0, change_points_thresholded ,size)
  for(k in 1:(counter + 1)) {
    xbar = sum(log_ratios[(bars[k] + 1):bars[k + 1]] / (bars[k + 1] - bars[k]))
    segments(bars[k] + 1, xbar, bars[k + 1], xbar, col = 'red') 
  }
  
  barplot(all_posteriors)
  par(opar)
#   list(change_points, max_posteriors, all_posteriors)
  list(windows_sum, max_posteriors, change_points)
} #end function

ratios <- read.table("tests/test.txt", header = FALSE)
window <- 20
threshold <- 0.55
results <- seqbbs(ratios, window = window, threshold = threshold)




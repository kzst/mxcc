#-----------------------------------------------------------------------------#
#                                                                             #
#R Package for Maxwell Control Charts                             #
#                                                                             #
#  Written by: Zahid Khan, Zsolt T. Kosztyan                                  #
#              Department of Quantitative Methods                             #
#              University of Pannonia, Hungary                                #
#              kosztyan.zsolt@gtk.uni-pannon.hu                               #
#                                                                             #
# Last modified: April 2026                                                  #
#-----------------------------------------------------------------------------#

  #' @export
  mxewma <- function(y,
                     n = 4,
                     lambada = 0.20,
                     L = 3,
                     chart = c("V", "VSQ"),
                     summary = FALSE) {


    if (missing(y)) stop("Data vector/matrix 'y' is not specified")
    chart <- match.arg(chart)
    lambda <- lambada
    y <- stats::na.omit(y)
    n_int <- n * floor(length(y) / n)
    y <- y[1:n_int]
    yy <- matrix(y, ncol = n)
    m <- nrow(yy)

    v <- numeric(m)

    if (chart == "VSQ") {
      # ---- VSQ statistic ----
      for (j in 1:m) v[j] <- sqrt(sum(yy[j, ]^2) / (3 * n))

      a <- (sqrt(2) / sqrt(3 * n)) * (gamma((3 * n + 1) / 2) / gamma(3 * n / 2))
      sigma <- mean(v)

      # EWMA calculation
      En <- numeric(m)
      En[1] <- lambda * v[1] + (1 - lambda) * (a * sigma)
      for (i in 2:m) En[i] <- lambda * v[i] + (1 - lambda) * En[i - 1]

      # Control limits
      UCL <- numeric(m)
      LCL <- numeric(m)
      for (i in 1:m) {
        factor <- sqrt((lambda / (2 - lambda)) * (1 - (1 - lambda)^(2 * i)) * (1 - a^2))
        LCL[i] <- (a - L * factor) * sigma
        UCL[i] <- (a + L * factor) * sigma
      }
      CL <- a * sigma

    } else {
      # ---- V statistic ----
      for (j in 1:m) v[j] <- sum(yy[j, ]^2) / (3 * n)
      sigma <- mean(v)

      # EWMA calculation
      En <- numeric(m)
      En[1] <- lambda * v[1] + (1 - lambda) * (sigma^2)
      for (i in 2:m) En[i] <- lambda * v[i] + (1 - lambda) * En[i - 1]

      # Control limits
      UCL <- numeric(m)
      LCL <- numeric(m)
      for (i in 1:m) {
        factor <- sqrt((2/(3*n))*(lambda / (2 - lambda)) * (1 - (1 - lambda)^(2 * i)))
        LCL[i] <- (1 - L * factor) * (sigma^2)
        UCL[i] <- (1 + L * factor) * (sigma^2)
      }
      CL <- sigma^2
    }

    # Plot EWMA chart
    ii <- seq_len(m)
    plot(ii, En, type = "o", lty = 1, pch = 16,
         xlab = "Sample Number",
         ylab = paste("EWMA", chart, "Statistic"),
         main = (""),
         mgp = c(2,1,0),
         xlim = c(0, m + 1),
         ylim = range(c(En, UCL, LCL)),
         cex = 0.8)
    lines(ii, UCL, lty = 2, col = "green", lwd = 3)
    lines(ii, rep(CL, m), lty = 3)
    lines(ii, LCL, lty = 2, col = "green", lwd = 3)

    # Prepare output
    output <- list(
      LCL = LCL,
      UCL = UCL,
      CL = CL,
      Plotting_stat = En,
      Statistic = v,
      m = m,
      n = n,
      sigma = sigma,
      lambda = lambda,
      L = L,
      chart = chart
    )
    class(output) <- "mxewma"

    # Print summary if requested
    if (summary) {
      cat("EWMA Control Chart Summary\n")
      cat("Chart Type:", chart, "\n")
      cat("Subgroups (m):", m, " Sample size (n):", n, "\n")
      cat("Lambda:", lambda, " L:", L, "\n")
      cat("Estimated Sigma:", round(sigma, 4), "\n\n")
      cat("EWMA Statistics:\n"); print(round(En, 4))
      cat("\nLCL:\n"); print(round(LCL, 4))
      cat("\nUCL:\n"); print(round(UCL, 4))
    }

    return(invisible(output))
  }

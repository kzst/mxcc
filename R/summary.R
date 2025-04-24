#-----------------------------------------------------------------------------#
#                                                                             #
#R Package for Maxwell Control Charts                             #
#                                                                             #
#  Written by: Zahid Khan, Zsolt T. Kosztyan                                  #
#              Department of Quantitative Methods                             #
#              University of Pannonia, Hungary                                #
#              kosztyan.zsolt@gtk.uni-pannon.hu                               #
#                                                                             #
# Last modified: February 2025                                                  #
#-----------------------------------------------------------------------------#
#A Brief Summary of Control Chart Parameters for real (mxrpc) and simulated  (mxspc) data
#' @export
summary <- function(object, ...) {
  UseMethod("summary", object)
}

#' @export
summary.mxrpc <- function(object, ...) {
  cat("Summary of Control Chart Parameters:\n")
  cat("Subgroup Number (m):", object$m, "\n")
  cat("Sample Size (n):", object$n, "\n")
  cat(ifelse(object$limit == "PCL", "Lower Probability Limit (LPL):", "Lower Control Limit (LCL):"), round(object$LCL, 4), "\n")
  cat("Central Line (CL):", round(object$CL, 4), "\n")
  cat(ifelse(object$limit == "PCL", "Upper Probability Limit (UPL):", "Upper Control Limit (UCL):"), round(object$UCL, 4), "\n")
  cat("Estimated Sigma value:", round(object$sig, 4), "\n")
  cat("Limit Type:", object$limit, "\n")
  cat("Chart Type:", object$chart, "\n")
  cat("\nSummary Statistics for Plotting Statistic (V values):\n")
  print(summary(object$v))
  cat("\nSummary Statistics for Real Data:\n")
  print(summary(object$data))
}

#' @export
summary.mxspc <- function(object, ...) {
  cat("Summary of Control Chart Parameters:\n")
  cat("Subgroup Number (m):", object$m, "\n")
  cat("Sample Size (n):", object$n, "\n")
  cat(ifelse(object$limit == "PCL", "Lower Probability Limit (LPL):", "Lower Control Limit (LCL):"), round(object$LCL, 4), "\n")
  cat("Central Line (CL):", round(object$CL, 4), "\n")
  cat(ifelse(object$limit == "PCL", "Upper Probability Limit (UPL):", "Upper Control Limit (UCL):"), round(object$UCL, 4), "\n")
  cat("Estimated Sigma value:", round(object$sig, 4), "\n")
  cat("Limit Type:", object$limit, "\n")
  cat("Chart Type:", object$chart, "\n")
  cat("\nSummary Statistics for Plotting Statistic (V values):\n")
  print(summary(object$v))
  cat("\nSummary Statistics for Simulated Data (x values):\n")
  print(summary(object$a))
}


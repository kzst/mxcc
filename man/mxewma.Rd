\name{mxewma}
\alias{mxewma}
\title{
EWMA Control Chart for V and VSQ Statistics
}
\description{
The \code{mxewma} function constructs an Exponentially Weighted Moving Average (EWMA) control chart
for monitoring the V and VSQ statistics of a process.
It calculates subgroup statistics, EWMA plotting statistics, control limits, and optionally displays
a summary including estimated sigma.
}
\usage{
mxewma(y, n = 4, lambada = 0.20, L = 3, chart = c("V", "VSQ"), summary = FALSE)
}

\arguments{
  \item{y}{
A numeric vector of process observations.
}
  \item{n}{
Subgroup size. Default is 4.
}
  \item{lambada}{
EWMA smoothing parameter (0 < \code{lambda} <= 1). Default is 0.20.
}
  \item{L}{
Control limit multiplier. The value of \code{L} should be selected based on \code{lambda} according to established quality control references. Default is 3.
}
  \item{chart}{
Type of chart to construct: "V"  or "VSQ" for square-rooted statistic. Default is c("V", "VSQ").
}
  \item{summary}{
Logical value indicating whether to print a summary of the chart, including estimated sigma. Default is \code{FALSE}.
}
}

\details{
This function divides the input data \code{y} into subgroups of size \code{n},
and computes the V or VSQ statistic for each subgroup based on the
single-parameter Maxwell distribution.

The Maxwell distribution is assumed as the underlying model for the process,
and the subgroup statistics are constructed accordingly. The EWMA statistic is
then computed to monitor small and moderate shifts in the process.The appropriate value of the multiplier L
can be obtained from the cited references for selected values of lambda
and sample size (n) to achieve an in-control Average Run Length (ARL0) of 370.

Dynamic control limits (LCL and UCL) are calculated at each subgroup using the
EWMA formulation. For the VSQ chart, a constant
\eqn{a = (\sqrt{2}/\sqrt{3n}) (\Gamma((3n+1)/2)/\Gamma(3n/2))}
is used to adjust the center line and control limits.

The estimated sigma is obtained as the mean of the subgroup statistics.
The value of \code{L} should be selected corresponding to the smoothing
parameter \code{lambda} based on established EWMA design tables available
in the literature.

When \code{summary = TRUE}, the function prints a textual summary including
the EWMA statistics, control limits, and the estimated sigma value.
}

\value{
An invisible list containing the following components:
\itemize{
  \item \code{LCL}: Lower Control Limit vector.
  \item \code{UCL}: Upper Control Limit vector.
  \item \code{CL}: Center Line.
  \item \code{Plotting_stat}: EWMA statistics for each subgroup.
  \item \code{Statistic}: Original subgroup statistic (V or VSQ).
  \item \code{m}: Number of subgroups.
  \item \code{n}: Subgroup size.
  \item \code{sigma}: Estimated sigma (mean of subgroup statistics).
  \item \code{lambda}: EWMA smoothing parameter.
  \item \code{L}: Control limit multiplier.
  \item \code{chart}: Chart type selected.
}
If \code{summary = TRUE}, a textual summary is also printed.
}

\references{
Hossain, M.P. and Riaz, M. (2021) "On designing a new VEWMA control chart for efficient process monitoring". Computers & Industrial Engineering, 162, 107751. <doi:10.1016/j.cie.2021.107751>

Khan, Z., Saghir, A., Katona, A. and Kosztyán, Z.T. (2025) "EWMA control chart framework for efficient Maxwell quality characteristic monitoring: An application to the aerospace industry". Computers & Industrial Engineering, 200, 110753. <doi:10.1016/j.cie.2024.110753>
}
\author{
Zahid Khan
}

\seealso{
\code{\link{mxrpc}}, \code{\link{plot}}, \code{\link{summary}}
}

\examples{
y <- c(
1.5391137, 2.2764737, 1.7794763, 2.3362682, 1.6382138, 1.6943709, 1.3135375, 0.5767501,
1.4950770, 2.2973635, 1.9949414, 2.6226897, 0.8730001, 2.1368288, 0.7282481, 1.7775871,
2.3898929, 2.3293989, 1.8709563, 2.5531137, 1.8444661, 1.7220544, 0.5774853, 1.6022162,
1.0952978, 0.4020387, 1.5097896, 1.8443609, 1.5732371, 0.8237968, 1.6290878, 1.0708337,
1.6846298, 0.8824187, 1.1953376, 1.6381853, 0.7522207, 3.2879978, 0.9592238, 0.7797138
)
mxewma(y, n = 4, lambada = 0.2, L = 3, chart = "VSQ", summary = FALSE)

mxewma(y, n = 4, lambada = 0.2, L = 3, chart = "V", summary = FALSE)
}

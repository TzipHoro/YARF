% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/YARF_asynchronicity.R
\name{YARF_progress_reports}
\alias{YARF_progress_reports}
\title{Prints out a messages reflecting the progress of the YARF model construction until completion}
\usage{
YARF_progress_reports(yarf_mod, time_delay_in_seconds = 10,
  plot_oob_error = FALSE, trail_pts = 5)
}
\arguments{
\item{yarf_mod}{The yarf model object}

\item{time_delay_in_seconds}{Frequency of messages in seconds. Default is \code{10} seconds.}

\item{plot_oob_error}{Create a plot of oob error rate as trees are built to assess 
convergence of the YARF model}

\item{trail_pts}{If \code{plot_oob_error} is \code{TRUE}, this optional parameter
if non null will plot a secondary window with only the trees in the last \code{trail_trees}
number of progress updates defined by \code{time_delay_in_seconds} 
versus error which allows assessing the convergence more closely.
Default is \code{5}.}
}
\description{
Prints out a messages reflecting the progress of the YARF model construction until completion
}
\author{
Kapelner
}


\docType{methods}
\name{plot,criterion-method}
\alias{plot,criterion-method}
\alias{plot.criterion}
\title{Plot method for object criterion}
\usage{
\\S4method{plot}{criterion}(x, y, xvar="lambda1",
log.scale=ifelse(xvar=="lambda1",TRUE,FALSE), plot=TRUE, ...)
}
\arguments{
  \item{x}{an object with class criterion}

  \item{y}{used for S4 compatibility.}

  \item{xvar}{a character among \code{"lambda1"},
  \code{"fraction"} or \code{"df"} for the x-axis (only
  \code{lambda1} is available when the criterion object
  come from a call to \code{cv.spring}).}

  \item{log.scale}{a boolean to put the x-axis on a
  log-scale.}

  \item{plot}{is \code{TRUE}, a plot is produce. Otherwise
  a \code{ggplot} object is sent back.}

  \item{...}{Arguments to be passed to methods, such as
  \link{graphical parameters} (see \code{\link{par}}).
  Many methods will accept the following arguments:
  \describe{ \item{\code{type}}{what type of plot should be
  drawn.  Possible types are \itemize{ \item \code{"p"} for
  \bold{p}oints, \item \code{"l"} for \bold{l}ines, \item
  \code{"b"} for \bold{b}oth, \item \code{"c"} for the
  lines part alone of \code{"b"}, \item \code{"o"} for both
  \sQuote{\bold{o}verplotted}, \item \code{"h"} for
  \sQuote{\bold{h}istogram} like (or \sQuote{high-density})
  vertical lines, \item \code{"s"} for stair \bold{s}teps,
  \item \code{"S"} for other \bold{s}teps, see
  \sQuote{Details} below, \item \code{"n"} for no plotting.
  } All other \code{type}s give a warning or an error;
  using, e.g., \code{type = "punkte"} being equivalent to
  \code{type = "p"} for S compatibility.  Note that some
  methods, e.g. \code{\link{plot.factor}}, do not accept
  this.  } \item{\code{main}}{an overall title for the
  plot: see \code{\link{title}}.} \item{\code{sub}}{a sub
  title for the plot: see \code{\link{title}}.}
  \item{\code{xlab}}{a title for the x axis: see
  \code{\link{title}}.} \item{\code{ylab}}{a title for the
  y axis: see \code{\link{title}}.} \item{\code{asp}}{the
  \eqn{y/x} aspect ratio, see \code{\link{plot.window}}.} }
  }
}
\value{
a invisible ggplot object.
}
\description{
A plotting method for the class \code{criterion}.
}
\seealso{
\code{\linkS4class{criterion}}.
}


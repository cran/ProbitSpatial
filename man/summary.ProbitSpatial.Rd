% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ProbitSpatial.R
\name{summary.ProbitSpatial}
\alias{summary.ProbitSpatial}
\title{Spatial probit model summaries.}
\usage{
\method{summary}{ProbitSpatial}(object, covar = FALSE, ...)
}
\arguments{
\item{object}{an object of class \code{ProbitSpatial}.}

\item{covar}{should the statistics be computed with the matrix of
variance of the parametes or not. Default is FALSE, hence Likelihood-ratio 
    statistics are printed.}

\item{...}{further arguments}
}
\value{
This functions does not return any value.
}
\description{
Print the results of a \code{ProbitSpatial} model.
}
\details{
The \code{summary} function prints
\describe{
  \item{Model}{Featurs on the model and dataset.}
  \item{Time}{Estimation time.}
  \item{Statistics}{Standard errors of the estimated parameters. If 
\code{covar=TRUE}, it uses the matrix of variance of the parameters, else the 
likelihood ratio test.}
  \item{Accuracy}{Confusion Matrix and accuracy of the estimated model.}
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{yget_args}
\alias{yget_args}
\title{Extract vars in the parent frame, then filtered vars by a given function's formals}
\usage{
yget_args(..., .f = NULL)
}
\arguments{
\item{...}{filter all the args by the formals of .f, return all matched,}

\item{.f}{.filter, is.function() and is the filter function, the returned args list will match
the formals of the .f}
}
\value{
like `list(...)[intersect(list(...)%>%names,formals(.f)%>%names)]`
}
\description{
Get the arguments whose name match the given function(.f)'s formal arg names in the input args
lists
}

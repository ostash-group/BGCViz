% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_format_transformation.R
\name{sempi_to_csv}
\alias{sempi_to_csv}
\title{sempi_to_csv}
\usage{
sempi_to_csv(project_archive, write_to = getwd())
}
\arguments{
\item{project_archive}{- path to project.zip file, downloaded from SEMPI}

\item{write_to}{- path where to write generated csv file}
}
\value{
csv file in specified location
}
\description{
Function, which transforms Track.db file into dataframe, which could be then written to csv.
Download project folder from SEMPI and supply as `project_archive` argument to a function
}
\examples{
\dontrun{
sempi_to_csv(<zip-file>)
}
}

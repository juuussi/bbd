% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis_functions.R
\name{summarize_data}
\alias{summarize_data}
\title{Summarize data}
\usage{
summarize_data(data, row.checks = c("duplicated_row"), col.checks = c("na",
  "infinite", "nan", "empty_string", "whitespace",
  "leading_or_trailing_whitespace", "byte_sequence_character",
  "unicode_replacement_character", "linebreak", "excel_formula_error",
  "comma_as_decimal_separator", "duplicated_column"), parallel = FALSE)
}
\arguments{
\item{data}{a data.frame to process}

\item{row.checks}{row specific checks to run.}

\item{col.checks}{column specific checks to run.}

\item{parallel}{logical flag indicating if the tasks should be run in parallel. The default value is \code{FALSE}.}
}
\value{
An object of class "summarize_data" with the following components: (TBA)
}
\description{
Summarize and check \code{data.frame} for potential data errors.
}
\examples{
summarize_data(iris)

# Add values that will cause warnings in checks
tmp <- iris
tmp[3,4] <- " "
tmp[4,2] <- ""
tmp[1,2] <- NA
tmp[8,2] <- " test"
tmp[1,2] <- Inf
tmp[1,4] <- NaN
tmp[18,1] <- "#NAME?"

summarize_data(tmp)
}


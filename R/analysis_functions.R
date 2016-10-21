
#' Test of correlation for multiple variables between paired samples
#'
#' @description Test of correlations for multiple variables between paired samples, using \code{cor.test} function.
#'
#' @param data a data frame with the data for all variables to be correlated.
#' @param x.variables character vector of column names for variables to correlate to \code{variables.y}. If \code{variables.y} is left empty, correlate all \code{variables.x} to each other.
#' @param y.variables character vector of column names for variables to correlate to \code{variables.x}. If left \code{NULL}, correlate all \code{variables.x} to each other.
#' @param parallel logical flag indicating if the tasks should be run in parallel. The default value is \code{FALSE}.
#' @param parallel.inorder logical flag indicating if parallelry run tasks should be combined in the same order that they were submitted. If the order is not important, then setting \code{parallel.inorder} to \code{FALSE} can result in improved performance. The default value is \code{FALSE}.
#' @param ... further arguments to be passed to \code{cor.test} method (e.g. \code{method}, \code{conf.level}, \code{alternative}, \code{exact}).
#'
#' @return A data.frame with correlation results for each pair.
#' @export
#'
#' @seealso \code{\link{cor.test}} for correlation testing function and its arguments.
#'
#' @examples
#' correlations(data=iris, x.variables=c("Sepal.Length", "Sepal.Width",
#'              "Petal.Length", "Petal.Width"), parallel=FALSE)
#'
#' correlations(data=iris, x.variables=c("Sepal.Length", "Sepal.Width"),
#'              y.variables=c("Petal.Length", "Petal.Width"), parallel=FALSE)
#'
correlations <- function(data, x.variables, y.variables=NULL, parallel=FALSE, parallel.inorder=FALSE, ...) {

  # Toggle for logging/availability of futile.logger package
  f_logging <- requireNamespace("futile.logger", quietly = TRUE)

  if (f_logging) futile.logger::flog.debug("Initializing correlation analysis.")

  if (is.null(x.variables) || !x.variables %in% colnames(data)) {
    stop("x.variables NULL or not found in data column names.")
  }
  if (!is.null(y.variables) && !y.variables %in% colnames(data)) {
    stop("y.variables not found in data column names.")
  }

  # Select serial or parallel operator for foreach-function
  if (parallel) {
    `%do_operator%` <- foreach::"%dopar%"
  } else {
    `%do_operator%` <- foreach::"%do%"
  }

  `%nesting_operator%` <- foreach::"%:%"

  # Create list of variable pairs to correlate
  if (is.null(y.variables)) {
    # x.variables paired with each other
    pairs <- foreach::foreach(i=1:(length(x.variables)-1), .combine="c", .inorder=parallel.inorder) %nesting_operator%
      foreach::foreach(j=(i+1):length(x.variables)) %do_operator% {
        c(x.variables[i], x.variables[j])
      }
  } else {
    # x.variables paired with y.variables
    pairs <- foreach::foreach(i=1:length(x.variables), .combine="c", .inorder=parallel.inorder) %nesting_operator%
      foreach::foreach(j=1:length(y.variables), .inorder=parallel.inorder) %do_operator% {
        c(x.variables[i], y.variables[j])
      }
  }

  total.correlations <- length(pairs)

  if (f_logging) futile.logger::flog.debug(paste0("data: [ ", nrow(data), " x ", ncol(data), " ] Number of x.variables: ", length(x.variables), ". Number of y.variables: ", length(y.variables), ". Number of total correlations: ", total.correlations))
  if (f_logging) futile.logger::flog.trace(paste0("x.variables: '", paste(x.variables, collapse="', '"), "'. y.variables: ", paste(y.variables, collapse="', '"), "'"))
  if (f_logging) futile.logger::flog.trace(paste0("Correlation pairs: ", pairs))

  results <- foreach::foreach(i=1:length(pairs), .combine="rbind", .inorder=parallel.inorder) %do_operator% {
    x <- pairs[[i]][1]
    y <- pairs[[i]][2]

    n <- length(which(complete.cases(data[,c(x, y)])))
    if (f_logging) futile.logger::flog.trace(paste0("Correlating \'", x, "\' to \'", y, "\'. Correlation #: ", i, "/", total.correlations))
    result.row <- NULL
    tryCatch({
      c <- cor.test(x=data[,x], y=data[,y], ...)
      result.row <- data.frame(Variable1=x, Variable2=y, N=n, L_CI=c$conf.int[1], U_CI=c$conf.int[2], Estimate=c$estimate, P=c$p.value)
    }, error = function(msg) {
      warning.msg <- paste0("Problem while correlation variables: \'", x, "\' and \'", y, "\'. ", msg)
      warning(warning.msg)
      if (f_logging) futile.logger::flog.warn(warning.msg)
    }
    )
    if (is.null(result.row)) {
      result.row <- data.frame(Variable1=x, Variable2=y, N=n, L_CI=NA, U_CI=NA, Estimate=NA, P=NA)
    }
    result.row
  }
  results$PFDR <- p.adjust(results$P, method="BH")
  rownames(results) <- NULL
  if (f_logging) futile.logger::flog.debug("Finished correlation analysis.")
  results
}

#' Summarize data
#'
#' @description Summarize and check \code{data.frame} for potential data errors.
#'
#' @param data a data.frame to process
#' @param checks checks to perform.
#' @param parallel logical flag indicating if the tasks should be run in parallel. The default value is \code{FALSE}.
#'
#' @return An object of class "summarize_data" with the following components: (TBA)
#'
#' @export summarize_data
#'
#' @examples
#' summarize_data(iris)
#'
#' # Add values that will cause warnings in checks
#' tmp <- iris
#' tmp[3,4] <- " "
#' tmp[4,2] <- ""
#' tmp[1,2] <- NA
#' tmp[8,2] <- " test"
#' tmp[1,2] <- Inf
#' tmp[1,4] <- NaN
#' tmp[18,1] <- "#NAME?"
#'
#' summarize_data(tmp)

summarize_data <- function(data, checks=c("na", "infinite", "nan", "empty_string", "whitespace", "leading_or_trailing_whitespace", "byte_sequence_character", "unicode_replacement_character", "linebreak", "excel_formula_error"), parallel=FALSE) {

  # Select serial or parallel operator for foreach-function
  if (parallel) {
    `%do_operator%` <- foreach::"%dopar%"
  } else {
    `%do_operator%` <- foreach::"%do%"
  }

  classes <- sapply(X=data, FUN=class)
  classes_df <- as.data.frame(table(classes))
  class_freq <- classes_df$Freq
  names(class_freq) <- classes_df$classes

  checks <- foreach::foreach (i=1:length(checks), .combine=c, .inorder=FALSE) %do_operator% {
    check <- checks[i]
    if (check == "na") {
      result <- list(na=which(colSums(sapply(X=data, FUN=function(x) { is.na(x) })) > 0))
    } else if (check == "infinite") {
      result <-  list(infinite=which(colSums(sapply(X=data, FUN=function(x) { !is.na(x) & x == Inf })) > 0))
    } else if (check == "nan") {
      result <-  list(nan=which(colSums(sapply(X=data, FUN=function(x) { !is.na(x) & is.nan(x) })) > 0))
    } else if (check == "empty_string") {
      result <- list(empty_string=which(colSums(sapply(X=data, FUN=function(x) { !is.na(x) & x == "" })) > 0))
    } else if (check == "whitespace") {
      result <- list(whitespace=which(sapply(X=data, FUN=function(x) { length(grep("^[[:space:]]+$", x))}) > 0))
    } else if (check == "leading_or_trailing_whitespace") {
      result <- list(leading_or_trailing_whitespace=which(sapply(X=data, FUN=function(x) { length(grep("(^\\s+\\S+)|(\\S+\\s+$)", x, perl=TRUE))}) > 0))
    } else if (check == "byte_sequence_character") {
      result <- list(byte_sequence_character=which(sapply(X=data, FUN=function(x) { length(grep("[\\]x", x))}) > 0))
    } else if (check == "unicode_replacement_character") {
      result <- list(unicode_replacement_character=which(sapply(X=data, FUN=function(x) { length(grep("\xEF\xBF\xBD", x))}) > 0))
    } else if (check == "linebreak") {
      result <- list(linebreak=which(sapply(X=data, FUN=function(x) { length(grep("(\r)|(\n)", x))}) > 0))
    } else if (check == "excel_formula_error") {
      result <- list(excel_formula_error=which(sapply(X=data, FUN=function(x) { length(grep("^((#DIV/0!)|(#N/A)|(#NAME\\?)|(#NULL!)|(#NUM!)|(#REF!)|(#VALUE!)|(#GETTING_DATA))$", x))}) > 0))
    }
    result
  }

  structure(list(dimensions=dim(data), classes=classes, class_freq=class_freq, checks=checks), class="summarized.data")
}

#' Printing summarize_data
#'
#' @description Print a \code{summarized.data} object.
#'
#' @usage
#' ## S3 method for class 'summarized.data'
#'
#' @param x object of class \code{summarize_data}
#'
#' @export print.summarized.data
#'
#' @seealso \code{\link{summarize_data}}
#'
#' @examples
#' print(summarize_data(iris))
#'
print.summarized.data <- function(x, ...) {
  cat(paste0("Dimensions: ",  x$dimensions[1], " (rows) x ", x$dimensions[2], " (cols). Total observations: ", x$dimensions[1] * x$dimensions[2], "\n"))
  cat(paste0("Column classes:\n - ", paste0(names(x$class_freq), ": ", x$class_freq, collapse=", "), "\n"))
  warning_text_displayed <- FALSE
  for (name in names(x$checks)) {
    warning_hits <- x$checks[[name]]
    if (length(warning_hits) > 0) {
      if (!warning_text_displayed) {
        cat("Warnings (amount of columns with possible errors):\n")
        warning_text_displayed <- TRUE
      }
      cat(paste0(" - ", name, ": ", length(warning_hits), collapse=" "), "\n")
    }
  }
}

#' Convert summarize_data object into a molten data frame.
#'
#' @description This function melts summarize_data objects into a long-format \code{data.frame}
#' for easier processing and visualization.
#'
#' @param x summarized.data object to melt
#' @param parallel logical flag indicating if the tasks should be run in parallel. The default value is \code{FALSE}.
#'
#' @return molten \code{data.frame}
#' @export melt_summarized_data
#'
#' @examples
#' # Add random errors to data
#' tmp <- mtcars
#' tmp[sample(nrow(tmp),sample(1:10, 1)), sample(ncol(tmp), 3)] <- NA
#' tmp[sample(nrow(tmp),sample(1:10, 1)), sample(ncol(tmp), 2)] <- ""
#'
#' s <- summarize_data(tmp)
#' m <- melt_summarized_data(s)
#'
#' \dontrun{
#'
#' # Visualize the classes and warnings with histograms
#' library(ggplot2)
#' p <- ggplot(data=m, aes(x=Column)) + geom_histogram(col="black", aes(fill=Type))
#' p <- p + facet_wrap(facets = "Label", ncol=1, scales="free_y")
#' p <- p + theme_minimal() + scale_fill_manual(values=c("steelblue", "sienna1"))
#' p
#' }


melt_summarized_data <- function(x, parallel=FALSE) {

  # Select serial or parallel operator for foreach-function
  if (parallel) {
    `%do_operator%` <- foreach::"%dopar%"
  } else {
    `%do_operator%` <- foreach::"%do%"
  }


  df <- data.frame(Column=1:length(x$classes), Type="class", Name=x$classes, Label=paste0("Class: ", x$classes), Value=1)

  df2 <- foreach::foreach(i=1:length(x$checks), .combine=rbind, .inorder=FALSE) %do_operator% {
    name <- names(x$checks[i])
    c <- x$checks[[i]]

    result.rows <- NULL
    if (length(c) > 0) {
      result.rows <- data.frame(Column=c, Type="warning", Name=name, Label=paste0("Warning: ", name), Value=1)
    }
    result.rows
  }

  rbind(df, df2)

}


#' Test of correlation for multiple variables between paired samples
#'
#' @description Test of correlations for multiple variables between paired samples, using \code{cor.test} function.
#'
#' @param data a data frame with the data for all variables to be correlated.
#' @param x.variables character vector of column names for variables to correlate to \code{variables.y}. If \code{variables.y} is left empty, correlate all \code{variables.x} to each other.
#' @param y.variables character vector of column names for variables to correlate to \code{variables.x}. If left \code{NULL}, correlate all \code{variables.x} to each other.
#' @param parallel logical flag indicating if the tasks should be run in parallel. The default value is \code{TRUE}.
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
correlations <- function(data, x.variables, y.variables=NULL, parallel=TRUE, parallel.inorder=FALSE, ...) {

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

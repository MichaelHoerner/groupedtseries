#' Fit grouped GARCH model to time series
#'
#' Fit a Generalized Autoregressive Conditional Heteroscedastic GARCH(p,q)
#' time series model to the data by computing the maximum-likelihood
#' estimates of the conditionally normal model
#' to a group of univariate time series by finding coefficients
#' that serve best the whole group.
#' @param x a two-dimensional numeric vector of a grouped time series with accompanying group index in second dimension.
#' One-dimensional if input shall be only a single time series.
#' @param order a vector c(p,q) that indicates the order p of ARCH and q of GARCH part
#' @param series name for the series. Defaults to \code{deparse(substitute(x))}
#' @param control a list of control parameters as set up by \code{garch.control}
#' @param maxiter gives the maximum number of log-likelihood function
#' evaluations \code{maxiter} and the maximum number of iterations
#' \code{2*maxiter} the optimizer is allowed to compute.
#' @param trace logical. Trace optimizer output?
#' @param start If given this numeric vector is used as the initial estimate
#' of the GARCH coefficients.  Default initialization is to set the
#' GARCH parameters to slightly positive values and to initialize the
#' intercept such that the unconditional variance of the initial GARCH
#' is equal to the variance of \code{x}
#' @param grad character indicating whether analytical gradients or a numerical
#' approximation is used for the optimization.
#' @param abstol absolute function convergence tolerance.
#' @param reltol relative function convergence tolerance.
#' @param xtol coefficient-convergence tolerance.
#' @param falsetol false convergence tolerance.
#' @param \dots additional arguments for \code{\link{qr}} when computing
#' the asymptotic covariance matrix.
#'
#' @details \code{garch} uses a Quasi-Newton optimizer to find the maximum
#' likelihood estimates of the conditionally normal model.  The first
#' max(p, q) values are assumed to be fixed.  The optimizer uses a hessian
#' approximation computed from the BFGS update.  Only a Cholesky factor
#' of the Hessian approximation is stored.  For more details see Dennis
#' et al. (1981), Dennis and Mei (1979), Dennis and More (1977), and
#' Goldfarb (1976).  The gradient is either computed analytically or
#' using a numerical approximation.
#' @return Object of class \code{GARCH}.
#' \item{order}{the order of the fitted model.}
#' \item{coef}{coef estimated GARCH coefficients for the fitted model across all grouped time series.}
#' \item{n.likeli}{the negative log-likelihood function evaluated at the
#' coefficient estimates (apart from some constant).}
#' \item{n.used}{the number of observations of \code{x}.}
#' \item{residuals}{the series of residuals.}
#' \item{fitted.values}{the bivariate series of conditional standard
#' deviation predictions for \code{x}.}
#' \item{series}{the name of the series \code{x}.}
#' \item{frequency}{the frequency of the series \code{x}.}
#' \item{call}{the call of the \code{garch} function.}
#' \item{vcov}{outer product of gradient estimate of the asymptotic-theory
#' covariance matrix for the coefficient estimates.}
#'
#'
#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export



grouped_garch <- function (x, order = c(1, 1), series = NULL, control = garch.control(...), ...)
{
  index.segments <- function() {

    if (length(unique(x[,2])) == 1) {
      index_start_segment <- 1
      index_end_segment <- nrow(x)
    } else {
      index_start_segment <- c(1,1+which(diff(x[,2])!=0))
      index_end_segment <- c(index_start_segment[-1] - 1, nrow(x))
    }

    if (any(is.na(x[,1]))) {
      index_na <- which(is.na(x[,1]))
      index_end_na <- index_na - 1
      index_start_na <- index_na + 1

      #cancel all entries of consecutive na's
      index_end_na <- index_end_na[!index_end_na %in% index_na]
      index_start_na <- index_start_na[!index_start_na %in% index_na]

      #delete indizes if they are out of bound
      if (max(index_start_na) > nrow(x))
        index_start_na <- index_start_na[-which(index_start_na > nrow(x))]
      if (min(index_end_na) < 1)
        index_end_na <- index_end_na[-which(index_end_na < 1)]
    } else {
      index_na <- NULL
      index_end_na <- NULL
      index_start_na <- NULL
    }

    #cancel all entries for start/end patient where overlapping with na
    index_end_segment <- index_end_segment[!index_end_segment %in% index_na]
    index_start_segment <- index_start_segment[!index_start_segment %in% index_na]

    #concat index vectors and delete duplicates
    index_start <- unique(sort(c(index_start_na, index_start_segment)))
    index_end <- unique(sort(c(index_end_na, index_end_segment)))

    #create indizes without NA gaps
    index_start_wo_na <- index_start
    index_end_wo_na <- index_end
    if (length(index_start) > 1) {
      for (i in 2:length(index_start)) {
        incr <- index_start[i] - index_end[i-1] - 1
        index_start_wo_na[i:length(index_start)] <- index_start_wo_na[i:length(index_start)] - incr
        index_end_wo_na[i:length(index_start)] <- index_end_wo_na[i:length(index_start)] - incr
      }
    }

    ts_segments <- data.frame(cbind(index_start, index_end, index_start_wo_na, index_end_wo_na, index_end - index_start + 1))
    colnames(ts_segments) <- c("index_start", "index_end", "index_start_wo_na", "index_end_wo_na", "length")

    return(ts_segments)
  }

    if (NCOL(x) > 2) {
      stop("x is not a vector or univariate time series with indizes.
               There is at least a third column.")
    } else if (NCOL(x) == 2 && (any(is.na(x[,2])) || any(x[,2] == 0))) {
      stop("0 or NA are no valid entries for time series index.")
    } else if (NCOL(x) == 1){

      x_wo_na <- na.omit(x)
      index_x <- c(rep(1, length(x_wo_na)))
      x_wo_na <- cbind(x_wo_na, index_x)
      ts_segments <- index.segments()
    } else {
      x_wo_na <- x[!is.na(x[,1]),1:2]
      ts_segments <- index.segments()
    }


    if(!is.vector(order)) stop("order is not a vector")
    switch(control$grad,
           analytical = (agrad <- TRUE),
           numerical = (agrad <- FALSE))
    if(is.null(series)) series <- deparse(substitute(x))
    ists <- is.ts(x_wo_na)
    x_wo_na <- as.ts(x_wo_na[,1])
    xfreq <- frequency(x_wo_na)

    if(ists) xtsp <- tsp(x_wo_na)
    x_wo_na <- as.matrix(x_wo_na)

    n <- nrow(x)
    n_wo_na <- nrow(x_wo_na)
    n_segments <- nrow(ts_segments)
    n_max_segment <- max(ts_segments$length)

    e <- double(n)
    e_wo_na <- double(n_wo_na)

    ncoef <- order[1]+order[2]+1
    hess <- matrix(0.0, ncoef, ncoef)
    small <- 0.05
    coef <- control$start
    if(is.null(coef))
        coef <- c(var(x_wo_na)*(1.0-small*(ncoef-1)),rep.int(small,ncoef-1))
    if(!is.vector(coef)) stop("coef is not a vector")
    if(ncoef != length(coef)) stop("incorrect length of coef")
    nlikeli <- 1.0e+10

    fit <- .C("fit_grouped_garch",
              as.vector(x_wo_na, mode = "double"),
              as.integer(n_wo_na),
              as.integer(ts_segments[,3]),
              as.integer(ts_segments[,4]),
              as.integer(ts_segments[,5]),
              as.integer(n_segments),
              coef = as.vector(coef, mode = "double"),
              as.integer(order[1]),
              as.integer(order[2]),
              as.integer(control$maxiter),
	            as.double(control$abstol),
	            as.double(control$reltol),
	            as.double(control$xtol),
	            as.double(control$falsetol),
              nlikeli = as.double(nlikeli),
              as.integer(agrad),
              as.integer(control$trace),
              PACKAGE="groupedtseries")

    pred <- .C("pred_grouped_garch",
               as.vector(x_wo_na, mode = "double"),
               e_wo_na = as.vector(e_wo_na, mode = "double"),
               as.integer(ts_segments[,3]),
               as.integer(ts_segments[,4]),
               as.integer(ts_segments[,5]),
               as.integer(n_segments),
               as.integer(n_wo_na),
               as.vector(fit$coef, mode = "double"),
               as.integer(order[1]),
               as.integer(order[2]),
               as.integer(FALSE),
               PACKAGE = "groupedtseries")

    com.hess <- .C("ophess_grouped_garch",
                   as.vector(x_wo_na, mode = "double"),
                   as.integer(ts_segments[,3]),
                   as.integer(ts_segments[,4]),
                   as.integer(ts_segments[,5]),
                   as.integer(n_segments),
                   as.integer(n_wo_na),
                   as.vector(fit$coef, mode = "double"),
                   hess = as.matrix(hess),
                   as.integer(order[1]),
                   as.integer(order[2]),
                   PACKAGE="groupedtseries")
    rank <- do.call("qr", c(list(x = com.hess$hess), control$qr))$rank
    if(rank != ncoef) {
	vc <- matrix(NA, nrow = ncoef, ncol = ncoef)
        warning("singular information")
    }
    else
        vc <- solve(com.hess$hess)

    temp_sigt <- sqrt(pred$e_wo_na)
    sigt <- as.double(rep(NA, n))

    for (i in 1:n_segments) {
      if(ts_segments[i,5] > max(order[1],order[2])) {
        temp_sigt[ts_segments[i,3]:(ts_segments[i,3]+max(order[1],order[2]) -1)] <- rep.int(NA, max(order[1],order[2]))
        sigt[ts_segments[i,1]:ts_segments[i,2]] <- temp_sigt[ts_segments[i,3]:ts_segments[i,4]]
      }
    }

    f <- cbind(sigt,-sigt)
    colnames(f) <- c("sigt","-sigt")

    e <- as.vector(x[,1])/sigt

    if(ists) {
        attr(e, "tsp") <-  attr(f, "tsp") <- xtsp
        attr(e, "class") <- attr(f, "class") <- "ts"
    }

    names(order) <- c("p","q")

    coef <- fit$coef

    nam.coef <- "a0"

    if(order[2] > 0)
        nam.coef <- c(nam.coef, paste("a", seq(order[2]), sep = ""))
    if(order[1] > 0)
        nam.coef <- c(nam.coef, paste("b", seq(order[1]), sep = ""))
    names(coef) <- nam.coef
    colnames(vc) <- rownames(vc) <- nam.coef
    grouped_garch <- list(order = order,
                  coef = coef,
                  n.likeli = fit$nlikeli,
                  n.used = n_wo_na,
                  residuals = e,
                  fitted.values = f,
                  series = series,
                  frequency = xfreq,
                  call = match.call(),
		  vcov = vc)
    class(grouped_garch) <- "garch"
    return(grouped_garch)
}

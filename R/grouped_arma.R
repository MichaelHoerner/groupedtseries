#' Fit grouped ARMA model to time series
#'
#' Fit an ARMA model with the order c(p,q) to a group of univariate time series
#' by finding coefficients that serve best the whole group.
#'
#' @param x a two-dimensional numeric vector of a grouped time series with accompanying group index in second dimension.
#' One-dimensional if input shall be only a single time series.
#' @param order a vector c(p,q) that indicates the order p of the AR part and q of the MA part
#' @param quick.estimation logical value (TRUE, FALSE) indicating whether a quick estimation
#' mode shall be used
#' @param lag a list with components \code{ar} and \code{ma}.
#' Each component is an integer vector, specifying the AR and MA lags that are
#' included in the model. If both, \code{order} and \code{lag}, are
#' given, only the specification from \code{lag} is used.
#' @param coef If given this numeric vector is used as the initial estimate
#' of the ARMA coefficients. The preliminary estimator suggested in
#' Hannan and Rissanen (1982) is used for the default initialization.
#' @param include.intercept Should the model contain an intercept?
#' @param series name for the series. Defaults to \code{deparse(substitute(x))}.
#' @param qr.tol the \code{tol} argument for \code{\link{qr}} when computing
#' the asymptotic standard errors of \code{coef}.
#' @param \dots additional arguments for \code{\link{optim}} when fitting
#' the model.
#'
#'
#' @return Object of class \code{arma}.
#'    \item{lag}{the lag specification of the fitted model.}
#'    \item{coef}{estimated ARMA coefficients for the fitted model across all grouped time series.}
#'    \item{css}{the conditional sum-of-squared errors.}
#'    \item{n.used}{the number of observations of \code{x}.}
#'    \item{residuals}{the series of residuals.}
#'    \item{fitted.values}{the fitted series.}
#'    \item{series}{the name of the series \code{x}.}
#'    \item{frequency}{the frequency of the series \code{x}.}
#'    \item{call}{the call of the \code{arma} function.}
#'    \item{vcov}{estimate of the asymptotic-theory covariance matrix for the
#'      coefficient estimates.}
#'    \item{convergence}{The \code{convergence} integer code from
#'        \code{\link{optim}}.}
#'    \item{include.intercept}{Does the model contain an intercept?}
#'
#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export
#'
grouped_arma <- function (x, order = c(1, 1), lag = NULL, coef = NULL,
                                        include.intercept = TRUE, series = NULL, qr.tol = 1e-07,
                                        quick.estimation = FALSE, ...)
{
  #create sequence from 1:N
  seqN <- function(N) {
    if (0 == length(N))
      NULL
    else if (N <= 0)
      NULL
    else seq(N)
  }

  #target function
  err <- function(coef) {
    u <- double(n)
    r_squared_term <- 0

    for (i in 1:n_segments){

      u_temp <- double(ts_segments$length[i])

      if (ts_segments$length[i] > max.order) {

        u_temp[seqN(max.order)] <- 0
        u_temp <- .C("arma", as.vector(x[ts_segments$index_start[i]:ts_segments$index_end[i]],
                                       mode = "double"), u = as.vector(u_temp),
                     as.vector(coef, mode = "double"), as.integer(lag$ar),
                     as.integer(lag$ma), as.integer(ar.l), as.integer(ma.l),
                     as.integer(max.order), as.integer(ts_segments$length[i]),
                     as.integer(include.intercept),
                     PACKAGE = "groupedtseries")$u

      } else {

        u_temp[1:length(u_temp)] <- 0

      }

      u[ts_segments$index_start[i]:ts_segments$index_end[i]] <- u_temp

      #Alternative - return minimum compound R2
      if (length(unique(x[(ts_segments$index_start[i]+max.order):ts_segments$index_end[i]]))==1) {
        r_squared_term <-
          r_squared_term + (ts_segments$length[i]-max.order)*sum((u_temp)^2) /
          sum((x[(ts_segments$index_start[i]+max.order):ts_segments$index_end[i]] - mean(na.remove(x)))^2)

        # print(paste0("r_squared_term above: ", r_squared_term))

      } else {
        r_squared_term <-
          r_squared_term + (ts_segments$length[i]-max.order)*sum((u_temp)^2) /
          sum((x[(ts_segments$index_start[i]+max.order):ts_segments$index_end[i]]-mean(x[(ts_segments$index_start[i]+max.order):ts_segments$index_end[i]]))^2)

        # print(paste0("r_squared_term below: ", r_squared_term))
        # print(paste0("mean : ", mean(x[(ts_segments$index_start[i]+max.order):ts_segments$index_end[i]])))
      }

    }
    #return of squared error - to be minimized in optimization
    #return(sum(u^2))

    #Alternative - return minimum compound R2
    return(r_squared_term)

  }

  #return residual error
  resid <- function(coef) {
    u <- double(n)
    u[1:length(u)] <- NA
    for (i in 1:n_segments){

      u_temp <- double(ts_segments$length[i])

      if (ts_segments$length[i] > max.order) {

        #print(paste0("x ", x[ts_segments$index_start[i]:ts_segments$index_end[i]]))
        #print(paste0("coef ", coef))

        u_temp[seqN(max.order)] <- 0
        u_temp <- .C("arma", as.vector(x[ts_segments$index_start[i]:ts_segments$index_end[i]],
                                       mode = "double"), u = as.vector(u_temp),
                     as.vector(coef, mode = "double"), as.integer(lag$ar),
                     as.integer(lag$ma), as.integer(ar.l), as.integer(ma.l),
                     as.integer(max.order), as.integer(ts_segments$length[i]),
                     as.integer(include.intercept),
                     PACKAGE = "groupedtseries")$u

        u_temp[seqN(max.order)] <- NA

      } else {

        u_temp[1:length(u_temp)] <- NA

      }

      u[ts_segments$index_start[i]:ts_segments$index_end[i]] <- u_temp

    }

    return(u)
  }

  #initialize ARMA coefficients - if quick.estimatin == TRUE, this will serve as only
  #calculation of coefficients
  arma.init <- function() {

    x_wo_na <- na.remove(as.numeric(x))

    k <- max.order
    e <- na.omit(drop(ar.ols(x_wo_na, order.max = k, aic = FALSE,
                             demean = FALSE, intercept = include.intercept)$resid))
    ee <- embed(as.vector(e), max.order + 1)
    xx <- embed(x_wo_na[-(1:k)], max.order + 1)
    if (include.intercept == TRUE) {
      if (is.null(lag$ar))
        coef <- lm(xx[, 1] ~ ee[, lag$ma + 1])$coef
      else if (is.null(lag$ma))
        coef <- lm(xx[, 1] ~ xx[, lag$ar + 1])$coef
      else coef <- lm(xx[, 1] ~ xx[, lag$ar + 1] + ee[,
                                                      lag$ma + 1])$coef
      coef <- c(coef[-1], coef[1])
    } else {
      if (is.null(lag$ar))
        coef <- lm(xx[, 1] ~ ee[, lag$ma + 1] - 1)$coef
      else if (is.null(lag$ma))
        coef <- lm(xx[, 1] ~ xx[, lag$ar + 1] - 1)$coef
      else coef <- lm(xx[, 1] ~ xx[, lag$ar + 1] + ee[,
                                                      lag$ma + 1] - 1)$coef
    }
    return(coef)
  }


  if (!is.null(order) && !is.null(lag))
    warning("order is ignored")
  if (is.null(order) && is.null(lag))
    stop("order or lag must be given")
  if (is.null(lag) && !is.null(order))
    lag <- list(ar = seqN(order[1]), ma = seqN(order[2]))

  lag$ar <- unique(lag$ar)
  lag$ma <- unique(lag$ma)
  max.order <- max(unlist(lag), 0)
  ar.l <- length(lag$ar)
  ma.l <- length(lag$ma)

  n <- nrow(x)

  if (NCOL(x) > 2) {
    stop("x is not a vector or univariate time series with indizes.
         There is at least a third column.")
  } else if (NCOL(x) == 2 && (any(is.na(x[,2])) || any(x[,2] == 0))) {
    stop("0 or NA are no valid entries for time series index.")
  } else if (NCOL(x) == 1){
    ts_segments <- data.frame(index_start = 1, index_end = n, length = n)
  } else {
    ts_segments <- index_segments(x, max.order)
  }

  if (is.null(series))
    series <- deparse(substitute(x))

  ists <- is.ts(x)
  x <- as.ts(x[,1])
  xfreq <- frequency(x)

  if (ists)
    xtsp <- tsp(x)

  n_wo_na <- sum(ts_segments$length)
  n_segments <- nrow(ts_segments)
  n_max_segment <- max(ts_segments$length)

  if (!is.null(unlist(lag)))
    if ((min(unlist(lag)) < 1) || (max(unlist(lag)) > (n -
                                                       1)))
      stop("invalid lag")

  ncoef <- length(unlist(lag)) + as.numeric(include.intercept)

  if (is.null(coef)) {
    if (!is.null(unlist(lag))) {
      coef <- arma.init()

    #prevent NA values in coef
      coef[which(is.na(coef))] <- 0
    } else {
      coef <- 0
    }
  }

  if (length(coef) != ncoef)
    stop("invalid coef")

  if (quick.estimation == FALSE) {
    md <- optim(coef, err, gr = NULL, hessian = TRUE, ...)
    coef <- md$par
    rank <- qr(md$hessian, qr.tol)$rank
    if (rank != ncoef) {
      vc <- matrix(NA, nrow = ncoef, ncol = ncoef)
      warning("singular Hessian")
    }
    else {
      vc <- 2 * md$value/n * solve(md$hessian)
      if (any(diag(vc) < 0))
        warning("Hessian negative-semidefinite")
    }
  } else {
    vc <- matrix(NA, nrow = ncoef, ncol = ncoef)
  }

  e <- resid(coef)

  f <- x - e

  if (ists) {
    attr(e, "tsp") <- xtsp
    attr(e, "class") <- "ts"
    attr(f, "tsp") <- xtsp
    attr(f, "class") <- "ts"
  }
  nam.ar <- if (!is.null(lag$ar))
    paste("ar", lag$ar, sep = "")
  else NULL
  nam.ma <- if (!is.null(lag$ma))
    paste("ma", lag$ma, sep = "")
  else NULL
  nam.int <- if (include.intercept)
    "intercept"
  else NULL
  nam.coef <- c(nam.ar, nam.ma, nam.int)
  names(coef) <- nam.coef
  colnames(vc) <- rownames(vc) <- nam.coef

  if (quick.estimation == FALSE) {
    group_arma <- list(coef = coef, n.used = n_wo_na, residuals = e,
                       fitted.values = f, series = series, frequency = xfreq,
                       call = match.call(), vcov = vc, lag = lag, include.intercept = include.intercept,
                       css = md$value, convergence = md$convergence)
  } else {
    group_arma <- list(coef = coef, n.used = n_wo_na, residuals = e,
                       fitted.values = f, series = series, frequency = xfreq,
                       call = match.call(), vcov = vc, lag = lag, include.intercept = include.intercept,
                       css = 1, convergence = NA)
  }

  class(group_arma) <- "arma"
  return(group_arma)
}

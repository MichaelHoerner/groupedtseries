grouped_garch <-
function (x, order = c(1, 1), series = NULL, control = grouped_garch.control(...), ...)
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

    ts_segments <- data.frame(cbind(index_start, index_end, index_end - index_start + 1))
    colnames(ts_segments) <- c("index_start", "index_end", "length")

    return(ts_segments)
  }

    if (NCOL(x) > 2) {
      stop("x is not a vector or univariate time series with indizes.
               There is at least a third column.")
    } else if (NCOL(x) == 2 && (any(is.na(x[,2])) || any(x[,2] == 0))) {
      stop("0 or NA are no valid entries for time series index.")
    } else if (NCOL(x) == 1){
      index_x <- c(rep(1, length(x)))
      x <- cbind(x, index_x)
      ts_segments <- index.segments()
    } else {
      ts_segments <- index.segments()
    }



    if(!is.vector(order)) stop("order is not a vector")
    switch(control$grad,
           analytical = (agrad <- TRUE),
           numerical = (agrad <- FALSE))
    if(is.null(series)) series <- deparse(substitute(x))
    ists <- is.ts(x)
    x <- as.ts(x[,1])
    xfreq <- frequency(x)

    if(ists) xtsp <- tsp(x)
    x <- as.matrix(x)

    n <- nrow(x)
    n_wo_na <- sum(ts_segments$length)
    n_segments <- nrow(ts_segments)
    n_max_segment <- max(ts_segments$length)

    e <- double(n)
    e_wo_na <- double(n_wo_na)

    ncoef <- order[1]+order[2]+1
    hess <- matrix(0.0, ncoef, ncoef)
    small <- 0.05
    coef <- control$start
    if(is.null(coef))
        coef <- c(var(na.omit(x))*(1.0-small*(ncoef-1)),rep.int(small,ncoef-1))
    if(!is.vector(coef)) stop("coef is not a vector")
    if(ncoef != length(coef)) stop("incorrect length of coef")
    nlikeli <- 1.0e+10
    fit <- .C("fit_grouped_garch",
              as.vector(x, mode = "double"),
              as.integer(n),
              as.integer(ts_segments[,1]),
              as.integer(ts_segments[,2]),
              as.integer(ts_segments[,3]),
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
               as.vector(x, mode = "double"),
               e = as.vector(e, mode = "double"),
               as.integer(ts_segments[,1]),
               as.integer(ts_segments[,2]),
               as.integer(ts_segments[,3]),
               as.integer(n_segments),
               as.integer(n),
               as.vector(fit$coef, mode = "double"),
               as.integer(order[1]),
               as.integer(order[2]),
               as.integer(FALSE),
               PACKAGE = "groupedtseries")
    # com.hess <- .C("ophess_grouped_garch",
    #                as.vector(na.omit(x), mode = "double"),
    #                as.integer(n_wo_na),
    #                as.vector(fit$coef, mode = "double"),
    #                hess = as.matrix(hess),
    #                as.integer(order[1]),
    #                as.integer(order[2]),
    #                PACKAGE="groupedtseries")
    rank <- 0 #do.call("qr", c(list(x = com.hess$hess), control$qr))$rank
    if(rank != ncoef) {
	vc <- matrix(NA, nrow = ncoef, ncol = ncoef)
        warning("singular information")
    }
    else
        vc <- solve(com.hess$hess)

    sigt <- sqrt(pred$e)
    sigt[1:max(order[1],order[2])] <- rep.int(NA, max(order[1],order[2]))

    f <- cbind(sigt,-sigt)
    colnames(f) <- c("sigt","-sigt")
    e_wo_na <- as.vector(x)/sigt
    if(ists) {
        attr(e_wo_na, "tsp") <-  attr(f, "tsp") <- xtsp
        attr(e_wo_na, "class") <- attr(f, "class") <- "ts"
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
                  residuals = e_wo_na,
                  fitted.values = f,
                  series = series,
                  frequency = xfreq,
                  call = match.call(),
		  vcov = vc)
    class(grouped_garch) <- "grouped_garch"
    return(grouped_garch)
}

#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export



index_segments <- function(x) {

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

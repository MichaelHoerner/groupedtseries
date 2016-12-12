#' @useDynLib groupedtseries
#' @importFrom Rcpp sourceCpp
#' @importFrom stats is.ts as.ts
#' @export

create_X_AR <- function(y,taille,AR_lags,MA_lags, ts_segments) {

number_ts_segments <- nrow(ts_segments)

X_AR <- matrix(0,taille-AR_lags*number_ts_segments,1+AR_lags+MA_lags)
X_AR[,1] <- matrix(1,taille-AR_lags*number_ts_segments,1)

y_AR <- vector(length = taille-AR_lags*number_ts_segments, mode = "double")

ts_segments <- cbind(ts_segments, matrix(0, number_ts_segments, 3))
colnames(ts_segments) <- c("index_start", "index_end", "index_start_wo_na",
                           "index_end_wo_na", "length", "index_AR_start", "index_AR_end", "length_AR")


for (seg in 1:number_ts_segments){
  if(seg > 1){
    ts_segments$index_AR_start[seg] <- ts_segments$index_start_wo_na[seg]-AR_lags*(seg-1)
  } else {
    ts_segments$index_AR_start[seg] <- ts_segments$index_start_wo_na[seg]
  }
  ts_segments$index_AR_end[seg] <- ts_segments$index_end_wo_na[seg]-AR_lags*(seg)
  ts_segments$length_AR[seg] <- ts_segments$length[seg]-AR_lags
  y_AR[ts_segments$index_AR_start[seg]:ts_segments$index_AR_end[seg]] <-
    y[(ts_segments$index_start[seg]+AR_lags):ts_segments$index_end[seg],1]

  for (  i in 1 : AR_lags ) {
    X_AR[ts_segments$index_AR_start[seg]:ts_segments$index_AR_end[seg],(AR_lags+2-i)] <-
      y[(ts_segments$index_start[seg]+i-1):(ts_segments$index_end[seg]-AR_lags+i-1),1]
  }
}

 #

returnvalues <- list(y_AR, X_AR, ts_segments)

return(returnvalues)

}


#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

#' print stdbscan class
#'
#' @export
#' @param x `stdbscan` 関数の結果オブジェクトを指定します。
#' @param ... `print` 関数で使用できる引数を指定します。
print.stdbscan <- function(x, ...) {
  num_d <- if (x$metric$neighbortype == "random") {
    nrow(x$results$geo)
  } else {
    nrow(x$results$geo) * length(x$results$time)
  }
  message("ST-DBSCAN clustering for ", num_d, " objects, ",
          length(x$results$time), " time length.")
  stdbscan_call(x)
}

#' print hstdbscan class
#'
#' @export
#' @param x `hstdbscan` 関数の結果オブジェクトを指定します。
#' @param ... `print` 関数で使用できる引数を指定します。
print.hstdbscan <- function(x, ...) {
  num_d <- if (x$metric$neighbortype == "random") {
    nrow(x$results$geo)
  } else {
    nrow(x$results$geo) * length(x$results$time)
  }
  message("Hierarchical ST-DBSCAN clustering for ", num_d, " objects, ",
          length(x$results$time), " time length.")
  stdbscan_call(x)
  message("\n\n", "Can use the `cutclust` function to split it into `k` clusters")
}


stdbscan_call <- function(x) {
  message("Parameters: ",
          "eps1 = ", x$eps1, ", eps2 = ", x$eps2, ", minPts = ", x$minPts)
  message("Using ", x$metric$metric, " distances, neighbor's metric = ", x$metric$neighbortype,
          ", ST-DBSCAN type = ", x$metric$dbscantype)
  clst_num <- unique(x$cluster)
  noise_num <- grep("Noise", x$cluster)
  if (length(noise_num) == 0) {
    cn <- length(clst_num)
    nn <- 0
  } else {
    cn <- length(clst_num) - 1
    nn <- length(noise_num)
  }  
  message("The clustering contains ", cn, " cluster(s) and ", nn, " noise points.")
  message("D (Δeps): ", paste("\n  ", names(x$delta_eps), paste0("(", x$delta_eps, ")"), collapse = ", "))
}

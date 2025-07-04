#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL


#' plot stclust data
#'
#' @export
#' @param x `stclust` 関数の結果オブジェクトを指定します。
#' @param ... `plot` 関数で指定できる引数を指定します。
plot.stclust <- function(x, ...) {
  plot(x$geo, ...)
}


#' plot cluster dendrogram
#'
#' @export
#' @param x `hstdbscan` 関数の結果オブジェクトを指定します。
#' @param ... `plot` 関数で指定できる引数を指定します。
plot.hstdbscan <- function(x, ...) {
  plot(x$hc, ...)
}


#' plot ST-DBSCAN
#'
#' `stdbscan` 関数で計算した結果オブジェクトを表示します。
#' 
#' @export
#' @param x `stdbscan` 関数の結果オブジェクトを指定します。
#' @param ... `plot` 関数で指定できる引数を指定します。
plot.stdbscan <- function(x, ...) {
  d <- if (x$metric$neighbortype == "spatial") {
    do.call(cbind, lapply(x$results$value, function(m) {
      as.integer(m$cluster)
    }))
  } else {
    res <- subset(x$results$value,
                  select = c(time, id, cluster))
    res$cluster <- as.integer(res$cluster)
    reshape(res, 
            direction = "wide",
            v.names = "cluster",
            timevar = "time",
            idvar = "id")[-1]
  }
  d <- as.matrix(t(d))
  image(x$results$time,
        x$results$geo$id, 
        d,
        xlab = "Time",
        ylab = "ID (geospatial)",
        main = "ST-DBSCAN Clustering")
}

#' Draw rectangle using Hierarchical ST-DBSCAN
#'
#' @export
#' @param x `hstdbscan` 関数の結果オブジェクトを指定します。
#' @param k クラスター分割数を指定します。
#' @param ... その他、`rect.hclust` で指定できる引数を設定します。
rect_hstdbscan <- function(x, k, ...) {
  rect.hclust(x$hc, k, ...)
}

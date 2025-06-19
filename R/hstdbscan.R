#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom cluster daisy
## usethis namespace: end
NULL

#' calc Hierarchical ST-DBSCAN clustering.
#'
#' 階層的 (Hierarchical) ST-DBSCAN法でクラスタリングを行います。
#'
#' @export
#' @param x `stdbscan` 関数の `x` 引数の説明を参照してください。
#' @param ... `stdbscan` 関数の `x` 引数以外の引数を指定します。
#' @param dist_option 階層的クラスタリングの距離の計算条件を指定します。ここでは `cluster::daisy` 関数の引数を `dist_option` 引数に `list` 形式で指定します。
#' @param hclust_option クラスタリングの条件を指定します。 `base::hclust` 関数で利用できる引数を `list` 形式で指定します。
#' @return
#' `stdbscan` の結果に加え、以下の項目が表含まれています。
#' * dist: 距離の計算結果テーブルが格納されています。
#' * hc: クラスタリングのツリー情報が含まれています。
#' @examples
#' x <- seq(130, 140, by = 1)
#' y <- seq(30, 40, by = 1)
#' t <- as.POSIXct("2024-01-01 00:00:00", tz = "JST") + 3600 * seq(0, 23, by = 6)
#' geo <- sf::st_as_sf(expand.grid(x, y), coords = c("Var1", "Var2"), crs = 4326)
#' 
#' D <- abs(runif(nrow(geo) * length(t))) * 100
#' clust <- stdbscan(x = geo, time = t, eps1 = 144, eps2 = 3600 * 6, minPts = 6,
#'                   vals = list(list(D = D,
#'                                    delta_eps = 20)),
#'                   metric = "geo", neighbortype = "spatial", dbscantype = "grid")
#' cutclust(clust, 4)
hstdbscan <- function(x,
                      ...,
                      dist_option = list(metric = "gower"),
                      hclust_option = list(method = "ward.D2")) {
  clust <- stdbscan(x, ...)
  d <- if (clust$metric$neighbortype == "spatial") {
    reshape(subset(do.call(rbind, lapply(clust$results$value, function(m) m)),
                   select = c(id, time, cluster)),
            direction = "wide",
            v.names = "cluster",
            timevar = "time",
            idvar = "id")
  } else {
    subset(clust$results$value, select = c(id, cluster))
  }
  rownames(d) <- seq(1, nrow(clust$results$geo), by = 1)
  def <- do.call(cluster::daisy, c(list(x = d[-1]), dist_option))
  hc <- do.call(hclust, c(list(d = def), hclust_option))
  ret <- c(clust,
           list(hc = hc,
                dist = def))
  class(ret) <- "hstdbscan"
  return(ret)
}


#' `hstdbscan` 関数で計算したクラスターの結果を分割します。
#'
#' @export
#' @param x `hstdbscan` の結果オブジェクトを指定します。
#' @param k クラスターの分割数を整数で指定します。
#' @param ... その他、 `cutree` 関数で使用できる引数を指定することができます。
#' @return
#' 階層的クラスタリング結果のラベル、位置情報が含まれる data.frame が返ります。
#' @examples
#' x <- seq(130, 140, by = 1)
#' y <- seq(30, 40, by = 1)
#' t <- as.POSIXct("2024-01-01 00:00:00", tz = "JST") + 3600 * seq(0, 23, by = 6)
#' geo <- sf::st_as_sf(expand.grid(x, y), coords = c("Var1", "Var2"), crs = 4326)
#' 
#' D <- abs(runif(nrow(geo) * length(t))) * 100
#' clust <- stdbscan(x = geo, time = t, eps1 = 144, eps2 = 3600 * 6, minPts = 6,
#'                   vals = list(list(D = D,
#'                                    delta_eps = 20)),
#'                   metric = "geo", neighbortype = "spatial", dbscantype = "grid")
#' cutclust(clust, 4)
cutclust <- function(x, k, ...) {
  hc <- x$hc
  new <- cutree(tree = hc, k = k, ...)
  cbind(cluster = as.factor(new),
        x$results$geo)
}

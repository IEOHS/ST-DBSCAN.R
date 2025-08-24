#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import sf
#' @import ggplot2
#' @importFrom tidyr pivot_longer
## usethis namespace: end
NULL

#' plot for time-series data about sparseDFM model, and check models.
#'
#' @export
#' @param x stclust クラスのデータオブジェクトを指定します。
#' @param time_label 時系列情報を指定します。
#' @param ... geom_path で指定できる引数値を指定します。
validPlot <- function(x, time_label = seq(1, ncol(x$data)), ...) {
  p1 <- spClustPlot(x, ...)
  p2 <- tsFactorPlot(x, time_label)
  gridExtra::grid.arrange(p1, p2, ncol = 2, widths = c(1, 1.3))
}

spClustPlot <- function(x, ...) {
  p1 <- ggplot(subset(x$geo, select = cluster)) +
    geom_sf(aes(fill = cluster, col = cluster), ...)
  p2 <- tidyr::pivot_longer(subset(x$geo, select = -cluster),
                            cols = -geometry) |>
    ggplot() +
    geom_sf(aes(fill = value, col = value)) +
    facet_wrap(name ~ ., ncol = 1)
  gridExtra::grid.arrange(p1, p2, ncol = 1, heights = c(1, 2))
}

tsFactorPlot <- function(x, time_label = seq(1, ncol(x$data)), ...) {
  force(time_label)
  ptable <- do.call(rbind, lapply(x$ts_model, function(m) {
    obs <- scale(m$ts_model$data$X)
    pred <- scale(m$ts_model$state$factors)
    tb <- tidyr::pivot_longer(as.data.frame(cbind(date = time_label,
                                                  pred, obs)),
                              cols = -date)
    tb$col <- ifelse(grepl("pred|F\\d", tb$name), "black", "gray")
    tb$cluster <- m$cluster
    tb
  }))

  ggplot(ptable) +
    geom_path(aes(x = date, y = value, group = name), col = ptable$col, ...) +
    facet_grid(cluster ~ .) +
    labs(y = "normalized value (per cluster)", x = "Elapsed Time")
}


#' plot stclust data
#'
#' @export
#' @param x `stclust` 関数の結果オブジェクトを指定します。
#' @param bgmap 背景用の地図を sf クラスのデータで指定します。
#' @param ... `plot` 関数で指定できる引数を指定します。
plot.stclust <- function(x, bgmap = NULL, ...) {
  map <- if (!is.null(bgmap)) {
    rbind(x$geo, sf::st_sf(cluster = NA, sp_check = NA, ts_check = NA, geometry = sf::st_geometry(bgmap)))
  } else {
    x$geo
  }
  ncol_num <- if (length(levels(x$geo$cluster)) + 4 > 12) {
    3
  } else {
    floor((length(levels(x$geo$cluster)) + 4) / 4)
  }
  cols <- sf::sf.colors(n = max(4, length(levels(x$geo$cluster))), categorical = TRUE)
  legend_icon <- c("***", "**", "*", ".")
  plot(map, ...)
  legend("bottomright",
         legend = c(paste("Sign:", legend_icon),
                    paste("Cluster:", levels(x$geo$cluster))),
         col = c(cols[1:4], cols),
         pch = 15,
         title = "contour color",
         ncol = ncol_num)
}
# plot.stclust <- function(x, ...) {
#   plot(x$geo, ...)
# }


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

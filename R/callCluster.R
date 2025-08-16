#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom dtwclust tsclust
#' @import sf
## usethis namespace: end
NULL

#' run any clustering method.
#'
#' 複数のクラスタリング方法でクラスター解析を実施します。
#'
#' @export
#' @param x 行方向がユニークな位置、列方向が時系列となる2次元配列のデータを指定します。
#' @param geo 位置情報を `sf` または `sfc` 形式で指定します。
#' @param times 時間情報を指定します。
#' @param cuts クラスター分割数を整数で指定します。複数同時に指定する場合はベクトル (`c(1, 2, 3, ...)` のように指定してください。)
#' @param method クラスタリング方法を指定します。ウォード法 ("ward"), k-means法 ("kmeans"), k-shape法 ("kshape"), 動的時間伸縮法 ("dtw"), 階層的ST-DBSCAN法 ("hstdbscan") から一つ指定します。
#' @param ... クラスタリングの引数オプションを指定します。
#' @return
#' * cluster
#' cuts で指定した分割数 (k) 毎にクラスタリングラベルが作成され、位置情報とマージされて sf 形式のデータで格納されます。
#' * model
#' クラスタリングを実行した結果が格納されています。 `type = "partitonal"` の場合は、cuts 毎にモデルが作成され、リスト形式で格納されています。
#' * type
#' クラスタリングの方法が階層的 ("hierarchical") か、非階層的 ("partitional") か示されます。
#' * method
#' クラスタリング方法で指定した method 引数の値が格納されます。
#' @examples
#' x <- seq(130, 140, by = 1)
#' y <- seq(30, 40, by = 1)
#' t <- as.POSIXct("2024-01-01 00:00:00", tz = "JST") + 3600 * seq(0, 23, by = 6)
#' geo <- sf::st_as_sf(expand.grid(x, y), coords = c("Var1", "Var2"), crs = 4326)
#' D <- abs(runif(nrow(geo) * length(t))) * 100
#' D <- Reduce(cbind, split(D, do.call(c, Map(rep, t, nrow(geo)))), init = NULL)
#' clust <- geostclust(x = D, geo = geo, times = t, method = "ward")
#' clust <- geostclust(x = D, geo = geo, times = t, method = "kmeans")
#' clust <- geostclust(x = D, geo = geo, times = t, method = "kshape")
#' clust <- geostclust(x = D, geo = geo, times = t, method = "dtw")
#' clust <- geostclust(geo = geo, times = t, method = "hstdbscan", eps1 = 144, eps2 = 3600 * 6, minPts = 6,
#'                     vals = list(list(D = D,
#'                                      delta_eps = 20)))
geostclust <- function(x = double(),
                       geo,
                       times,
                       cuts = c(4, 8, 16, 32), 
                       method = c("ward", "kmeans", "kshape", "dtw", "hstdbscan", "others"),
                       ...) {
  
  method <- match.arg(method)
  if (!any(class(geo) %in% c("sf", "sfc"))) {
    stop("`geo` object must be 'sf' or 'sfc' class.")
  }
  geo <- sf::st_geometry(geo)
  if (method != "hstdbscan") {
    d <- matrix(x, nrow = length(geo), ncol = length(times))
  }
  cluster_label <- paste(method, cuts, "cluster", sep = "_")
  
  if (method == "ward") {
    model <- hclust(dist(d), method = "ward.D2", ...)
    clust <- list(label = Map(assign, cluster_label, Map(cutree, list(model), k = cuts)),
                  model = model)
  } else if (method == "kmeans") {
    model <- Map(assign, cluster_label, lapply(cuts, function(n) {
      kmeans(d, centers = n, ...)
    }))
    clust <- list(label = Map(function(m) m$cluster, model),
                  model = model)
  } else if (method == "kshape") {
    model <- Map(assign, cluster_label, lapply(cuts, function(n) {
      dtwclust::tsclust(series = d,
                        type = "partitional",
                        k = n,
                        distance = "sbd",
                        centroid = "shape", ...)
    }))
    clust <- list(label = Map(function(m) m@cluster, model),
                  model = model)
  } else if (method == "dtw") {
    model <- dtwclust::tsclust(series = d,
                               type = "hierarchical",
                               distance = "dtw_basic",
                               ...)
    clust <- list(label = Map(assign, cluster_label, Map(cutree, list(model), k = cuts)),
                  model = model)
  } else if (method == "hstdbscan") {
    model <- hstdbscan(x = geo, time = times, metric = "geo", neighbortype = "spatial", dbscantype = "grid", ...)
    clust <- list(label = Map(assign,
                              cluster_label,
                              lapply(cuts, function(n) {
                                cutclust(model, n)$cluster
                              })),
                  model = model)
  }

  list(cluster = sf::st_as_sf(do.call(data.frame, clust$label), geometry = sf::st_geometry(geo)),
       model = clust$model,
       type = if (method %in% c("ward", "dtw", "hstdbscan")) "hierarchical" else "partitional",
       method = method)
  
}

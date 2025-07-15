#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom sf st_as_sf
#' @importFrom sf st_coordinates
#' @importFrom sf st_geometry
#' @importFrom sf st_sf
#' @importFrom spdep dnearneigh
## usethis namespace: end
NULL


#' calc ST-DBSCAN clustering.
#'
#' ST-DBSCAN法でクラスタリングを行います。
#' article{birant2007st,
#'  title={ST-DBSCAN: An algorithm for clustering spatial--temporal data},
#'  author={Birant, Derya and Kut, Alp},
#'  journal={Data \& knowledge engineering},
#'  volume={60},
#'  number={1},
#'  pages={208--221},
#'  year={2007},
#'  publisher={Elsevier}
#'}
#'
#' @export
#' @param x matrix[, c(x, y)], data.frame(x, y), sf, sfc を指定します。
#' matrixの場合は、2列の行列を指定し、1列目にはx軸 (long), 2列目にはy軸 (lat) の情報を指定します。
#' data.frameの場合は、`x` 列にx軸, `y` 列にy軸の情報を指定します。
#' sf及びsfcの場合はそのまま指定できますが、sfデータにgeometryのみが利用されます。
#' neighbortype = "random" の場合、'x' に 時間情報を含めることができます。その場合、列名は 'time' としてください。
#' @param time 時間軸情報をベクトルで指定します。数値、POSIXct等、任意のデータ型を指定します。
#' @param eps1 空間的な距離の閾値を指定します。これは `dbscan::dbscan` 関数で指定する `eps` 引数の値に該当します。
#' @param eps2 時間的な距離の閾値をしています。
#' `time` 引数を数値で指定した場合は、隣接時間となる値を数値で指定し、 `POSIXct` 形式で指定した場合は秒 (sec) で指定します。
#' @param minPts クラスターと認識する最低の隣接点数を整数で指定します。
#' @param vals 空間・時間軸に該当するデータを `list` 形式でしています。また、このデータのクラスターを判別するための閾値 (Δeps) も指定します。
#' データは次の例のように複数指定することができます
#' list(list(D = vector(), delta_eps = double()), list(D = vector(), delta_eps = double()), ...)
#' @param metric 空間的距離を算出するための方法を指定します。
#' * "euclidean": ユークリッド距離 (`sqrt((x1 - x1) ^ 2 + (y1 - y1) ^ 2 )`) で距離が計算されます。
#' * "geo": 地理空間上の座標 (long, lat) に基づき地点間の距離が計算されます。計算には `spdep::dnearneigh' が利用されます。
#' @param neighbortype 隣接点が空間的に同一か、完全にランダムな配置であるか指定します。
#' * "spatial": netCDFデータや定点データのように、位置が時間に応じて変わらないデータを使用する場合に指定します。
#' * "random": 地震発生源のような、時間と空間上でランダムなデータを利用する場合に指定します。
#' @param dbscantype クラスターの判定の際、netCDF形式で配布されている気象データのようなグリッドデータでは、常に一定数の隣接点が存在します。
#' このような場合、ST-DBSCANのもともとの計算アルゴリズムでは必ずクラスターが生成されるため、意図しない計算結果となることが想定されます。
#' `dbscantype` 引数では、"default" を指定すると従来のアルゴリズムを踏襲し、"grid" を指定すると隣接点と中心点との値差がΔepsとなる隣接点を抽出し、その数が `minPts` 以上であればクラスターと認識する計算が実行されます。
#' @param FUN 集計用の計算関数を指定します。
#' @return
#' * cluster: valsに指定したデータ数と対応するクラスタリングの結果ラベルが含まれています。
#' * eps1: `eps1` 引数の値
#' * eps2: `eps2` 引数の値
#' * minPts: `minPts` 引数の値
#' * metric: 距離の計算、クラスタリングの条件が含まれています。
#'   - metric: `metric` 引数の値
#'   - neighbortype: `neighbortype` 引数の値
#'   - dbscantype: `dbscantype` 引数の値
#' * neighborlist: データの隣接点リスト (`nb` class) が含まれています。
#' * delta_eps: `vals` 引数の `delta_eps` で指定した値
#' * results: クラスタリングの結果が位置情報、時間情報、クラスタリングラベルに分かれて格納されています。
#'   - geo: 位置情報
#'   - time: 時間情報
#'   - value: 位置情報に該当するid, 時間情報, vals引数で指定した各データ値, クラスターラベルがデータフレーム形式で格納されています。`neighbortype = "spatial"` の場合は、時間毎に list 化されて格納されています。
#' @examples
#' x <- seq(130, 140, by = 1)
#' y <- seq(30, 40, by = 1)
#' t <- as.POSIXct("2024-01-01 00:00:00", tz = "JST") + 3600 * seq(0, 23, by = 6)
#' geo <- sf::st_as_sf(expand.grid(x, y), coords = c("Var1", "Var2"), crs = 4326)
#' 
#' ## for ramdom
#' D <- abs(runif(nrow(geo))) * 100
#' clust <- stdbscan(x = cbind(geo, time = rep(t, nrow(geo))[1:nrow(geo)]),
#'                   eps1 = 144, eps2 = 3600 * 6, minPts = 6,
#'                   vals = list(list(D = D,
#'                                    delta_eps = 20)),
#'                   metric = "geo", neighbortype = "random", dbscantype = "default")
#' 
#' ## for spatial
#' D <- abs(runif(nrow(geo) * length(t))) * 100
#' clust <- stdbscan(x = geo, time = t, eps1 = 144, eps2 = 3600 * 6, minPts = 6,
#'                   vals = list(list(D = D,
#'                                    delta_eps = 20)),
#'                   metric = "geo", neighbortype = "spatial", dbscantype = "grid")
stdbscan <- function(x,
                     time = NULL, 
                     eps1 = double(),
                     eps2 = double(),
                     minPts = integer(),
                     vals = list(list(D = double(),
                                      delta_eps = double())),
                     metric = c("euclidean", "geo"),
                     neighbortype = c("spatial", "random"),
                     dbscantype = c("grid", "default"),
                     FUN = mean) {

  stopifnot(!is.null(x))
  stopifnot(length(eps1) > 0)
  stopifnot(length(eps2) > 0)
  stopifnot(length(minPts) > 0)
  metric <- match.arg(metric)
  neighbortype <- match.arg(neighbortype)
  dbscantype <- match.arg(dbscantype)

  if (neighbortype == "random" && is.null(time)) {
    time <- if (any(class(x) == "matrix")) {
      x[, "time"]
    } else {
      x$time
    } 
  }
  if (any(class(x) == "sf")) {
    x <- sf::st_geometry(x)
  }
  
  message("===== Start ST-DBSCAN method =====")
  message("\n1. Calculation Neighbor List")
  nb <- try(if (any(class(x) == "matrix")) {
    stnb(x = x[, "x"], y = x[, "y"], time = time, eps1 = eps1, eps2 = eps2,
         method = metric, neighbortype = neighbortype)
  } else if (any(class(x) == "data.frame")) {
    stnb(x = x$x, y = x$y, time = time, eps1 = eps1, eps2 = eps2,
         method = metric, neighbortype = neighbortype)
  } else {
    stnb(x = x, time = time, eps1 = eps1, eps2 = eps2,
         method = metric, neighbortype = neighbortype)
  })
  if (class(nb) == "try-error") {
    stop("Calc Error: Check your `x` data.")
  }

  message("\n2. Calculation Cluster")
  cluster <- try(as.factor(st_dbscan(nb = nb, vals = vals, type = dbscantype, minPts = minPts, FUN = FUN)))
  if (class(cluster) == "try-error") {
    stop("Calc Error: Check your `vals` data.")
  }

  message("\nCompleted.")
  obs_name <- paste0("Obs_", seq(1, length(vals), by = 1))
  raw_data <- local({
    d <- do.call(data.frame, lapply(vals, function(m) {
      m$D
    }))
    names(d) <- obs_name
    ret <- list()
    num <- if (any(class(x) == "sf") || any(class(x) == "data.frame") || any(class(x) == "matrix")) {
      nrow(x)
    } else {
      length(x)
    }
    if (neighbortype == "spatial") {
      d <- data.frame(id = seq(1, num, by = 1),
                      time = do.call(c, Map(rep, time, num)),
                      d)
      d$cluster <- cluster
      for (i in seq(1, length(time), by = 1)) {
        ret[[i]] <- d[seq(1, num, by = 1) + num * (i - 1), ]
      }
      names(ret) <- as.character(time)
    } else {
      ret <- data.frame(id = seq(1, num, by = 1),
                        time = time,
                        d,
                        cluster = cluster)
    }
    return(ret)
  })
  if (neighbortype == "random") {
    time <- sort(unique(time))
  }
  
  ret_table <- if (any(class(x) == "matrix")) {
    data.frame(id = seq(1, nrow(x), by = 1),
               x = x[, 1],
               y = x[, 2])
  } else if (any(class(x) == "data.frame")) {
    data.frame(id = seq(1, nrow(x), by = 1),
               x = x$x,
               y = x$y)
  } else {
    sf::st_sf(id = seq(1, length(x), by = 1),
              geometry = x)
  }
  ret <- list(cluster = cluster,
              eps1 = eps1,
              eps2 = eps2,
              minPts = minPts,
              metric = list(metric = metric,
                            neighbortype = neighbortype,
                            dbscantype = dbscantype),
              neighborlist = nb,
              delta_eps = local({
                d <- lapply(vals, function(m) {
                  m$delta_eps
                })
                names(d) <- obs_name
                return(d)
              }),
              results = list(geo = ret_table,
                             time = time,
                             value = raw_data))
  class(ret) <- "stdbscan"
  return(ret)
}

#' calc ST-DBSCAN clustering and return cluster label.
#'
#' ST-DBSCAN法でクラスタリングを行います。
#' https://www.sciencedirect.com/science/article/pii/S0169023X06000218
#' 
#' article{birant2007st,
#'  title={ST-DBSCAN: An algorithm for clustering spatial--temporal data},
#'  author={Birant, Derya and Kut, Alp},
#'  journal={Data \& knowledge engineering},
#'  volume={60},
#'  number={1},
#'  pages={208--221},
#'  year={2007},
#'  publisher={Elsevier}
#'}
#'
#' @param nb neighbor list
#' @param minPts minimum number of points
#' @param vals as delta_eps list
#' @param type
#' * "grid": For Grid data. Adjust nb size before calc.
#' * "default": For normal(random point) data. not Adjustment.
#' @param FUN 集計用の計算関数を指定します。
st_dbscan <- function(nb = NULL,
                      vals = list(list(D = NULL, ## val list
                                       delta_eps = 5)),
                      type = c("grid", "default"),
                      minPts = 10, FUN = mean) {
  type <- match.arg(type)
  stopifnot(is.list(vals))
  stopifnot(!is.null(nb))
  
  ## set cluster number
  cluster <- 0

  ## set cluster column
  ## label <- rep(NA, times = length(nb))
  label <- rep("Noise", times = length(nb))
  
  ## st-dbscan
  message("\nStart Clustering:  ", date())
  for (i in seq(1, length(nb))) {
    if (label[i] == "Noise") {
      neighbours <- if (type == "grid") {
        ## 中心点と隣接点との差がΔeps以下のデータのみ抽出
        rneighbors(nb = nb[[i]],
                   vals = vals,
                   mean_val = lapply(vals, function(m) {
                     m$D[i]
                   }),
                   label = NULL, mode = "all", FUN = FUN)
      } else if (type == "default") {
        nb[[i]]
      }

      ## minPtsより隣接数が少ない場合は "Noise" 判定
      if (minPts <= length(neighbours)) {
        ## up to cluster number
        cluster <- cluster + 1

        ## set cluster
        label[c(i, neighbours)] <- as.character(cluster)
        message("\tCreate Cluster: ", as.character(cluster))

        ## check node and Add Cluster.
        while (length(neighbours) != 0) {
          nodeList <- c()
          for (ni in neighbours) {
            clustnum <- which(label == as.character(cluster))
            node <- rneighbors(nb = nb[[ni]],
                               vals = vals,
                               mean_val = lapply(vals, function(m) {
                                 m$D[clustnum]
                               }),
                               label = label,
                               mode = "step", FUN = FUN)
            nodeList <- c(nodeList, node)
            
            label[node] <- ifelse(label[node] == "Noise",
                                  as.character(cluster),
                                  label[node])
          }
          neighbours <- sort(unique(nodeList))
        }
      }
    }
  }
  message("\n", date(), "  Completed.")
  return(label)
}


#' retrive neighbors function.
#'
#' 隣接点を探します。
#'
#' @param nblist `stnb` 関数の出力結果を指定します。
#' @param vals `list(list(D = ..., delta_eps))` のように指定します。
# @param num アクセス番号を指定します。この番号は 1 ~ `nb` の長さの間の整数値を指定します。
#' @param mean_val 平均値または平均値を算出するための数値ベクトルを指定します。
#' @param label クラスタリングのラベルを指定します。NULLの場合は隣接点がそのまま出力され、labelベクトルを指定すると、"Noise" となっている番号をピックアップして出力します。
#' @param mode
#' * "all": 全ての計算を一括で行います。各値と平均値との差をベクトル演算するため、計算回数は少ない反面クラスター平均からの距離 (Δeps) は不正確なります。平均値が変わらない場合に利用可能です。
#' * "step": 逐次計算を実行し、各値と平均値との差の計算を毎回実行し、クラスター平均からの距離 (Δeps) に正確に返します。
#' @param FUN 集計用の計算関数を指定します。
rneighbors <- function(nblist,
                       vals,
                       #num,
                       mean_val = list(),
                       label = NULL,
                       mode = "step",
                       FUN = mean) {
  #nblist <- nb[[num]]
  base <- list()
  for (i in seq(1, length(vals), by = 1)) {
    m <- vals[[i]]
    base[[i]] <- rneighbors_one(nb = nblist,
                                x = m$D[nblist],
                                deps = m$delta_eps,
                                mean_val = mean_val[[i]],
                                mode = mode,
                                FUN = FUN)
  }
  ret <- if (length(vals) > 1) {
    intersects(base)
  } else if (length(vals) == 1) {
    base[[1]]
  }
  if (!is.null(label)) {
    ret[label[ret] == "Noise"]
  } else {
    return(ret)
  }
}


rneighbors_one <- function(nb, x, deps, mean_val,
                           mode = "all", FUN = mean) {
  if (mode == "all") {
    if (length(mean_val) > 1) {
      #mean_val <- mean(mean_val)
      mean_val <- FUN(mean_val)
    }
    return(nb[abs(x - mean_val) <= deps])
  } else if (mode == "step") {
    stopifnot(length(mean_val) > 1)
    ret <- rep(NA, length(x))
    #local_mean_value <- mean(mean_val)
    local_mean_value <- FUN(mean_val)
    for (i in seq(1, length(x), by = 1)) {
      ch <- abs(x[i] - local_mean_value) <= deps
      if (ch) {
        ret[i] <- TRUE
        mean_val <- c(mean_val, x[i])
        #local_mean_value <- mean(mean_val)
        local_mean_value <- FUN(mean_val)
      } else {
        ret[i] <- FALSE
      }
    }
    return(nb[ret])
  }
}


#' create neighbor list
#'
#' 隣接点行列を作成します。
#'
#' @param x double(), sf, sfc
#' @param y double()
#' @param time numeric, date, ISOdate, POSIXct
#' @param eps0 dist x and y (lower)
#' @param eps1 dist x and y (uppwer)
#' @param eps2 dist time. set seconds.
#' @param method use spdep::dnearneigh.
#' * "geo": calc for geographics. 
#' * "euclidean": calc for euclidean area.
#' @param neighbortype
#' * "spatial": For Grid data. Adjust nb size before calc.
#' * "random": For normal(random point) data. not Adjustment.
#' 
stnb <- function(x = double(),
                 y = double(),
                 time = NULL,
                 eps0 = 0,
                 eps1 = NULL,
                 eps2 = NULL,
                 method = c("geo", "euclidean"),
                 neighbortype = c("spatial", "random")
                 ) {
  stopifnot(!is.null(x))
  stopifnot(!is.null(time))
  stopifnot(!is.null(eps1))
  stopifnot(!is.null(eps2))
  method <- match.arg(method)
  neighbortype <- match.arg(neighbortype)
  
  longlat <- ifelse(method == "geo", TRUE, FALSE)
  if (any(class(x) %in% c("sf", "sfc"))) {
    g <- spdep::dnearneigh(x = x,
                           d1 = eps0, d2 = eps1)
  } else {
    g <- spdep::dnearneigh(x = cbind(x, y),
                           d1 = eps0, d2 = eps1,
                           longlat = longlat)
  }
  timediff <- lapply(time, function(tn) {
    which(abs(time - tn) <= eps2)
  })
  
  if (neighbortype == "spatial") {
    if (any(class(x) == "sf")) {
      len <- nrow(x)
    } else {
      len <- length(x)
    }
    ret <- do.call(c, lapply(timediff, function(num) {
      lapply(g, function(gd) {
        unique(do.call(c, Map(`+`, list(gd), (num - 1) * len)))
      })
    }))
  } else if (neighbortype == "random") {
    ret <- Map(intersect, g, timediff)
  }

  return(structure(Map(as.integer, ret),
                   class = "nb"))
}


intersects <- function(x = list()) {
  if (length(x) == 1) {
    return(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else {
    intersects(c(list(intersect(x[[1]], x[[2]])), x[-c(1:2)]))
  }
}


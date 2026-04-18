#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom sf st_as_sf
#' @importFrom sf st_coordinates
#' @importFrom sf st_geometry
#' @importFrom sf st_sf
#' @importFrom spdep dnearneigh
#' @importFrom spdep knearneigh
#' @importFrom spdep knn2nb
#' @importFrom spdep poly2nb
#' @importFrom spdep nbdists
#' @importFrom dbscan dbscan
#' @importFrom dbscan frNN
#' @importFrom units set_units
## usethis namespace: end
NULL


#' set stdbscan's value
#'
#' 時間・空間軸に該当するデータと、このデータのクラスターを判別するための閾値 (Δeps) を指定します。
#'
#' @export
#' @param name データを判別するための名前を文字列で指定します。
#' @param v 時間・空間軸のデータをベクトル形式で指定します。
#' @param deps クラスターを判別するための閾値 (Δeps) を数値で指定します。
#' @return
#' * v: 時間・空間軸のデータ
#' * deps: クラスターを判別するための閾値 (Δeps) 
#' attributes
#' * name: データ名
#' @examples
#' setOneVals(name = "ox", v = abs(runif(100)) * 100, deps = 2)
setOneVals <- function(name = NULL,
                       v = NULL,
                       deps = NULL) {
  stopifnot(any(is.vector(v), is.matrix(v)))
  if (is.null(v)) {
    v <- 0
  }
  if (is.null(deps)) {
    deps <- 1
  }
  if (is.null(name)) {
    name <- "stdbscan raw data"
  }
  return(list(structure(list(v = v,
                             deps = deps),
                        class = "raw_stdbscan",
                        label = name)))
}
#' print raw stdbscan data
#' @export
#' @param x `setOneVals() ` 関数で作成したデータを指定します。
#' @examples
#' setOneVals(name = "ox", v = abs(runif(100)) * 100, deps = 2)
print.raw_stdbscan <- function(x) {
  message("===== stdbscan raw data =====")
  message("name: ", attr(x, "label"), ", Δeps: ", x$deps, ", data: ")
  str(x$v)
}

#' set some stdbscan's value
#'
#' 複数の時間・空間軸に該当するデータと、このデータのクラスターを判別するための閾値 (Δeps) を作成します。
#'
#' @export
#' @param ... `setOneVlas()` で指定する name, v, deps 引数の値を順番に指定します。
#' 複数のデータを指定する場合は、name, v, deps, name, v, deps, ...のように順番に指定します。
#' @examples
#' setVals("test1", rnorm(100), 2, "test2", matrix(rnorm(100), ncol = 10), 4.1)
setVals <- function(...) {
  x <- list(...)
  if (length(x) == 0) {
    return(NULL)
  } else {
    stopifnot(length(x) %% 3 == 0)
    ret <- do.call(c, lapply(seq(1, length(x), by = 3), function(m) {
      setOneVals(name = x[[m]],
                 v = x[[m + 1]],
                 deps = x[[m + 2]])
    }))
    if (length(unique(lapply(ret, function(m) {length(m$v)}))) > 1) {
      stop("Contains different numbers of data. All must be aligned to the same number of data.")
    }
    return(ret)
  }
}

#' Obtains data at an arbitrary time position from the data created by the setVals function.
#'
#' setVals関数で作成したデータから、任意の時間位置のデータを取得します。
#'
#' @export
#' @param x setVals関数で作成されたデータを指定します。
#' @param timelength 各時間のデータが同数の場合、全時間の長さを指定します。
#' @param time 各時間のデータが同数ではない場合、setVals関数のv変数と同数の時間情報ベクトルを指定します
#' @param pos 抽出する時間の位置または時間の値を指定します。
#' timelengthに値を指定した場合は、何番目の時間のデータを抽出するか指定します。このとき、x$vのデータはncol = timelength のmatrixに変換され、pos 列目のデータを返すようになります。(as.matrix(x$v, ncol = timelength)[, pos])
#' time引数にデータを指定した場合は、抽出する時間値そのものを指定します。データに一致する情報だけを取得します。(x$v[time == pos])
#' @keywords internal
#' @examples
#' vals <- setVals("test1", rnorm(100), 2, "test2", matrix(rnorm(100), ncol = 10), 4.1)
#' getVals(vals, timelength = 10, pos = 2)
#' getVals(vals, time = sample(1:10, 100, replace = TRUE), pos = 2)
getVals <- function(x, timelength = NULL, time = NULL, pos) {
  lapply(x, function(m) {
    m$v <- if (all(!is.null(timelength), is.null(time))) {
      matrix(m$v, ncol = timelength)[, pos]
    } else if (all(is.null(timelength), !is.null(time))) {
      m$v[which(time == pos)]
    }
    return(m)
  })  
}

#' Calculate adjacent points of time information
#'
#' 時間情報ベクトルの隣接情報を作成します。
#'
#' @export
#' @param x sfcクラスのデータを指定します。
#' @param t1 隣接点の最低値を指定します。
#' @param t2 隣接点の最大値を指定します。
#' @examples
#' td <- sample(seq(as.Date("2025-01-01"), as.Date("2025-02-01"), by = "1 day"), 100, replace = TRUE)
#' tnearneigh(x = td, t1 = 2, t2 = 4)
tnearneigh <- function(x, t1 = 0, t2) {
  stopifnot(!is.null(x))
  stopifnot(any(class(t1) %in% c("POSIXct", "POSIXt", "Date", "numeric")))
  stopifnot(any(class(t2) %in% c("POSIXct", "POSIXt", "Date", "numeric")))
  timediff <- lapply(seq(1, length(x)), function(i) {
    d <- abs(x[[i]] - x)
    n <- which(t1 <= d & d <= t2)
    n[n != i]
  })
  return(structure(timediff, class = "nb"))
}

#' calc neighbors.
#'
#' @param x sf class data
#' @param fun Specifies the method used to calculate the adjacency points.
#' * dist: use spdep::dnearneigh function.
#' * knn: use spdep::knearneigh and spdep::knn2nb function.
#' * poly: use spdep::poly2nb function.
#' * eublidean: use spdep::dnearneigh function and x converted sf to data.frame.
#' @param opt set function options.
calc_nb <- function(x, fun = c("dist", "knn", "poly", "euclidean"), opt) {
  if (fun == "dist") {
    message("call: spdep::dnearneigh")
    # spdep::dnearneigh(x = x, d1 = d1, d2 = d2, ...)
    return(do.call(spdep::dnearneigh, c(list(x = x), opt)))
  } else if (fun == "knn") {
    message("call: spdep::knearneigh")
    # spdep::knearneigh(x = x, k = k, ...)
    return(spdep::knn2nb(do.call(spdep::knearneigh, c(list(x = x), opt))))
  } else if (fun == "poly") {
    message("call: spdep::poly2nb")
    # spdep::poly2nb(pl = x, ...)
    return(do.call(spdep::poly2nb, c(list(pl = x), opt)))
  } else if (fun == "euclidean") {
    message("call: spdep::dnearneigh (longlat = FALSE)")
    #return(do.call(dbscan::frNN, c(list(x = sf::st_coordinates(x)), opt)))
    return(do.call(spdep::dnearneigh, c(list(x = sf::st_coordinates(x)), opt)))
  }
}

#' create neighbor list
#'
#' 隣接点行列を作成します。
#'
#' @export
#' @param x `sfc` クラスの位置情報を指定します。
#' @param time numeric, date, ISOdate, POSIXct
#' @param eps0 dist x and y (lower) [km]
#' @param eps1 dist x and y (uppwer) [km]
#' @param eps2 dist time. set seconds.
#' @param neighbortype
#' * "spatial": For Grid data. Adjust nb size before calc.
#' * "random": For normal(random point) data. not Adjustment.
#' @keywords internal
#' @examples
#' eps1 <- 144
#' minPts <- 5
#' num <- 1000
#' vals <- setVals("test1", rnorm(num), 0.2,
#'                 "test2", matrix(rnorm(num), ncol = 10), 0.2)
#' ## ramdom pointdata.
#' x <- sf::st_as_sf(data.frame(lon = runif(num, min = 130, max = 140),
#'                              lat = runif(num, min = 30, max = 40)),
#'                   coords = c("lon", "lat"), crs = 4326) |>
#'   sf::st_geometry()
#' (nb <- stnb(x = x, eps1 = eps1, time = NULL, neighbortype = "random"))
#' eps2 <- 3600 * 3
#' t <- seq(as.POSIXct("2020-01-01 03:00:00"), as.POSIXct("2020-01-01 12:00:00"), by = "3 hour")
#' vals <- setVals("test1", rnorm(num * length(t)), 0.2,
#'                 "test2", matrix(rnorm(num * length(t)), ncol = 10), 0.2)
#' x <- sf::st_as_sf(data.frame(lon = runif(num * length(t), min = 130, max = 140),
#'                              lat = runif(num * length(t), min = 30, max = 40)),
#'                   coords = c("lon", "lat"), crs = 4326) |>
#'   sf::st_geometry()
#' (nb <- stnb(x = x, eps1 = eps1, eps2 = eps2, time = sample(t, num, replace = TRUE), neighbortype = "random"))
#' 
#' ## spatial point
#' lon <- seq(130, 140, by = 1)
#' lat <- seq(30, 40, by = 1)
#' x <- sf::st_as_sf(expand.grid(lon, lat), coords = c("Var1", "Var2"), crs = 4326) |>
#'   sf::st_geometry()
#' (nb <- stnb(x = x, eps1 = eps1, eps2 = eps2, time = sample(t, length(x), replace = TRUE), neighbortype = "spatial"))
stnb <- function(x = NULL,
                 time = NULL,
                 eps0 = 0,
                 eps1 = NULL,
                 eps2 = NULL,
                 k = NULL,
                 neighbortype = c("random", "spatial"),
                 fun = "dist",
                 ...) {
  stopifnot(any(class(x) %in% "sfc"))
  stopifnot(!is.null(eps1))
  neighbortype <- match.arg(neighbortype)
  fun <- match.arg(fun, c("dist", "knn", "poly", "euclidean"))
  
  if (fun == "dist") {
    opt <- list(d1 = eps0,
                d2 = eps1,
                ...)
  } else if (fun == "knn") {
    opt <- list(k = k, ...)
  } else if (fun == "poly") {
    opt <- list(...)
  } else if (fun == "euclidean") {
    opt <- list(eps = eps1, ...)
  }
  
  timevar <- sort(unique(time))
  # if (any(is.null(time), length(timevar) == 1)) {
  #   sp_nb <- calc_nb(x = x, fun = fun, opt = opt)
  #   ts_nb <- NA # tnearneigh(x = x, t2 = eps2)
  # } else if (all(!is.null(time), length(timevar) > 1)) {
  #   stopifnot(!is.null(eps2))
  #   ts_nb <- tnearneigh(timevar, t2 = eps2)
  #   if (all(neighbortype == "random", parallel)) {
  #     spdata <- split(x, time)[as.character(timevar)]
  #     sp_nb <- lapply(spdata, function(m) {
  #       calc_nb(x = m, fun = fun, opt = opt)
  #     })
  #     ts_nb <- tnearneigh(timevar, t2 = eps2)
  #   } else if (all(neighbortype == "spatial", parallel)) {
  #     sp_nb <- calc_nb(x = x, fun = fun, opt = opt)
  #     ts_nb <- tnearneigh(timevar, t2 = eps2)
  #   } else if (all(neighbortype == "random", !parallel)) {
  #     sp_nb <- structure(Map(intersect,
  #                            calc_nb(x = x, fun = fun, opt = opt),
  #                            tnearneigh(time, t2 = eps2)),
  #                        class = "nb")
  #     ts_nb <- NA
  #   } else if (all(neighbortype == "spatial", !parallel)) {
  #     d_nb <- calc_nb(x = x, fun = fun, opt = opt)
  #     sp_nb <- do.call(c, lapply(tnearneigh(timevar, t2 = eps2), function(m) {
  #       ln <- length(d_nb) * (m - 1)
  #       lapply(d_nb, function(aa) {
  #         as.integer(unique(do.call(c, Map(`+`, list(aa), ln))))
  #       })
  #     }))
  #     sp_nb <- structure(sp_nb, class = "nb")
  #     ts_nb <- NA
  #   }
  # }
  if (any(is.null(time), length(timevar) == 1)) {
    sp_nb <- calc_nb(x = x, fun = fun, opt = opt)
  } else if (all(!is.null(time), length(timevar) > 1)) {
    stopifnot(!is.null(eps2))
    # ts_nb <- tnearneigh(timevar, t2 = eps2)
    if (neighbortype == "random") {
      sp_nb <- structure(Map(intersect,
                             calc_nb(x = x, fun = fun, opt = opt),
                             tnearneigh(time, t2 = eps2)),
                         class = "nb")      
    } else if (neighbortype == "spatial") {
      d_nb <- calc_nb(x = x, fun = fun, opt = opt)
      sp_nb <- do.call(c, lapply(tnearneigh(timevar, t2 = eps2), function(m) {
        ln <- length(d_nb) * (m - 1)
        lapply(d_nb, function(v) {
          as.integer(unique(do.call(c, Map(`+`, list(v), ln))))
        })
      }))
      sp_nb <- structure(sp_nb, class = "nb")
    }
  }
  if (all(!is.null(time), length(timevar) > 1, neighbortype == "spatial")) {
    ln <- length(x)
    msp <- lapply(sp_nb, function(m) {
      ret <- m %% ln
      ifelse(ret == 0, ln, ret)
    })
    sp_dist <- do.call(c, lapply(split(msp, rep(timevar, each = length(x))), function(m) {
      class(m) <- "nb"
      spdep::nbdists(nb = m, coords = x)
    }))
    names(sp_dist) <- NULL
    # sp_dist <- spdep::nbdists(nb = sp_nb, coords = rep(x, length(timevar)))
  } else {
    sp_dist <- spdep::nbdists(nb = sp_nb, coords = x)
  }
  return(structure(list(dist = sp_dist,
                        id = sp_nb,
                        eps = eps1,
                        metric = fun,
                        neighbortype = neighbortype,
                        coord = x,
                        time = time,
                        timevar = timevar),
                   class = c("frNN", "NN") #class = "stnb"
                   ))
}
# stnb <- function(x = NULL,
#                  time = NULL,
#                  eps0 = 0,
#                  eps1 = NULL,
#                  eps2 = NULL,
#                  neighbortype = c("random", "spatial"),
#                  parallel = FALSE) {
#   stopifnot(any(class(x) %in% "sfc"))
#   stopifnot(!is.null(eps1))
#   neighbortype <- match.arg(neighbortype)

#   ## convert m -> km
#   eps0 <- eps0 * 10 ^ 3 
#   eps1 <- eps1 * 10 ^ 3

#   timevar <- sort(unique(time))
#   if (any(is.null(time), length(timevar) == 1)) {
#     sp_nb <- spdep::dnearneigh(x = x,
#                                d1 = eps0,
#                                d2 = eps1)
#     ts_nb <- NA # tnearneigh(x = x, t2 = eps2)
#   } else if (all(!is.null(time), length(timevar) > 1)) {
#     stopifnot(!is.null(eps2))
#     ts_nb <- tnearneigh(timevar, t2 = eps2)
#     if (all(neighbortype == "random", parallel)) {
#       spdata <- split(x, time)[as.character(timevar)]
#       sp_nb <- lapply(spdata, function(m) {
#         spdep::dnearneigh(x = m,
#                           d1 = eps0,
#                           d2 = eps1)
#       })
#       ts_nb <- tnearneigh(timevar, t2 = eps2)
#     } else if (all(neighbortype == "spatial", parallel)) {
#       sp_nb <- spdep::dnearneigh(x = x,
#                                  d1 = eps0,
#                                  d2 = eps1)
#       ts_nb <- tnearneigh(timevar, t2 = eps2)
#     } else if (all(neighbortype == "random", !parallel)) {
#       sp_nb <- structure(Map(intersect,
#                              spdep::dnearneigh(x = x,
#                                                d1 = eps0,
#                                                d2 = eps1),
#                              tnearneigh(time, t2 = eps2)),
#                          class = "nb")
#       ts_nb <- NA
#     } else if (all(neighbortype == "spatial", !parallel)) {
#       d_nb <- spdep::dnearneigh(x = x,
#                                 d1 = eps0,
#                                 d2 = eps1)
#       sp_nb <- structure(do.call(c, lapply(tnearneigh(timevar, t2 = eps2), function(m) {
#         unique(do.call(c, Map(`+`, list(d_nb), (length(d_nb) * (m - 1)))))
#       })), class = "nb")
#       ts_nb <- NA
#     }
#   }
#   return(structure(list(sp_neigh = sp_nb,
#                         ts_neigh = ts_nb,
#                         neighbortype = neighbortype,
#                         coord = x,
#                         time = time,
#                         timevar = timevar),
#                    class = "stnb"))
# }



#' retrive neighbors function.
#'
#' 隣接点を探します。
#'
#' @param vals `list(list(D = ..., delta_eps))` のように指定します。
#' @param nb 隣接点を vector 形式で指定します。
#' @param cluster_pos valsにおける同一クラスターへのインデックスを指定します。
#' @param FUN 集計用の計算関数を指定します。
#' @keywords internal
#' @note
#' eps1 <- 144
#' eps2 <- 3600 * 6
#' minPts <- 5
#' vals <- setVals("test1", rnorm(100), 0.2,
#'                 "test2", matrix(rnorm(100), ncol = 10), 0.2)
#' ## ramdom point and non temporal data.
#' x <- sf::st_as_sf(data.frame(lon = runif(100, min = 130, max = 140),
#'                              lat = runif(100, min = 30, max = 40)),
#'                   coords = c("lon", "lat"), crs = 4326) |>
#'   sf::st_geometry()
#' nb <- stnb(x = x, eps1 = eps1/1000, eps2 = eps2, time = NULL, neighbortype = "random")
#' message("neighbors: ")
#' nb$id[[1]][[1]]
#' message("cluster neighbors: ")
#' rneighbors(vals = vals, nb = nb$id[[1]][[1]], cluster_pos = 1, FUN = mean)
rneighbors <- function(vals, nb, cluster_pos = NULL, FUN) {
  rep_num <- Reduce(intersect, lapply(vals, function(m) {
    d <- FUN(m$v[cluster_pos])
    which((m$v[nb] - d) <= m$deps)
  }))
  return(nb[rep_num])
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
#' @param dbscantype
#' * "grid": For Grid data. Adjust nb size before calc.
#' * "default": For normal(random point) data. not Adjustment.
#' @param FUN 集計用の計算関数を指定します。
#' @param na_value 元データにNA値が含まれている場合に代替する値を指定します。
#' なお、データは代替して計算されますが、NA値のクラスタリング結果は Noise 判定となります。
#' Noise判定にしたくない場合はこの機能を利用せず、引数指定時にNA値を除去または代替したデータを設定するようにしてください。
#' @keywords internal
#' @note
#' eps1 <- 144 ## km
#' eps2 <- 3600 * 6
#' minPts <- 20
#' FUN <- mean
#' num <- 10000
#' vals <- setVals("test1", rnorm(num, mean = 10, sd = 1), 0.05,
#'                 "test2", matrix(rnorm(num, mean = 100, sd = 10), ncol = 10), 5)
#' ## ramdom point and non temporal data.
#' x <- sf::st_as_sf(data.frame(lon = runif(num, min = 130, max = 140),
#'                              lat = runif(num, min = 30, max = 40)),
#'                   coords = c("lon", "lat"), crs = 4326) |>
#'   sf::st_geometry()
#' nb <- stnb(x = x, eps1 = eps1, eps2 = eps2, time = NULL, neighbortype = "random")
#' cl <- calc_stdbscan(nb = nb$id, vals = vals, dbscantype = "default", minPts = minPts)
calc_stdbscan <- function(nb = NULL,
                          vals = setVals(),
                          dbscantype = c("default", "grid"),
                          minPts = 5,
                          FUN = mean,
                          na_value = -9999) {
  dbscantype <- match.arg(dbscantype)
  stopifnot(is.list(vals))
  stopifnot(!is.null(nb))
  
  ## set cluster number
  cluster <- 0

  ## set cluster column
  # if contain NA value the create "cluter_NA" value (= -9999) at that position.
  # otherwise, create "Noise" label.(= 0)
  label <- ifelse(apply(do.call(cbind,
                                Map(function(x) as.vector(x$v),
                                    vals)), 1,
                        function(x) any(is.na(x))),
                  -9999, 0)
  if (all(do.call(c, lapply(vals, function(m) is.matrix(m$v))))) {
    label <- matrix(label, nrow = nrow(vals[[1]]$v), ncol = ncol(vals[[1]]$v))
  }
  
  ## convert na_value
  vals <- lapply(vals, function(m) {
    m$v <- ifelse(is.na(m$v), na_value, m$v)
    return(m)
  })
  
  ## st-dbscan
  pb <- txtProgressBar(min = 0, max = length(nb), style = 3)
  for (i in seq(1, length(nb))) {
    setTxtProgressBar(pb, i)
    if (any(label[i] == na_value, 0 < label[i])) {
      next
    }
    
    neigh <- if (dbscantype == "grid") {
      rneighbors(vals = vals, nb = nb[[i]], cluster_pos = nb[[i]], FUN = FUN)
    } else if (dbscantype == "default") {
      nb[[i]]
    }
    neigh <- neigh[which(label[neigh] == 0)]
    if (length(neigh) <= minPts) {
      next
    }

    ## up to cluster number
    cluster <- cluster + 1

    ## set cluster
    label[c(i, neigh)] <- cluster
    
    ## check node and Add Cluster.
    while (length(neigh) != 0) {
      clustpos <- which(label == cluster)
      merge_neigh <- Reduce(function(x, y) setdiff(x, y), list(unique(unlist(nb[neigh])),
                                                               neigh,
                                                               which(label > 0)))
      node <- rneighbors(vals = vals,
                         nb = merge_neigh,
                         cluster_pos = clustpos,
                         FUN = FUN)
      label[node] <- ifelse(label[node] == 0 | label[node] == -9999,
                            cluster,
                            label[node])
      neigh <- sort(unique(node))
      }
  }
  close(pb)
  
  # ## Rcpp::sourceCpp("./src/stdbscan.cpp")
  # label <- stdbscan_core(nb = nb, vals = vals, dbscantype = dbscantype,
  #                        minPts = minPts, na_value = na_value)

  clust_label_level <- unique(label)
  clust_label_level <- clust_label_level[which(clust_label_level != na_value & clust_label_level > 0)]
  att <- lapply(clust_label_level, function(m) {
    num <- which(label == m)
    valnames <- do.call(c, lapply(vals, function(vls) attr(vls, "label")))
    list(name = m,
         pos = num,
         values = Map(assign,
                      valnames,
                      lapply(vals, function(vls) vls$v[num]))
         )
  })


  ## Convert "Noise" label
  label <- ifelse(label == 0 | label == na_value, "Noise", as.character(label))
  return(structure(label,
                   cluster = unique(label),
                   deps = Map(assign,
                              do.call(c, lapply(vals, function(vls) attr(vls, "label"))),
                              lapply(vals, function(vls) vls$deps)),
                   FUN = substitute(FUN),
                   info = att))
}


#' calc neighbors and ST-DBSCAN.
#'
#' @export
#' @param x 位置情報データをsf, sfc, matrix, data.frameのい何れかで指定します。
#' matrix, data.frame で指定する場合、経度をx, 緯度をyとして列名 (colnames, names) 、データを作成してください。
#' また、matrix, data.frameで指定した場合は、内部でsfcクラスのデータに変換されます。
#' @param time 時間情報を指定します。
#' `neighbortype = "random"` の場合、x (位置情報) の数とtimeの長さが一致している必要があります。 (`length(x) == length(time)`)
#' `neighbortype = "spatial"` の場合、`length(x) * length(time) ` がデータのサイズになります。
#' @param eps1 空間近傍地点のしきい値を距離 (km) で指定します。
#' @param eps2 時間近傍地点のしきい値を数値で指定します。
#' time変数に数値を指定指定した場合は、eps2も数値で指定します。
#' time変数をDate型で指定した場合は、日差 (day) として指定します。
#' time変数をPOSIXct型で指定した場合は、秒差 (sec) として指定します。
#' @param minPts クラスター判定するための最小近傍地点数を整数で指定します。
#' @param vals クラスタリングする属性値を指定します。
#' ここでは複数属性のデータを指定することができ、`setVals("data-name", valus(vector), Δeps)` として指定することが出来ます。
#' 指定の方法の詳細は `setOneVals()` 関数を参照してください。
#' @param neighbortype 空間隣接地点の計算方法を指定します。
#' "random" と指定した場合、位置情報がランダムであるデータの計算を実施します。
#' このとき、位置情報 (x) の地点数と 属性値 (vals) のデータ数が同数である必要があります。
#' timeを指定する場合は、timeの長さも位置情報の地点数と同数である必要があります。
#' (length(x) == length(time) == length(vals[[1]]$v))
#' "spatial" と指定した場合、位置情報 (x) の地点数と時間情報 (time) の数を乗じた値が属性値のデータ数と等しくなる必要がります。
#' (length(x) * length(time) == length(vals[[1]]$v))。
#' @param dbscantype 空間隣接地点数の補足設定を行います。
#' WRFに代表される気象シミュレーションのように、領域を格子状に分割して計算された (netCDF) データでは、地点の位置情報は同じで時間情報だけが異なるデータ構造になっています。
#' このようなデータでは一定以上の隣接地点数を見込むことができるため、必ずクラスターが作成されるという弊害が存在します。
#' そのため、 `dbscantype = "grid"` と指定した場合に限り、クラスターの作成条件を本来のST-DBSCANより厳しくします。
#' 具体的に、ST-DBSCANのアルゴリズムでは 隣接地点数が minPts 以上である場合に無条件でクラスター化しますが、このオプションを指定することで属性値の計算結果を判定条件とするようにします。
#' 特に指定がない場合、`neighbortype == spatial` の場合、自動で `dbscantype = "grid"` として計算されます。
#' @param FUN 同一のクラスターを判別するための方法を定義します。
#' 文献では、クラスターの属性値平均値との差 (Δeps) を確認し、新たな地点を追加しています。標準の設定ではこの内容に従い、 `kean` 関数をが指定されています。
#' クラスターへの追加方法を変更する場合 (例えば中央値との差を見るようにする場合) は、 `FUN = mesian` のように指定します。
#' @param na_value 属性情報 vals にNA値が含まれている場合、そのままではクラスタリングはエラーになります。
#' そのため、一時的に属性情報を na_valueで上書きしてクラスタリングの計算を行います。
#' この値は属性情報とは関係がない値にすることが望ましく、例えば、属性情報がプラスの値であれば na_valueは大きなマイナスの値にすると Nose として判定され舞うs。
#' 標準では `na_value = -9999` に指定されています。
#' @param nb_fun 隣接地点の計算方法を指定します。
#' * dist: 一定の距離 (km) に基づく隣接点を指定します。
#' * knn: 近い方から k 地点を隣接点とします。別途、 `k = ` 引数で条件を指定して下さい。
#' * poly: ポリゴン同士が接している場合に隣接点とみなします。この場合、 x 引数には sfc_MULTIPOLYGON 属性のデータを指定する必要があります。
#' * euclidean: ユークリッド距離による座標間の直接的な距離を元に隣接点を計算します。
#' @param k `nb_fun = "knn"` の場合に指定します。整数値で近傍点数を指定して下さい。
#' @param ... 隣接点の計算オプションを指定します。以下の関数で使用できる引数オプションを指定できます。
#' spdep::dnearneigh, spdep::knearneigh, spdep::poly2nb, dbscan::frNN
#' 
#' @note
#' eps1 <- 144 ## km
#' eps2 <- 3600 * 6
#' time <- NULL
#' minPts <- 20
#' FUN <- mean
#' num <- 1000
#' 
#' ## only geospatial data (dbscan).
#' x <- sf::st_as_sf(data.frame(lon = runif(num, min = 130, max = 140),
#'                              lat = runif(num, min = 30, max = 40)),
#'                   coords = c("lon", "lat"), crs = 4326) |>
#'   sf::st_geometry()
#' clust <- stdbscan(x = x, eps1 = eps1, minPts = minPts)
#' 
#' ## ramdom point and non temporal data.
#' vals <- setVals("test1", rnorm(num, mean = 10, sd = 1), 0.05,
#'                 "test2", matrix(rnorm(num, mean = 100, sd = 10), ncol = 10), 5)
#' x <- sf::st_as_sf(data.frame(lon = runif(num, min = 130, max = 140),
#'                              lat = runif(num, min = 30, max = 40)),
#'                   coords = c("lon", "lat"), crs = 4326) |>
#'   sf::st_geometry()
#' clust <- stdbscan(x = x, eps1 = eps1, minPts = minPts,
#'                   vals = vals, neighbortype = "random")
#' 
#' # random point and temporal data.
#' time <- as.POSIXct("2024-01-01 00:00:00", tz = "JST") + 3600 * seq(0, 23, by = 6)
#' vals <- setVals("test1", rnorm(num * length(time), mean = 10, sd = 1), 0.5,
#'                 "test2", matrix(rnorm(num * length(time), mean = 100, sd = 10), ncol = 10), 5)
#' x <- sf::st_as_sf(data.frame(lon = runif(num * length(time), min = 130, max = 140),
#'                              lat = runif(num * length(time), min = 30, max = 40)),
#'                   coords = c("lon", "lat"), crs = 4326) |>
#'   sf::st_geometry()
#' clust <- stdbscan(x = x, time = rep(time, each = num), eps1 = eps1, eps2 = eps2, minPts = minPts,
#'                   vals = vals, neighbortype = "random")
#' 
#' ## grid-base spatial point and temporal data.
#' time <- as.POSIXct("2024-01-01 00:00:00", tz = "JST") + 3600 * seq(0, 23, by = 6)
#' vals <- setVals("test1", rnorm(length(x) * length(time), mean = 10, sd = 1), 0.5,
#'                 "test2", matrix(rnorm(length(x) * length(time), mean = 100, sd = 10), ncol = length(time)), 5)
#' x <- sf::st_as_sf(expand.grid(seq(130, 140, by = 1), seq(30, 40, by = 1)),
#'                   coords = c("Var1", "Var2"), crs = 4326) |>
#'   sf::st_geometry()
#' clust <- stdbscan(x = x, time = time, eps1 = eps1, eps2 = eps2, minPts = 4,
#'                   vals = vals, neighbortype = "spatial")
stdbscan <- function(x,
                     time = NULL, 
                     eps1 = double(),
                     eps2 = double(),
                     minPts = integer(),
                     vals = setVals(),
                     neighbortype = "random",
                     dbscantype = NULL,
                     FUN = mean,
                     na_value = -9999,
                     nb_fun = "dist",
                     k = integer(),
                     ...) {
  ## check data
  stopifnot(!is.null(x))
  stopifnot(length(minPts) > 0)
  neighbortype <- match.arg(neighbortype, c("random", "spatial"))
  if (all(is.null(dbscantype), neighbortype == "spatial", !is.null(time))) {
    dbscantype <- "grid"
  } else {
    dbscantype <- "default"
  }
  dbscantype <- match.arg(dbscantype, c("default", "grid"))

  ## check neighbortype
  ## convert x to sf class.
  if (any("matrix" %in% class(x), class(x) == "data.frame")) {
    x <- sf::st_as_sf(x, coords = c("x", "y"), crs = 4326)
  } 
  x <- sf::st_geometry(x)

  if (is.null(vals)) {
    scantype <- "dbscan"
  } else {
    scantype <- "stdbscan"
    if (is.null(time)) {
      stopifnot(length(x) == length(vals[[1]]$v))
    } else if (all(neighbortype == "random", !is.null(time))) {
      stopifnot(length(x) == length(time))
    } else if (all(neighbortype == "spatial", !is.null(time))) {
      stopifnot(length(x) * length(time) == length(vals[[1]]$v))
    }
  }

  if (nb_fun == "dist") {
    stopifnot(length(eps1) > 0)
    if (!is.null(time)) {
      stopifnot(length(eps2) > 0)
    }    
  } else if (nb_fun == "knn") {
    stopifnot(k > 0)
  } else if (nb_fun == "poly") {
    stopifnot("sfc_MULTIPOLYGON" %in% class(x))
  } else if (nb_fun == "euclidean") {
    stopifnot(length(eps1) > 0)
  }

  if (all(nb_fun %in% c("dist", "knn"), !"sfc_POINT" %in% class(x))) {
    warning("`stdbscan` requires POINT class data.")
    warning("Apply the `sf::st_controid` function and use the position of the polygon's center of gravity.")
    archive_x <- x
    x <- sf::st_centroid(x)
  }

  
  message("===== Start ST-DBSCAN method =====")
  message("\n1. Calculation Neighbor List")

  ## calc neighbors
  if (scantype == "stdbscan") {
    nb <- stnb(x = x, eps0 = 0, eps1 = eps1, eps2 = eps2, k = k, fun = nb_fun, 
               time = time, neighbortype = neighbortype, ...)
    # stnb(x = x, eps0 = 0, eps1 = eps1, eps2 = eps2, k = k, fun = "dist", time = time, neighbortype = "random")
  } else if (scantype == "dbscan") {
    if (all(nb_fun == "euclidean", !is.null(time))) {
      nb <- dbscan::frNN(x = sf::st_coordinates(x), eps = eps1, ...)
      nb$timevar <- NULL
    } else {
      nb <- stnb(x = x, eps0 = 0, eps1 = eps1, eps2 = eps2, k = k, fun = nb_fun,
                 time = time, neighbortype = neighbortype, ...)
    }    
  }
  
  ## calc ST-DBSCAN
  message("\n2. Calculation Cluster")
  if (scantype == "dbscan") {
    message("* Since `vals` is not specified, the calculation is performed by dbscan.")
    #cluster <- dbscan::dbscan(x = sf::st_coordinates(x), eps = eps1, minPts = minPts, ...)
    dbscan_out <- dbscan::dbscan(x = nb, minPts = minPts, ...)
    cluster <- as.character(dbscan_out$cluster)
  } else if (scantype == "stdbscan") {
    message("* Performs ST-DBSCAN calculations.")
    cluster <- calc_stdbscan(nb = nb$id, vals = vals, na_value = na_value,
                             dbscantype = dbscantype, minPts = minPts, FUN = FUN)
  }
  raw_data <- data.frame(id = seq(1, length(x)),
                         time = local({
                           if (all(neighbortype == "random", !is.null(time))) {
                             time
                           } else if (all(neighbortype == "random", is.null(time))) {
                             NA
                           } else if (neighbortype == "spatial") {
                             rep(nb$timevar, each = length(x))
                           }
                         }),
                         cluster = cluster)
  if (!is.null(vals)) {
    raw_data <- cbind(raw_data,
                      do.call(data.frame,
                              Map(function(m, n) assign(m, n),
                                  unlist(Map(function(m) attr(m, "label"), vals)),
                                  Map(function(m) as.vector(m$v), vals))))
  }
  if (neighbortype == "spatial") {
    raw_data <- split(raw_data, raw_data$time)
    names(raw_data) <- nb$timevar
  }  
  if (exists("archive_x")) {
    x <- archive_x
  }
  
  ret <- structure(list(cluster = as.character(cluster),
                        eps1 = units::set_units(eps1, "km"),
                        eps2 = local({
                          if (all(class(time) == "Date")) {
                            units::set_units(eps2, "d")
                          } else if (any(class(time) == "POSIXct")) {
                            units::set_units(eps2, "s")
                          } else {
                            eps2
                          }
                        }),
                        minPts = minPts,
                        FUN = substitute(FUN), 
                        metric = list(metric = local({
                                        if (nb_fun == "dist") {
                                          "geo(spdep::dnearneigh)"
                                        } else if (nb_fun == "knn") {
                                          "geo(spdep::knearneigh, spdep::knn2nb)"
                                        } else if (nb_fun == "poly") {
                                          "geo(spdep::poly2nb)"
                                        } else if (any(nb_fun == "euclidean", scantype == "dbscan")) {
                                          "euclidean(spdep::dnearneigh)"
                                        }
                                      }),
                                      neighbortype = neighbortype,
                                      dbscantype = dbscantype),
                        neighborlist = nb,
                        deps = attr(cluster, "deps"),
                        vals = vals,
                        cluster_info = attr(cluster, "info"),
                        results = list(geo = sf::st_sf(id = seq(1, length(x)),
                                                       geometry = x),
                                       time = nb$timevar,
                                       value = raw_data)),
                   class = "stdbscan")
  message("\nCompleted.")
  return(ret)
}

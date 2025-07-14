#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom sf st_geometry
#' @importFrom sf st_as_sf
#' @importFrom sparseDFM sparseDFM
#' @importFrom sparseDFM tuneFactors
#' @importFrom spdep knn2nb
#' @importFrom spdep joincount.multi
#' @importFrom spdep joincount.mc
#' @importFrom dplyr left_join
## usethis namespace: end
NULL


#' check list for Spatio-Temporal data.
#'
#' 時系列データの類似性及び空間データの集塊性の確認の基準を設定します。
#'
#' @export
#' @param ts `sparseDFM::sparseDFM` 関数で計算された結果について、R ^ 2値 (rsq) 及びRMSE (rmse) 値の基準を設定します。
#' `upper`, `middle`, `lower` の3段階の設定ができます。
#' @param sp `spdep::joincount.mc` 及び `spdep::joincount.multi` 関数で計算された結果について、z値 (z_value) 及び調整済みp値 (adj_p_value) の基準を設定します。
#' `upper`, `middle`, `lower` の3段階の設定ができます。
checkCond <- function(ts = list(upper = "rsq >= 0.8 & 6 >= rmse",
                                middle = "rsq >= 0.6 & 9 >= rmse",
                                lower = "rsq >= 0.4 & 12 >= rmse"),
                      sp = list(upper = "z_value > 0 & 0.001 >= adj_p_value",
                                middle = "z_value > 0 & 0.01 >= adj_p_value",
                                lower = "z_value > 0 & 0.05 >= adj_p_value")) {
  list(ts = ts, sp = sp)
}


#' evaluate the time series similarity and spatial aggregation of clustering results.
#'
#' hstdbscanなどの結果について、時系列性と空間的集塊性の評価を行います。
#'
#' @export
#' @param x クラスタリングに指定したデータセット (配列)。列: 位置情報、行: 時系列データとして作成
# @param checkmode クラスタリングの精度検証 について、時系列 ("temporal"), 空間 ("spatial"), またはこれらの両方を指定します。
#' @param cc `checkCond` 関数をで、時系列及び空間データの精度確認基準を指定します。
#' @param cluster クラスタリングの結果ラベルを指定します。内部処理ではfactor型に変換されます。
#' @param geo 位置情報
#' @param neighbor_method 隣接点を計算する関数を指定します。
#' @param neighbor_option `neighbor_method` で指定した関数の設定 (引数) を名前付きで指定します。
#' @param ts_num_factor 時系列の動的因子モデル (sparseDFM) の因子数を正の整数値で指定します。
#' 自動で最適な因子数を取得する場合は `NULL` を指定します。
#' @param set_crs `neighbor_method` で指定した関数を実行するために最適なCRS値を指定します。
#' @param nsim `spdep::joincount.mc` の繰り返し回数
#' @return
#' * results: 時系列類似性と空間的集塊性の検定結果の一覧表 (data.frame)
#' * check_mode: 時系列類似性 (temporal) と空間的集塊性 (spatial) の設定
#' * sp_model: Join-Count検定の結果一覧
#' ** joincount: Join-Count検定の結果
#' ** joincount_mc: Joun-Countのモンテカルロ・シミュレーション結果
#' ** neighbor_list: 隣接点情報
#' ** weight_list: 空間重み付き行列 (listw) 
#' * ts_model: 時系列データの動的因子モデルによる結果
#' * geo: 時系列類似性と空間的集塊性の結果を `geo` 引数で指定した位置情報にマージし、sf形式のデータとした
#' * check_conditions: 時系列類似性と空間的集塊性のチェック基準
#' * data: `x` 引数に指定したデータ
#' * cluster: クラスターラベル
#' * neighbor_method: 空間隣接点情報の計算方法
#' * nsim: モンテカルロ・シミュレーションの実行回数
#' @examples
#' x <- seq(130, 140, by = 1)
#' y <- seq(30, 40, by = 1)
#' t <- as.POSIXct("2024-01-01 00:00:00", tz = "JST") + 3600 * seq(0, 23, by = 6)
#' geo <- sf::st_as_sf(expand.grid(x, y), coords = c("Var1", "Var2"), crs = 4326)
#' 
#' D <- abs(runif(nrow(geo) * length(t))) * 100
#' clust <- hstdbscan(x = geo, time = t, eps1 = 144, eps2 = 3600 * 6, minPts = 6,
#'                   vals = list(list(D = D,
#'                                    delta_eps = 20)),
#'                   metric = "geo", neighbortype = "spatial", dbscantype = "grid")
#' stclust(x = matrix(D, ncol = length(t)), cluster = cutclust(clust, 5)$cluster, geo = clust$results$geo)
stclust <- function(x = NULL,
                    check_mode = c("spatial", "temporal"),
                    cc = checkCond(), 
                    cluster = NULL,
                    geo = NULL,
                    neighbor_method = spdep::knearneigh,
                    neighbor_option = list(k = 4),
                    ts_num_factor = 1L,
                    set_crs = 6670,
                    nsim = 999) {
  
  ## checkmode
  ## check_mode <- c("spatial", "temporal")
  check_mode <- match.arg(check_mode,
                          choices = c("spatial", "temporal"),
                          several.ok = TRUE)

  
  ## clustering labels
  cluster <- as.factor(cluster)

  ## geo data
  geo <- if (any(class(geo) %in% c("sf", "sfc"))) {
    sf::st_transform(sf::st_geometry(geo), crs = set_crs)
  } else if (any(class(geo) == "data.frame")) {
    sf::st_geometry(sf::st_as_sf(geo, coords = c("x", "y"),
                                 crs = set_crs))
  } else if (any(class(geo) == "matrix")) {
    local({
      g <- as.data.frame(geo)
      names(g) <- c("x", "y")
      sf::st_geometry(sf::st_as_sf(g, coords = c("x", "y"),
                                   crs = set_crs))
    })
  }
  
  if (any(check_mode == "spatial")) {
    message("1. Check cluster agglomeration properties in geospatial areas.")

    ## Create neighbor list
    neighbor_method <- substitute(neighbor_method)
    neighbor_option <- substitute(neighbor_option)

    ##nb <- do.call(neighbor_method, c(list(x = geo), neighbor_option))
    nb <- eval(bquote(do.call(.(neighbor_method),
                              c(list(x = geo), .(neighbor_option)))))
    if (class(nb) == "knn") {
      nb <- spdep::knn2nb(nb)
    }
    w <- spdep::nb2listw(nb, zero.policy = TRUE)

    ## Spatial Self-Correlation Analysis
    jc_data <- spdep::joincount.multi(cluster, listw = w)
    jc_mc_data <- spdep::joincount.mc(cluster, listw = w, nsim = nsim)

    sp_model <- cbind(cluster = levels(cluster),
                      data.frame(jc_data[paste(levels(cluster), levels(cluster), sep = ":"),
                                         c("Joincount", "z-value")]),
                      local({
                        p_value <- do.call(c, lapply(jc_mc_data, function(jd) {
                          jd$p.value
                        }))
                        adj_p_value <- p.adjust(p_value, method = "fdr")
                        data.frame(p_value = p_value,
                                   adj_p_value = adj_p_value)
                      }))
    names(sp_model)[3] <- "z_value"
    message("Complete 1.")
  }

  if (any(check_mode == "temporal")) {
    message("2. Evaluate the similarity of time series data.")
    ## time-series
    ts_model <- lapply(levels(cluster), function(m) {
      message("\tCheck cluster: ", m)
      row_index <- which(cluster == m)
      ts_data <- if (length(row_index) == 0) {
        NULL
      } else if (length(row_index) == 1) {
        matrix(x[row_index, ], ncol = 1)
      } else {
        ts(t(x[row_index, ]))
      }
      
      ## check r num
      if (is.null(ts_num_factor)) {
        # ts_num_factor <- as.integer(gsub("^.*is *", "",
        #                                  sparseDFM::tuneFactors(ts_data, plot = FALSE)))
        ts_num_factor <- get_num_factor(ts_data)
        message("num: ", ts_num_factor)
      }
      
      ## sparse model
      sparse_model <- if (any(class(ts_data) == "ts")) {
        tryCatch({
          message("fit sparse model: ", m)
          sparseDFM::sparseDFM(ts_data,
                               r = ts_num_factor,
                               alg = "EM")
        }, error = function(e) {
          if (ts_num_factor == 1L) {
            message("error, try 'r = best-num-factor'")
            sparseDFM::sparseDFM(ts_data,
                                 r = get_num_factor(ts_data),
                                 alg = "EM")
          } else if (ts_num_factor > 1L) {
            message("error, try 'r = 1'")
            sparseDFM::sparseDFM(ts_data,
                                 r = 1,
                                 alg = "EM")
          }
        }, finally = {
          message("complete. ")
        })
      } else if (is.null(ts_data)) {
        NULL
      } else {
        txt <- "The sparseDFM model cannot be constructed because there is only one training data."
        message(txt)
        txt
      }
      
      ## fit model
      fit <- if (class(sparse_model) == "sparseDFM") {
        fitted(sparse_model, alpha_index = "best", standardize = FALSE)
      } else if (is.null(sparse_model)) {
        NULL
      } else {
        ts_data
      }

      ## r.squared
      if (!is.null(fit)) {
        r2_d <- rsq(as.vector(ts_data),
                    as.vector(fit))
        rmse_d <- rmse(as.vector(ts_data),
                       as.vector(fit))
      } else {
        r2_d <- rmse_d <- NULL
      }
      ## return
      list(cluster = m,
           ts_model = sparse_model,
           obs_data = ts_data,
           fit_data = fit,
           rsq = r2_d,
           rmse = rmse_d)
    })
    message("Complete 2.")
  }

  ## merge sp_model and ts_model
  m_data <- if (all(check_mode == c("spatial", "temporal"))) {
    merge(as.data.frame(do.call(rbind, lapply(ts_model, function(m) {
      with(m, {
        data.frame(cluster = cluster,
                   rsq = rsq,
                   rmse = rmse)
      })
    }))), sp_model, by = "cluster", all = TRUE)
  } else if (check_mode == "spatial") {
    sp_model
  } else if (check_mode == "temporal") {
    as.data.frame(do.call(rbind, lapply(ts_model, function(m) {
      with(m, {
        data.frame(cluster = cluster,
                   rsq = rsq,
                   rmse = rmse)
      })
    })))
  }

  ## check clustering
  res_table <- within(m_data, {
    if (any(check_mode == "temporal")) {
      ts_check <- do.call(data.frame, lapply(cc$ts, function(f) {
        parse(text = f) |> eval()
      }))
      ts_check <- apply(ts_check, 1, function(r) {
        if (sum(r) == 0) {
          return(".")
        } else {
          strrep("*", sum(r))
        }
      })
    } else {
      ts_check <- NA
    }
    if (any(check_mode == "spatial")) {
      sp_check <- do.call(data.frame, lapply(cc$sp, function(f) {
        parse(text = f) |> eval()
      }))
      sp_check <- apply(sp_check, 1, function(r) {
        if (sum(r) == 0) {
          return(".")
        } else {
          strrep("*", sum(r))
        }
      })
    } else {
      sp_check <- NA
    }

  })

  ## return
  ret <- list(results = res_table,
              check_mode = check_mode,
              sp_model = local({
                if (any(check_mode == "spatial")) {
                  list(joincount = jc_data,
                       joincount_mc = jc_mc_data,
                       neighbor_list = nb,
                       weight_list = w)
                } else {
                  NA
                }
              }),
              ts_model = local({
                if (any(check_mode == "temporal")) {
                  ts_model
                } else {
                  NA
                }
              }),
              # sp_model = list(joincount = jc_data,
              #                 joincount_mc = jc_mc_data,
              #                 neighbor_list = nb,
              #                 weight_list = w),
              # ts_model = ts_model,
              geo = dplyr::left_join(x = sf::st_sf(cluster = cluster,
                                                   geometry = sf::st_transform(geo,
                                                                               crs = 4326)),
                                     y = subset(res_table,
                                                select = c(as.factor(cluster),
                                                           factor(sp_check, levels = c("***", "**", "*", ".")),
                                                           factor(ts_check, levels = c("***", "**", "*", ".")))),
                                     by = "cluster"),
              check_conditions = cc,
              data = x,
              cluster = cluster,
              neighbor_method = local({
                if (any(check_mode == "spatial")) {
                  deparse(bquote(do.call(.(neighbor_method),
                                         c(list(x = geo), .(neighbor_option)))))
                } else {
                  NA
                }
              }),
              # neighbor_method = deparse(bquote(do.call(.(neighbor_method),
              #                                          c(list(x = geo), .(neighbor_option))))),
              nsim = nsim)
  structure(ret, class = "stclust")
}

#' calc Root-Mean-Squared-Error (RMSE)
#'
#' 精度評価用に観測値と予測値との比較を行います。
#'
#' @param x 観測値を数値ベクトルで指定します。
#' @param y 予測値を数値ベクトルで指定します。
#' @examples
#' rmse(x = 1:10,
#'      y = rnorm(10))
rmse <- function(x = double(),
                 y = double()) {
  round(sqrt(mean((x - y) ^ 2, na.rm = TRUE)), digits = 2)
}


#' calc coefficient of determination
#'
#' R ^ 2値を計算します。
#'
#' @param x 観測値を数値ベクトルで指定します。
#' @param y 予測値を数値ベクトルで指定します。
#' @examples
#' rsq(x = 1:10,
#'    y = rnorm(10))
rsq <- function(x = double(),
                y = double()) {
  round(1 - mean((x - y) ^ 2, na.rm = TRUE) / var(x, na.rm = TRUE), digits = 2)
}

#' get num of sparseDFM model's factor.
#'
#' @param x set ts-data
get_num_factor <- function(x) {
  as.integer(gsub("^.*is *", "",
                  sparseDFM::tuneFactors(x, plot = FALSE)))
}

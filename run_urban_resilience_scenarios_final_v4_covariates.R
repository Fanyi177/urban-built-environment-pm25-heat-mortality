# ============================================================
# Urban environment scenario engine v4
# Long + short, Scenario 1/2/3, fixed PCA, bootstrap
# ============================================================

run_urban_resilience_scenarios_final_v4 <- function(
    data_list,
    time_scale = c("long", "short"),
    exposure_vars,
    urban_exposure,
    result_nonlinear_continuous_modification = NULL,
    future_projection_table = NULL,
    pca_projection_params = NULL,
    
    city_var = "city",
    province_var = "province",
    year_var = "dead_year",
    death_var = "n",
    pop_var = "total.population",
    
    gdp_quintile_var = "GDP_avg_quintile",
    # 必须使用从完整城市年度样本事先生成的城镇化率五分位，
    # 不应在每个分层文件或 bootstrap 样本内重新计算分位点。
    urban_quintile_var = "urbanization_rate_quintile",
    male_ratio_var = "Male_Ratio",
    age65_ratio_var = "popolation_age_65_Ratio",
    tmean_summer_var = "t2m_mean_summer",
    tmean_winter_var = "t2m_mean_winter",
    tsd_summer_var = "t2m_sd_summer",
    tsd_winter_var = "t2m_sd_winter",
    clim_zone_var = "climate_zone_name_code",
    
    dow_col = "dow",
    season_col = "season",
    rh_col = "lag03_rh",
    t2m_col = "lag03_t2m",
    time_col = "time",
    short_year_col = "year",
    ns_time_df = 7 * 8,
    k_smooth = 3,
    k_ti = 3,
    
    run_future = TRUE,
    run_future_short = FALSE,
    future_years = c(2030, 2050, 2070),
    
    run_threshold = TRUE,
    threshold_protective_q = 0.80,
    threshold_harmful_q = 0.20,
    slope_tol = 1e-4,
    
    run_policy_integrated_be = TRUE,
    run_compact_city = TRUE,
    run_complete_community = TRUE,
    run_active_mobility = TRUE,
    policy_enhancement_delta = 0.30,
    raw_indicator_lower_bound = 0,
    raw_indicator_upper_bound = 100,
    fixed_pca_clip = TRUE,
    
    compact_positive_vars = c(
      "Residential_per_totalArea",
      "Workplace_per_totalArea",
      "Commerce_services_per_totalArea",
      "Education_culture_per_totalArea",
      "Healthcare_per_totalArea",
      "Government_per_totalArea",
      "Leisure_tourism_per_totalArea",
      "Logistics_delivery_per_totalArea",
      "impervious",
      "Urban_public_transit_per_totalArea",
      "total_len_km_per_totalArea",
      "minor_len_km_per_totalArea"
    ),
    compact_negative_vars = character(0),
    
    complete_positive_vars = c(
      "Education_culture_per_totalArea",
      "Healthcare_per_totalArea",
      "Government_per_totalArea",
      "Commerce_services_per_totalArea",
      "Leisure_tourism_per_totalArea",
      "Logistics_delivery_per_totalArea",
      "Urban_public_transit_per_totalArea",
      "total_len_km_per_totalArea",
      "minor_len_km_per_totalArea",
      "forest_shrub_grassland",
      "water_wetland"
    ),
    complete_negative_vars = character(0),
    complete_residential_var = "Residential_per_totalArea",
    complete_workplace_var = "Workplace_per_totalArea",
    
    active_positive_vars = c(
      "Urban_public_transit_per_totalArea",
      "total_len_km_per_totalArea",
      "minor_len_km_per_totalArea",
      "Commerce_services_per_totalArea",
      "Education_culture_per_totalArea",
      "Healthcare_per_totalArea",
      "Leisure_tourism_per_totalArea",
      "forest_shrub_grassland",
      "water_wetland"
    ),
    active_negative_vars = c(
      "Motorized_transport_support_per_totalArea",
      "major_len_km_per_totalArea"
    ),
    
    n_boot = 300,
    cluster_var = "climate_zone_name_code",
    ntree = 300,
    seed = 2026,
    
    base_year = 2020,
    observed_years = NULL,
    death_reference_type = c("avg_annual_observed", "base_year_observed"),
    deaths_reference_base_year = NULL,
    annual_deaths_reference_fullpanel = NULL,
    
    output_dir = "result_project",
    output_city_snapshots = TRUE,
    show_message = TRUE
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(mgcv)
    library(splines)
    library(randomForest)
  })
  
  time_scale <- match.arg(time_scale)
  death_reference_type <- match.arg(death_reference_type)
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  safe_sd <- function(x) {
    out <- stats::sd(as.numeric(x), na.rm = TRUE)
    if (!is.finite(out) || out == 0) out <- 1
    out
  }
  
  safe_quan <- function(x, p) {
    as.numeric(stats::quantile(as.numeric(x), probs = p, na.rm = TRUE, names = FALSE, type = 8))
  }
  
  clip_bounds <- function(x, lower = raw_indicator_lower_bound, upper = raw_indicator_upper_bound) {
    pmin(pmax(as.numeric(x), lower), upper)
  }
  
  existing_vars <- function(df, vars) {
    intersect(vars, names(df))
  }
  
  sanitize_short_factors <- function(df) {
    factor_vars <- c(
      "dow",
      "season",
      ".year_for_exp_pred",
      "climate_zone_name_code"
    )
    
    factor_vars <- intersect(factor_vars, names(df))
    
    for (v in factor_vars) {
      df[[v]] <- factor(as.character(df[[v]]))
      attr(df[[v]], "contrasts") <- NULL
    }
    
    droplevels(df)
  }
  
  #
  sanitize_long_factors <- function(df) {
    factor_vars <- c(
      "GDP_avg_quintile",
      "urbanization_rate_quintile",
      "climate_zone_name_code",
      "dead_year_rf",
      "GDP_avg_quintile_rf",
      "urbanization_rate_quintile_rf",
      "climate_zone_name_code_rf"
    )
    
    factor_vars <- intersect(factor_vars, names(df))
    
    for (v in factor_vars) {
      df[[v]] <- factor(as.character(df[[v]]))
      attr(df[[v]], "contrasts") <- NULL
    }
    
    droplevels(df)
  }
  
  has_at_least_two_levels <- function(df, var) {
    if (is.null(df) || !(var %in% names(df))) return(TRUE)
    x <- df[[var]]
    x <- x[!is.na(x)]
    length(unique(as.character(x))) >= 2
  }
  
  make_dead_year_term <- function(df, default_k = 4) {
    if (is.null(df) || !("dead_year" %in% names(df))) {
      return(paste0("s(dead_year, k=", default_k, ")"))
    }
    
    yy <- unique(as.numeric(df$dead_year))
    yy <- yy[is.finite(yy)]
    nyr <- length(yy)
    
    if (nyr >= 5) {
      paste0("s(dead_year, k=", min(default_k, nyr - 1), ")")
    } else {
      "dead_year"
    }
  }
  
  
  
  safe_dev_expl <- function(gam_obj) {
    out <- tryCatch(summary(gam_obj)$dev.expl, error = function(e) NA_real_)
    as.numeric(out)
  }
  
  project_pca_score_fixed <- function(new_df, param, clip = TRUE) {
    vars <- param$vars
    miss <- setdiff(vars, names(new_df))
    if (length(miss) > 0) {
      stop("固定 PCA 投影失败，缺少变量: ", paste(miss, collapse = ", "))
    }
    
    X0 <- new_df %>%
      dplyr::select(all_of(vars)) %>%
      dplyr::mutate(across(everything(), ~ as.numeric(.)))
    X0 <- as.matrix(X0)
    
    center <- param$center[vars]
    scalev <- param$scale[vars]
    scalev[!is.finite(scalev) | scalev == 0] <- 1
    
    Xs <- sweep(X0, 2, center, "-")
    Xs <- sweep(Xs, 2, scalev, "/")
    
    loading <- param$loading[vars]
    raw <- as.numeric(Xs %*% loading)
    
    if (isTRUE(param$out_0_100)) {
      denom <- param$pc1_raw_max - param$pc1_raw_min
      
      if (!is.finite(denom) || denom == 0) {
        score <- rep((param$score_lower + param$score_upper) / 2, length(raw))
      } else {
        score <- (raw - param$pc1_raw_min) / denom *
          (param$score_upper - param$score_lower) + param$score_lower
      }
      
      if (isTRUE(clip)) {
        score <- pmin(pmax(score, param$score_lower), param$score_upper)
      }
    } else {
      score <- raw
    }
    
    score
  }
  
  add_fixed_pca_scores <- function(new_df, pca_projection_params, clip = TRUE) {
    out <- new_df
    blocks <- setdiff(names(pca_projection_params), "total_2stage")
    
    for (bb in blocks) {
      param <- pca_projection_params[[bb]]
      out[[param$score_col]] <- project_pca_score_fixed(out, param, clip = clip)
    }
    
    out
  }
  
  prepare_common_df <- function(dat) {
    dat <- as.data.frame(dat)
    
    if (time_scale == "long") {
      req <- c(
        city_var, province_var, year_var, death_var, pop_var,
        gdp_quintile_var, urban_quintile_var, male_ratio_var, age65_ratio_var,
        tmean_summer_var, tmean_winter_var,
        tsd_summer_var, tsd_winter_var, clim_zone_var
      )
      miss <- setdiff(req, names(dat))
      if (length(miss) > 0) {
        stop("长期数据缺少必要列: ", paste(miss, collapse = ", "),
             "。请先在完整城市年度数据中生成城镇化率五分位变量，并随原始文件/分层文件保存。")
      }
      
      out <- dat %>%
        mutate(
          .city = .data[[city_var]],
          .province = .data[[province_var]],
          n = as.numeric(.data[[death_var]]),
          dead_year = as.numeric(.data[[year_var]]),
          lpopulation = log(as.numeric(.data[[pop_var]])),
          GDP_avg_quintile = factor(.data[[gdp_quintile_var]]),
          urbanization_rate_quintile = factor(.data[[urban_quintile_var]]),
          Male_Ratio = as.numeric(.data[[male_ratio_var]]),
          popolation_age_65_Ratio = as.numeric(.data[[age65_ratio_var]]),
          t2m_mean_summer = as.numeric(.data[[tmean_summer_var]]),
          t2m_mean_winter = as.numeric(.data[[tmean_winter_var]]),
          t2m_sd_summer = as.numeric(.data[[tsd_summer_var]]),
          t2m_sd_winter = as.numeric(.data[[tsd_winter_var]]),
          climate_zone_name_code = factor(.data[[clim_zone_var]]),
          dead_year_rf = factor(.data[[year_var]]),
          GDP_avg_quintile_rf = factor(.data[[gdp_quintile_var]]),
          urbanization_rate_quintile_rf = factor(.data[[urban_quintile_var]]),
          climate_zone_name_code_rf = factor(.data[[clim_zone_var]])
        )
      
      return(out)
    }
    
    req <- c(
      city_var, province_var, death_var,
      dow_col, season_col, rh_col, t2m_col,
      time_col, clim_zone_var
    )
    miss <- setdiff(req, names(dat))
    if (length(miss) > 0) stop("短期数据缺少必要列: ", paste(miss, collapse = ", "))
    
    out <- dat %>%
      mutate(
        .city = .data[[city_var]],
        .province = .data[[province_var]],
        n = as.numeric(.data[[death_var]]),
        dow = .data[[dow_col]],
        season = .data[[season_col]],
        lag03_rh = as.numeric(.data[[rh_col]]),
        lag03_t2m = as.numeric(.data[[t2m_col]]),
        time = as.numeric(.data[[time_col]]),
        climate_zone_name_code = factor(.data[[clim_zone_var]])
      )
    
    if (short_year_col %in% names(out)) {
      out$.year_for_exp_pred <- factor(out[[short_year_col]])
    } else if (year_var %in% names(out)) {
      out$.year_for_exp_pred <- factor(out[[year_var]])
    } else {
      out$.year_for_exp_pred <- factor(floor((out$time - min(out$time, na.rm = TRUE)) / 365.25) + 1)
    }
    
    out
  }
  
  prep_modifier_df <- function(common_df, modifier_vec) {
    modifier_vec <- as.numeric(modifier_vec)
    mu <- mean(modifier_vec, na.rm = TRUE)
    sdv <- safe_sd(modifier_vec)
    
    out <- common_df %>%
      mutate(
        modifier_raw = modifier_vec,
        modifier_scaled = as.numeric((modifier_raw - mu) / sdv)
      )
    
    list(data = out, modifier_mean = mu, modifier_sd = sdv)
  }
  
  refresh_modifier_scaled <- function(df, mu, sdv) {
    df %>% mutate(modifier_scaled = as.numeric((modifier_raw - mu) / sdv))
  }
  
  compute_reference_deaths <- function(common_df) {
    if (!is.null(deaths_reference_base_year)) {
      base_ref <- as.numeric(deaths_reference_base_year)
    } else if (time_scale == "long" && "dead_year" %in% names(common_df)) {
      base_ref <- common_df %>%
        filter(dead_year == base_year) %>%
        summarise(x = sum(n, na.rm = TRUE)) %>%
        pull(x)
      if (length(base_ref) == 0 || !is.finite(base_ref)) base_ref <- NA_real_
    } else {
      base_ref <- sum(common_df$n, na.rm = TRUE)
    }
    
    if (!is.null(annual_deaths_reference_fullpanel)) {
      avg_ref <- as.numeric(annual_deaths_reference_fullpanel)
    } else if (time_scale == "long" && "dead_year" %in% names(common_df)) {
      yrs <- if (is.null(observed_years)) sort(unique(common_df$dead_year)) else sort(unique(as.numeric(observed_years)))
      avg_ref <- common_df %>%
        filter(dead_year %in% yrs) %>%
        group_by(dead_year) %>%
        summarise(deaths = sum(n, na.rm = TRUE), .groups = "drop") %>%
        summarise(x = mean(deaths, na.rm = TRUE)) %>%
        pull(x)
    } else {
      avg_ref <- sum(common_df$n, na.rm = TRUE)
    }
    
    ref_use <- if (death_reference_type == "base_year_observed") base_ref else avg_ref
    
    list(base_ref = base_ref, avg_ref = avg_ref, ref_use = ref_use)
  }
  
  fit_exposure_long_rf <- function(df_train, df_target, exposure_var) {
    rf_covars <- c(
      "modifier_raw",
      "dead_year_rf",
      "GDP_avg_quintile_rf",
      "urbanization_rate_quintile_rf",
      "t2m_mean_summer",
      "t2m_mean_winter",
      "t2m_sd_summer",
      "t2m_sd_winter",
      "climate_zone_name_code_rf"
    )
    
    rf_covars <- existing_vars(df_train, rf_covars)
    req <- unique(c(exposure_var, rf_covars))
    ok <- complete.cases(df_train[, req, drop = FALSE])
    
    if (sum(ok) < 10) return(NULL)
    
    df_fit <- df_train[ok, , drop = FALSE]
    fml <- as.formula(paste0(exposure_var, " ~ ", paste(rf_covars, collapse = " + ")))
    
    set.seed(seed)
    fit <- tryCatch(
      randomForest::randomForest(fml, data = df_fit, ntree = ntree),
      error = function(e) NULL
    )
    
    if (is.null(fit)) return(NULL)
    
    pred_fit <- predict(fit, newdata = df_fit)
    r2 <- tryCatch(cor(df_fit[[exposure_var]], pred_fit, use = "complete.obs")^2, error = function(e) NA_real_)
    rmse <- sqrt(mean((df_fit[[exposure_var]] - pred_fit)^2, na.rm = TRUE))
    
    pred_target <- tryCatch(
      predict(fit, newdata = df_target),
      error = function(e) rep(NA_real_, nrow(df_target))
    )
    
    out <- df_target
    out[[paste0(exposure_var, "_pred")]] <- as.numeric(pred_target)
    
    list(model = fit, predicted = out, model_type = "long_rf", r2 = r2, rmse = rmse)
  }
  
  fit_exposure_short_gam <- function(df_train, df_target, exposure_var) {
    df_train <- sanitize_short_factors(df_train)
    df_target <- sanitize_short_factors(df_target)
    
    req <- c(
      exposure_var,
      "modifier_raw",
      "lag03_t2m",
      "lag03_rh",
      "dow",
      "season",
      ".year_for_exp_pred",
      "climate_zone_name_code"
    )
    
    req <- existing_vars(df_train, req)
    ok <- complete.cases(df_train[, req, drop = FALSE])
    
    if (sum(ok) < 20) return(NULL)
    
    df_fit <- df_train[ok, , drop = FALSE]
    df_fit <- sanitize_short_factors(df_fit)
    
    fml <- as.formula(paste0(
      exposure_var, " ~ ",
      "s(modifier_raw, k=3) + ",
      "splines::ns(lag03_t2m, df=4) + ",
      "splines::ns(lag03_rh, df=3) + ",
      "dow + ",
      "season + ",
      ".year_for_exp_pred + ",
      "climate_zone_name_code"
    ))
    
    fit <- tryCatch(
      mgcv::bam(
        fml,
        data = df_fit,
        method = "fREML",
        discrete = TRUE,
        nthreads = 1,
        gc.level = 1,
        na.action = na.omit
      ),
      error = function(e) {
        message(
          "[SHORT EXPOSURE BAM ERROR] exposure=", exposure_var,
          " | error=", conditionMessage(e)
        )
        NULL
      }
    )
    
    if (is.null(fit)) return(NULL)
    
    pred_fit <- tryCatch(
      predict(fit, newdata = df_fit, type = "response"),
      error = function(e) rep(NA_real_, nrow(df_fit))
    )
    
    r2 <- tryCatch(
      cor(df_fit[[exposure_var]], pred_fit, use = "complete.obs")^2,
      error = function(e) NA_real_
    )
    
    rmse <- sqrt(mean((df_fit[[exposure_var]] - pred_fit)^2, na.rm = TRUE))
    
    pred_target <- tryCatch(
      predict(fit, newdata = df_target, type = "response"),
      error = function(e) {
        message(
          "[SHORT EXPOSURE PREDICT ERROR] exposure=", exposure_var,
          " | error=", conditionMessage(e)
        )
        rep(NA_real_, nrow(df_target))
      }
    )
    
    out <- df_target
    out[[paste0(exposure_var, "_pred")]] <- as.numeric(pred_target)
    
    list(
      model = fit,
      predicted = out,
      model_type = "short_bam",
      r2 = r2,
      rmse = rmse
    )
  }
  
  fit_exposure_model <- function(df_train, df_target, exposure_var) {
    if (time_scale == "long") {
      fit_exposure_long_rf(df_train, df_target, exposure_var)
    } else {
      fit_exposure_short_gam(df_train, df_target, exposure_var)
    }
  }
  
  add_exposure_prediction <- function(df_new, exposure_res, exposure_var) {
    if (is.null(exposure_res)) return(NULL)
    
    if (time_scale == "short") {
      df_new <- sanitize_short_factors(df_new)
    }
    
    pred <- if (inherits(exposure_res$model, "randomForest")) {
      tryCatch(
        predict(exposure_res$model, newdata = df_new),
        error = function(e) {
          message(
            "[EXPOSURE PREDICT ERROR] exposure=", exposure_var,
            " | error=", conditionMessage(e)
          )
          rep(NA_real_, nrow(df_new))
        }
      )
    } else {
      tryCatch(
        predict(exposure_res$model, newdata = df_new, type = "response"),
        error = function(e) {
          message(
            "[EXPOSURE PREDICT ERROR] exposure=", exposure_var,
            " | error=", conditionMessage(e)
          )
          rep(NA_real_, nrow(df_new))
        }
      )
    }
    
    df_new[[paste0(exposure_var, "_pred")]] <- as.numeric(pred)
    
    df_new
  }
  
  build_long_formula_linear <- function(exposure_var, df = NULL) {
    ep <- paste0(exposure_var, "_pred")
    
    rhs <- c(
      paste0(ep, " * modifier_scaled"),
      make_dead_year_term(df, default_k = 4)
    )
    
    if (has_at_least_two_levels(df, "GDP_avg_quintile")) {
      rhs <- c(rhs, "as.factor(GDP_avg_quintile)")
    }
    if (has_at_least_two_levels(df, "urbanization_rate_quintile")) {
      rhs <- c(rhs, "as.factor(urbanization_rate_quintile)")
    }

    rhs <- c(
      rhs,
      "Male_Ratio",
      "popolation_age_65_Ratio",
      "offset(lpopulation)",
      "t2m_mean_summer",
      "t2m_mean_winter",
      "t2m_sd_summer",
      "t2m_sd_winter"
    )
    
    if (has_at_least_two_levels(df, "climate_zone_name_code")) {
      rhs <- c(rhs, "s(climate_zone_name_code, bs='re')")
    }
    
    as.formula(paste0(
      death_var, " ~ ",
      paste(rhs, collapse = " + ")
    ))
  }
  build_long_formula_nonlin <- function(exposure_var, df = NULL) {
    ep <- paste0(exposure_var, "_pred")
    
    rhs <- c(
      paste0("s(", ep, ", k=4)"),
      paste0("s(", ep, ", by = modifier_scaled, k=4)"),
      "modifier_scaled",
      make_dead_year_term(df, default_k = 4)
    )
    
    if (has_at_least_two_levels(df, "GDP_avg_quintile")) {
      rhs <- c(rhs, "as.factor(GDP_avg_quintile)")
    }
    if (has_at_least_two_levels(df, "urbanization_rate_quintile")) {
      rhs <- c(rhs, "as.factor(urbanization_rate_quintile)")
    }

    rhs <- c(
      rhs,
      "Male_Ratio",
      "popolation_age_65_Ratio",
      "offset(lpopulation)",
      "t2m_mean_summer",
      "t2m_mean_winter",
      "t2m_sd_summer",
      "t2m_sd_winter"
    )
    
    if (has_at_least_two_levels(df, "climate_zone_name_code")) {
      rhs <- c(rhs, "s(climate_zone_name_code, bs='re')")
    }
    
    as.formula(paste0(
      death_var, " ~ ",
      paste(rhs, collapse = " + ")
    ))
  }
  
  build_short_formula_linear <- function(exposure_var) {
    ep <- paste0(exposure_var, "_pred")
    
    as.formula(paste0(
      death_var, " ~ ",
      ep, " * modifier_scaled + ",
      "as.factor(dow) + ",
      "as.factor(season) + ",
      "splines::ns(lag03_rh, 3) + ",
      "splines::ns(lag03_t2m, df=6) + ",
      "splines::ns(time, df=", ns_time_df, ") + ",
      "as.factor(climate_zone_name_code)"
    ))
  }
  
  build_short_formula_nonlin <- function(exposure_var) {
    ep <- paste0(exposure_var, "_pred")
    
    as.formula(paste0(
      death_var, " ~ ",
      "s(", ep, ", k=", k_smooth, ") + ",
      "s(", ep, ", by = modifier_scaled, k=", k_smooth, ") + ",
      "modifier_scaled + ",
      "as.factor(dow) + ",
      "as.factor(season) + ",
      "splines::ns(lag03_rh, 3) + ",
      "splines::ns(lag03_t2m, df=6) + ",
      "splines::ns(time, df=", ns_time_df, ") + ",
      "as.factor(climate_zone_name_code)"
    ))
  }
  
  fit_mortality_models <- function(df, exposure_var) {
    if (time_scale == "short") {
      df <- sanitize_short_factors(df)
    }
    
    if (time_scale == "long") {
      df <- sanitize_long_factors(df)
    }
    
    forms <- if (time_scale == "long") {
      list(
        Linear = build_long_formula_linear(exposure_var, df = df),
        Nonlinear = build_long_formula_nonlin(exposure_var, df = df)
      )
    } else {
      list(
        Linear = build_short_formula_linear(exposure_var),  
        Nonlinear = build_short_formula_nonlin(exposure_var)
      )
    }
    
    out <- list()
    
    for (nm in names(forms)) {
      out[[nm]] <- tryCatch(
        {
          if (time_scale == "short") {
            mgcv::bam(
              forms[[nm]],
              family = quasipoisson(link = "log"),
              data = df,
              method = "fREML",
              discrete = TRUE,
              gc.level = 1,
              nthreads = 1,
              na.action = na.omit
            )
          } else {
            mgcv::gam(
              forms[[nm]],
              family = quasipoisson(link = "log"),
              data = df,
              method = "REML",
              na.action = na.omit
            )
          }
        },
        error = function(e) {
          message(
            "[MORTALITY MODEL ERROR] time_scale=", time_scale,
            " | exposure=", exposure_var,
            " | model=", nm,
            " | error=", conditionMessage(e)
          )
          NULL
        }
      )
    }
    
    out <- out[!vapply(out, is.null, logical(1))]
    if (length(out) == 0) return(NULL)
    
    out
  }
  #
  make_pathway_counterfactuals <- function(df_base, df_cf_full, exposure_var) {
    ep <- paste0(exposure_var, "_pred")
    
    if (!(ep %in% names(df_base)) || !(ep %in% names(df_cf_full))) {
      return(list(full = df_cf_full))
    }
    
    # 1) full: U 改，E 也改
    df_full <- df_cf_full
    
    # 2) modifier_only: 只改 U，E 保持 baseline
    df_modifier_only <- df_cf_full
    df_modifier_only[[ep]] <- df_base[[ep]]
    
    # 3) exposure_only: 只改 E，U 保持 baseline
    df_exposure_only <- df_base
    df_exposure_only[[ep]] <- df_cf_full[[ep]]
    
    out <- list(
      full = df_full,
      modifier_only = df_modifier_only,
      exposure_only = df_exposure_only
    )
    
    if (time_scale == "short") {
      out <- lapply(out, sanitize_short_factors)
    }
    
    out
  }
  #
  
  predict_totals <- function(df_new, gam_list) {
    if (time_scale == "short") {
      df_new <- sanitize_short_factors(df_new)
    }
    
    vapply(gam_list, function(gm) {
      sum(predict(gm, newdata = df_new, type = "response"), na.rm = TRUE)
    }, numeric(1))
  }
  #
  make_point_rows <- function(
    dname,
    analysis_type,
    urban_var,
    exposure_var,
    scenario,
    baseline_tot,
    cf_tot,
    extra = data.frame()
  ) {
    common_models <- intersect(names(baseline_tot), names(cf_tot))
    
    out <- data.frame(
      disease = dname,
      time_scale = time_scale,
      analysis_type = analysis_type,
      urban = urban_var,
      exposure = exposure_var,
      scenario = scenario,
      model = common_models,
      baseline_total = as.numeric(baseline_tot[common_models]),
      scenario_total = as.numeric(cf_tot[common_models]),
      death_change = as.numeric(cf_tot[common_models] - baseline_tot[common_models]),
      AF = as.numeric((cf_tot[common_models] - baseline_tot[common_models]) / baseline_tot[common_models]),
      AF_percent = as.numeric((cf_tot[common_models] - baseline_tot[common_models]) / baseline_tot[common_models] * 100),
      stringsAsFactors = FALSE
    )
    
    if (nrow(extra) == 1) {
      out <- bind_cols(out, extra[rep(1, nrow(out)), , drop = FALSE])
    }
    
    out
  }
  
  improve_positive_vars <- function(df, vars, delta) {
    out <- df
    
    for (v in existing_vars(out, vars)) {
      out[[v]] <- clip_bounds(out[[v]] + delta * (raw_indicator_upper_bound - out[[v]]))
    }
    
    out
  }
  
  reduce_negative_vars <- function(df, vars, delta) {
    out <- df
    
    for (v in existing_vars(out, vars)) {
      out[[v]] <- clip_bounds(out[[v]] - delta * (out[[v]] - raw_indicator_lower_bound))
    }
    
    out
  }
  
  improve_res_work_balance <- function(df, delta) {
    out <- df
    rvar <- complete_residential_var
    wvar <- complete_workplace_var
    
    if (!all(c(rvar, wvar) %in% names(out))) return(out)
    
    r <- as.numeric(out[[rvar]])
    w <- as.numeric(out[[wvar]])
    
    r2 <- r
    w2 <- w
    
    idx_r_low <- is.finite(r) & is.finite(w) & r < w
    idx_w_low <- is.finite(r) & is.finite(w) & w < r
    
    r2[idx_r_low] <- r[idx_r_low] + delta * (w[idx_r_low] - r[idx_r_low])
    w2[idx_w_low] <- w[idx_w_low] + delta * (r[idx_w_low] - w[idx_w_low])
    
    out[[rvar]] <- clip_bounds(r2)
    out[[wvar]] <- clip_bounds(w2)
    
    out
  }
  
  build_policy_counterfactual_raw <- function(df, policy_framework, delta) {
    out <- df
    
    if (policy_framework == "compact_city") {
      out <- improve_positive_vars(out, compact_positive_vars, delta)
      out <- reduce_negative_vars(out, compact_negative_vars, delta)
    }
    
    if (policy_framework == "complete_community") {
      out <- improve_positive_vars(out, complete_positive_vars, delta)
      out <- reduce_negative_vars(out, complete_negative_vars, delta)
      out <- improve_res_work_balance(out, delta)
    }
    
    if (policy_framework == "active_mobility_supportive") {
      out <- improve_positive_vars(out, active_positive_vars, delta)
      out <- reduce_negative_vars(out, active_negative_vars, delta)
    }
    
    add_fixed_pca_scores(out, pca_projection_params, clip = fixed_pca_clip)
  }
  
  get_policy_registry <- function() {
    out <- data.frame(policy_framework = character(), stringsAsFactors = FALSE)
    
    if (isTRUE(run_compact_city)) {
      out <- bind_rows(out, data.frame(policy_framework = "compact_city"))
    }
    
    if (isTRUE(run_complete_community)) {
      out <- bind_rows(out, data.frame(policy_framework = "complete_community"))
    }
    
    if (isTRUE(run_active_mobility)) {
      out <- bind_rows(out, data.frame(policy_framework = "active_mobility_supportive"))
    }
    
    out
  }
  
  build_policy_raw_snapshot <- function(df_base, df_cf, dname, policy_framework, scenario, delta) {
    vars_policy <- switch(
      policy_framework,
      compact_city = unique(c(compact_positive_vars, compact_negative_vars)),
      complete_community = unique(c(
        complete_positive_vars,
        complete_negative_vars,
        complete_residential_var,
        complete_workplace_var
      )),
      active_mobility_supportive = unique(c(active_positive_vars, active_negative_vars)),
      character(0)
    )
    
    vars_policy <- existing_vars(df_base, vars_policy)
    if (length(vars_policy) == 0) return(data.frame())
    
    bind_rows(lapply(vars_policy, function(v) {
      data.frame(
        disease = dname,
        time_scale = time_scale,
        policy_framework = policy_framework,
        scenario = scenario,
        delta = delta,
        raw_indicator = v,
        baseline_mean = mean(as.numeric(df_base[[v]]), na.rm = TRUE),
        scenario_mean = mean(as.numeric(df_cf[[v]]), na.rm = TRUE),
        mean_change = mean(as.numeric(df_cf[[v]]) - as.numeric(df_base[[v]]), na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
  }
  
  build_policy_score_snapshot <- function(df_base, df_cf, dname, policy_framework, scenario, delta) {
    vars_score <- existing_vars(df_base, urban_exposure)
    if (length(vars_score) == 0) return(data.frame())
    
    bind_rows(lapply(vars_score, function(v) {
      data.frame(
        disease = dname,
        time_scale = time_scale,
        policy_framework = policy_framework,
        scenario = scenario,
        delta = delta,
        urban = v,
        baseline_mean = mean(as.numeric(df_base[[v]]), na.rm = TRUE),
        scenario_mean = mean(as.numeric(df_cf[[v]]), na.rm = TRUE),
        mean_change = mean(as.numeric(df_cf[[v]]) - as.numeric(df_base[[v]]), na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
  }
  build_exposure_prediction_snapshot <- function(df_base, df_cf, dname,
                                                 analysis_type,
                                                 policy_framework,
                                                 scenario,
                                                 delta,
                                                 urban_var,
                                                 exposure_var) {
    ep <- paste0(exposure_var, "_pred")
    
    if (!(ep %in% names(df_base)) || !(ep %in% names(df_cf))) {
      return(data.frame())
    }
    
    data.frame(
      disease = dname,
      time_scale = time_scale,
      analysis_type = analysis_type,
      policy_framework = policy_framework,
      scenario = scenario,
      delta = delta,
      urban = urban_var,
      exposure = exposure_var,
      exposure_pred_baseline_mean = mean(as.numeric(df_base[[ep]]), na.rm = TRUE),
      exposure_pred_scenario_mean = mean(as.numeric(df_cf[[ep]]), na.rm = TRUE),
      exposure_pred_mean_change = mean(as.numeric(df_cf[[ep]]) - as.numeric(df_base[[ep]]), na.rm = TRUE),
      exposure_pred_baseline_p25 = as.numeric(quantile(as.numeric(df_base[[ep]]), 0.25, na.rm = TRUE)),
      exposure_pred_scenario_p25 = as.numeric(quantile(as.numeric(df_cf[[ep]]), 0.25, na.rm = TRUE)),
      exposure_pred_baseline_p50 = as.numeric(quantile(as.numeric(df_base[[ep]]), 0.50, na.rm = TRUE)),
      exposure_pred_scenario_p50 = as.numeric(quantile(as.numeric(df_cf[[ep]]), 0.50, na.rm = TRUE)),
      exposure_pred_baseline_p75 = as.numeric(quantile(as.numeric(df_base[[ep]]), 0.75, na.rm = TRUE)),
      exposure_pred_scenario_p75 = as.numeric(quantile(as.numeric(df_cf[[ep]]), 0.75, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }
  #
  # ============================================================
  # City-level snapshot helpers for later economic valuation
  # ============================================================
  
  get_n_years_observed <- function(df) {
    if (!is.null(observed_years)) {
      yy <- unique(as.numeric(observed_years))
      yy <- yy[is.finite(yy)]
      if (length(yy) > 0) return(length(yy))
    }
    
    year_candidates <- unique(c(
      if (time_scale == "short") short_year_col else year_var,
      "year", "dead_year"
    ))
    
    for (v in year_candidates) {
      if (v %in% names(df)) {
        yy <- unique(as.numeric(df[[v]]))
        yy <- yy[is.finite(yy)]
        if (length(yy) > 0) return(length(yy))
      }
    }
    
    1
  }
  
  prep_city_id_cols <- function(df) {
    out <- df
    
    if (!(".city" %in% names(out)) && city_var %in% names(out)) {
      out$.city <- out[[city_var]]
    }
    
    if (!(".province" %in% names(out)) && province_var %in% names(out)) {
      out$.province <- out[[province_var]]
    }
    
    out
  }
  
  get_city_group_cols <- function(df) {
    cols <- c(".province", ".city", "province_code", "city_code")
    cols <- intersect(cols, names(df))
    
    if (!(".city" %in% cols)) {
      stop("city-level snapshot failed: .city/city not found")
    }
    
    cols
  }
  
  build_scenario_city_change_snapshot <- function(df_base,
                                                  df_cf,
                                                  gam_list,
                                                  dname,
                                                  analysis_type,
                                                  policy_framework,
                                                  scenario,
                                                  delta,
                                                  pathway,
                                                  urban_var,
                                                  exposure_var) {
    if (!isTRUE(output_city_snapshots)) return(data.frame())
    if (is.null(df_base) || is.null(df_cf) || is.null(gam_list)) return(data.frame())
    
    ep <- paste0(exposure_var, "_pred")
    if (!(ep %in% names(df_base)) || !(ep %in% names(df_cf))) return(data.frame())
    
    df_base <- prep_city_id_cols(df_base)
    df_cf <- prep_city_id_cols(df_cf)
    
    if (time_scale == "short") {
      df_base <- sanitize_short_factors(df_base)
      df_cf <- sanitize_short_factors(df_cf)
    }
    
    id_cols <- get_city_group_cols(df_base)
    id_cols <- intersect(id_cols, names(df_cf))
    
    n_years <- get_n_years_observed(df_base)
    if (!is.finite(n_years) || n_years <= 0) n_years <- 1
    
    dplyr::bind_rows(lapply(names(gam_list), function(model_nm) {
      gm <- gam_list[[model_nm]]
      
      b <- df_base
      cdf <- df_cf
      
      b$.pred_death <- tryCatch(
        as.numeric(predict(gm, newdata = b, type = "response")),
        error = function(e) rep(NA_real_, nrow(b))
      )
      
      cdf$.pred_death <- tryCatch(
        as.numeric(predict(gm, newdata = cdf, type = "response")),
        error = function(e) rep(NA_real_, nrow(cdf))
      )
      
      b_sum <- b %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(id_cols))) %>%
        dplyr::summarise(
          modifier_baseline = mean(as.numeric(modifier_raw), na.rm = TRUE),
          exposure_pred_baseline_mean = mean(as.numeric(.data[[ep]]), na.rm = TRUE),
          deaths_baseline_pred = sum(.pred_death, na.rm = TRUE),
          n_rows_baseline = dplyr::n(),
          .groups = "drop"
        )
      
      c_sum <- cdf %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(id_cols))) %>%
        dplyr::summarise(
          modifier_scenario = mean(as.numeric(modifier_raw), na.rm = TRUE),
          exposure_pred_scenario_mean = mean(as.numeric(.data[[ep]]), na.rm = TRUE),
          deaths_scenario_pred = sum(.pred_death, na.rm = TRUE),
          n_rows_scenario = dplyr::n(),
          .groups = "drop"
        )
      
      out <- b_sum %>%
        dplyr::full_join(c_sum, by = id_cols) %>%
        dplyr::mutate(
          disease = dname,
          time_scale = time_scale,
          analysis_type = analysis_type,
          policy_framework = policy_framework,
          scenario = scenario,
          delta = delta,
          pathway = pathway,
          urban = urban_var,
          exposure = exposure_var,
          model = model_nm,
          
          modifier_delta = modifier_scenario - modifier_baseline,
          modifier_delta_percent = dplyr::if_else(
            is.finite(modifier_baseline) & modifier_baseline != 0,
            modifier_delta / modifier_baseline * 100,
            NA_real_
          ),
          
          exposure_pred_delta_mean =
            exposure_pred_scenario_mean - exposure_pred_baseline_mean,
          
          death_change =
            deaths_scenario_pred - deaths_baseline_pred,
          
          AF = dplyr::if_else(
            is.finite(deaths_baseline_pred) & deaths_baseline_pred != 0,
            death_change / deaths_baseline_pred,
            NA_real_
          ),
          
          AF_percent = AF * 100,
          n_years_observed = n_years,
          annualized_death_change = death_change / n_years
        ) %>%
        dplyr::select(
          disease, time_scale, analysis_type, policy_framework,
          scenario, delta, pathway, urban, exposure, model,
          dplyr::all_of(id_cols),
          modifier_baseline, modifier_scenario,
          modifier_delta, modifier_delta_percent,
          exposure_pred_baseline_mean,
          exposure_pred_scenario_mean,
          exposure_pred_delta_mean,
          deaths_baseline_pred,
          deaths_scenario_pred,
          death_change,
          AF,
          AF_percent,
          n_years_observed,
          annualized_death_change,
          n_rows_baseline,
          n_rows_scenario
        )
      
      out
    }))
  }
  
  build_policy_raw_indicator_city_snapshot <- function(df_base,
                                                       df_cf,
                                                       dname,
                                                       policy_framework,
                                                       scenario,
                                                       delta) {
    if (!isTRUE(output_city_snapshots)) return(data.frame())
    
    vars_policy <- switch(
      policy_framework,
      compact_city = unique(c(compact_positive_vars, compact_negative_vars)),
      complete_community = unique(c(
        complete_positive_vars,
        complete_negative_vars,
        complete_residential_var,
        complete_workplace_var
      )),
      active_mobility_supportive = unique(c(active_positive_vars, active_negative_vars)),
      character(0)
    )
    
    vars_policy <- existing_vars(df_base, vars_policy)
    if (length(vars_policy) == 0) return(data.frame())
    
    df_base <- prep_city_id_cols(df_base)
    df_cf <- prep_city_id_cols(df_cf)
    
    id_cols <- get_city_group_cols(df_base)
    id_cols <- intersect(id_cols, names(df_cf))
    
    dplyr::bind_rows(lapply(vars_policy, function(v) {
      b <- df_base %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(id_cols))) %>%
        dplyr::summarise(
          indicator_baseline_std = mean(as.numeric(.data[[v]]), na.rm = TRUE),
          .groups = "drop"
        )
      
      cdf <- df_cf %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(id_cols))) %>%
        dplyr::summarise(
          indicator_scenario_std = mean(as.numeric(.data[[v]]), na.rm = TRUE),
          .groups = "drop"
        )
      
      direction <- dplyr::case_when(
        v %in% c(compact_positive_vars, complete_positive_vars, active_positive_vars) ~
          "positive_toward_upper",
        v %in% c(compact_negative_vars, complete_negative_vars, active_negative_vars) ~
          "negative_toward_lower",
        policy_framework == "complete_community" &&
          v %in% c(complete_residential_var, complete_workplace_var) ~
          "balance_improved",
        TRUE ~ "unspecified"
      )
      
      b %>%
        dplyr::full_join(cdf, by = id_cols) %>%
        dplyr::mutate(
          disease = dname,
          time_scale = time_scale,
          policy_framework = policy_framework,
          scenario = scenario,
          delta = delta,
          raw_indicator = v,
          direction = direction,
          indicator_delta_std =
            indicator_scenario_std - indicator_baseline_std,
          indicator_delta_percent_std = dplyr::if_else(
            is.finite(indicator_baseline_std) & indicator_baseline_std != 0,
            indicator_delta_std / indicator_baseline_std * 100,
            NA_real_
          )
        ) %>%
        dplyr::select(
          disease, time_scale, policy_framework, scenario, delta,
          dplyr::all_of(id_cols),
          raw_indicator, direction,
          indicator_baseline_std,
          indicator_scenario_std,
          indicator_delta_std,
          indicator_delta_percent_std
        )
    }))
  }
  
  build_policy_score_city_snapshot <- function(df_base,
                                               df_cf,
                                               dname,
                                               policy_framework,
                                               scenario,
                                               delta) {
    if (!isTRUE(output_city_snapshots)) return(data.frame())
    
    vars_score <- existing_vars(df_base, urban_exposure)
    if (length(vars_score) == 0) return(data.frame())
    
    df_base <- prep_city_id_cols(df_base)
    df_cf <- prep_city_id_cols(df_cf)
    
    id_cols <- get_city_group_cols(df_base)
    id_cols <- intersect(id_cols, names(df_cf))
    
    dplyr::bind_rows(lapply(vars_score, function(v) {
      b <- df_base %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(id_cols))) %>%
        dplyr::summarise(
          score_baseline = mean(as.numeric(.data[[v]]), na.rm = TRUE),
          .groups = "drop"
        )
      
      cdf <- df_cf %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(id_cols))) %>%
        dplyr::summarise(
          score_scenario = mean(as.numeric(.data[[v]]), na.rm = TRUE),
          .groups = "drop"
        )
      
      b %>%
        dplyr::full_join(cdf, by = id_cols) %>%
        dplyr::mutate(
          disease = dname,
          time_scale = time_scale,
          policy_framework = policy_framework,
          scenario = scenario,
          delta = delta,
          score_var = v,
          score_delta = score_scenario - score_baseline,
          score_delta_percent = dplyr::if_else(
            is.finite(score_baseline) & score_baseline != 0,
            score_delta / score_baseline * 100,
            NA_real_
          )
        ) %>%
        dplyr::select(
          disease, time_scale, policy_framework, scenario, delta,
          dplyr::all_of(id_cols),
          score_var,
          score_baseline,
          score_scenario,
          score_delta,
          score_delta_percent
        )
    }))
  }
  
  #
  
  smooth_curve_y <- function(curve_df) {
    curve_df <- curve_df %>%
      filter(!is.na(x)) %>%
      arrange(x)
    
    y_col <- if ("RR" %in% names(curve_df)) {
      "RR"
    } else if ("beta_per1" %in% names(curve_df)) {
      "beta_per1"
    } else {
      NULL
    }
    
    if (is.null(y_col)) return(NULL)
    
    x <- as.numeric(curve_df$x)
    y <- as.numeric(curve_df[[y_col]])
    
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]
    y <- y[ok]
    
    if (length(x) < 5) return(NULL)
    
    lo <- try(stats::loess(y ~ x, span = 0.5), silent = TRUE)
    y_sm <- if (inherits(lo, "try-error")) y else as.numeric(predict(lo, newdata = data.frame(x = x)))
    
    if (any(!is.finite(y_sm))) y_sm <- y
    
    list(x = x, y = y_sm)
  }
  
  classify_threshold_rule <- function(curve_df, observed_modifier) {
    cy <- smooth_curve_y(curve_df)
    if (is.null(cy)) return(NULL)
    
    x <- cy$x
    y <- cy$y
    
    dy <- diff(y)
    dy_sign <- ifelse(abs(dy) < slope_tol, 0, sign(dy))
    s <- dy_sign[dy_sign != 0]
    
    if (length(s) == 0) return(NULL)
    
    changes <- if (length(s) <= 1) 0 else sum(diff(s) != 0)
    first_s <- s[1]
    last_s <- s[length(s)]
    
    p20 <- safe_quan(observed_modifier, threshold_harmful_q)
    p80 <- safe_quan(observed_modifier, threshold_protective_q)
    
    y_at <- function(xx) {
      as.numeric(stats::approx(x, y, xout = xx, rule = 2)$y)
    }
    
    if (all(s < 0)) {
      return(data.frame(
        shape = "decreasing",
        threshold_type = "decreasing_raise_to_p80",
        action = "raise_to_threshold",
        scenario_x = p80,
        threshold_x = p80,
        RR_p20 = y_at(p20),
        RR_p80 = y_at(p80),
        stringsAsFactors = FALSE
      ))
    }
    
    if (all(s > 0)) {
      return(data.frame(
        shape = "increasing",
        threshold_type = "increasing_cap_to_p20",
        action = "cap_to_threshold",
        scenario_x = p20,
        threshold_x = p20,
        RR_p20 = y_at(p20),
        RR_p80 = y_at(p80),
        stringsAsFactors = FALSE
      ))
    }
    
    if (changes >= 1 && first_s < 0 && last_s > 0) {
      i_min <- which.min(y)
      
      return(data.frame(
        shape = "u_shape",
        threshold_type = "u_shape_minimum",
        action = "set_to_threshold",
        scenario_x = x[i_min],
        threshold_x = x[i_min],
        RR_p20 = y_at(p20),
        RR_p80 = y_at(p80),
        stringsAsFactors = FALSE
      ))
    }
    
    if (changes >= 1 && first_s > 0 && last_s < 0) {
      rr20 <- y_at(p20)
      rr80 <- y_at(p80)
      
      if (is.finite(rr20) && (!is.finite(rr80) || rr20 < rr80)) {
        return(data.frame(
          shape = "inverted_u",
          threshold_type = "inverted_u_low_side_p20",
          action = "cap_to_threshold",
          scenario_x = p20,
          threshold_x = x[which.max(y)],
          RR_p20 = rr20,
          RR_p80 = rr80,
          stringsAsFactors = FALSE
        ))
      } else {
        return(data.frame(
          shape = "inverted_u",
          threshold_type = "inverted_u_low_side_p80",
          action = "raise_to_threshold",
          scenario_x = p80,
          threshold_x = x[which.max(y)],
          RR_p20 = rr20,
          RR_p80 = rr80,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    NULL
  }
  
  apply_threshold_to_modifier <- function(df, rule) {
    out <- df
    x0 <- as.numeric(out$modifier_raw)
    sx <- as.numeric(rule$scenario_x[1])
    action <- as.character(rule$action[1])
    
    if (action == "raise_to_threshold") {
      out$modifier_raw <- pmax(x0, sx)
    } else if (action == "cap_to_threshold") {
      out$modifier_raw <- pmin(x0, sx)
    } else if (action == "set_to_threshold") {
      out$modifier_raw <- rep(sx, length(x0))
    }
    
    out
  }
  
  filter_curve_for_key <- function(mod_df, dname, exposure_var, urban_var) {
    if (is.null(mod_df) || nrow(mod_df) == 0) return(data.frame())
    
    out0 <- as.data.frame(mod_df)
    out <- out0
    
    if ("disease" %in% names(out)) {
      tmp <- out %>% filter(disease == dname)
      if (nrow(tmp) > 0) out <- tmp
    }
    
    if ("urban" %in% names(out)) {
      out <- out %>% filter(urban == urban_var)
    }
    
    if ("exposure_run" %in% names(out)) {
      out <- out %>% filter(exposure_run == exposure_var)
    } else if ("exposure_raw" %in% names(out)) {
      out <- out %>% filter(exposure_raw == exposure_var)
    } else if ("exposure" %in% names(out)) {
      tmp <- out %>% filter(exposure == exposure_var)
      if (nrow(tmp) > 0) out <- tmp
    }
    
    out
  }
  
  build_future_counterfactual <- function(df, future_year, proj_tbl, modifier_mean, modifier_sd) {
    if (time_scale == "long") {
      base <- df %>%
        group_by(.city) %>%
        slice_max(order_by = dead_year, n = 1, with_ties = FALSE) %>%
        ungroup()
    } else {
      base <- df
    }
    
    if (!is.null(proj_tbl) && nrow(proj_tbl) > 0) {
      pt <- proj_tbl
      
      if (!(".city" %in% names(pt)) && "city" %in% names(pt)) pt$.city <- pt$city
      
      pt <- pt %>%
        dplyr::select(.city, value_pred) %>%
        distinct(.city, .keep_all = TRUE)
      
      base <- base %>%
        left_join(pt, by = ".city") %>%
        mutate(modifier_raw = ifelse(!is.na(value_pred), value_pred, modifier_raw)) %>%
        dplyr::select(-any_of("value_pred"))
    }
    
    if (time_scale == "long") {
      base <- refresh_modifier_scaled(base, modifier_mean, modifier_sd)
    }
    
    base
  }
  
  bootstrap_sample <- function(df) {
    if (!is.null(cluster_var) && cluster_var %in% names(df)) {
      vals <- unique(df[[cluster_var]])
      vals <- vals[!is.na(vals)]
      
      if (length(vals) >= 2) {
        sampled <- sample(vals, length(vals), replace = TRUE)
        return(bind_rows(lapply(sampled, function(v) df[df[[cluster_var]] == v, , drop = FALSE])))
      }
    }
    
    if (".province" %in% names(df)) {
      vals <- unique(df$.province)
      vals <- vals[!is.na(vals)]
      
      if (length(vals) >= 2) {
        sampled <- sample(vals, length(vals), replace = TRUE)
        return(bind_rows(lapply(sampled, function(v) df[df$.province == v, , drop = FALSE])))
      }
    }
    
    if (".city" %in% names(df)) {
      vals <- unique(df$.city)
      vals <- vals[!is.na(vals)]
      
      if (length(vals) >= 2) {
        sampled <- sample(vals, length(vals), replace = TRUE)
        return(bind_rows(lapply(sampled, function(v) df[df$.city == v, , drop = FALSE])))
      }
    }
    
    df[sample(seq_len(nrow(df)), nrow(df), replace = TRUE), , drop = FALSE]
  }
  
  summarise_boot <- function(boot_df, reference_df) {
    if (nrow(boot_df) == 0) return(data.frame())
    
    if (!("pathway" %in% names(boot_df))) {
      boot_df$pathway <- "full"
    }
    
    boot_df <- boot_df %>%
      left_join(
        reference_df %>% dplyr::select(disease, ReferenceDeaths, death_reference_type),
        by = "disease"
      ) %>%
      mutate(
        annual_deaths_change_Linear = AF_Linear * ReferenceDeaths,
        annual_deaths_change_Nonlinear = AF_Nonlinear * ReferenceDeaths
      )
    
    boot_df %>%
      group_by(
        disease, time_scale, analysis_type, policy_framework, pathway,
        urban, exposure, scenario, delta,
        death_reference_type, ReferenceDeaths
      ) %>%
      summarise(
        n_boot_valid_Linear = sum(!is.na(AF_Linear)),
        AF_Linear_mean = mean(AF_Linear, na.rm = TRUE),
        AF_Linear_low = quantile(AF_Linear, 0.025, na.rm = TRUE),
        AF_Linear_high = quantile(AF_Linear, 0.975, na.rm = TRUE),
        death_change_Linear_mean = mean(death_change_Linear, na.rm = TRUE),
        death_change_Linear_low = quantile(death_change_Linear, 0.025, na.rm = TRUE),
        death_change_Linear_high = quantile(death_change_Linear, 0.975, na.rm = TRUE),
        annual_deaths_change_Linear_mean = mean(annual_deaths_change_Linear, na.rm = TRUE),
        annual_deaths_change_Linear_low = quantile(annual_deaths_change_Linear, 0.025, na.rm = TRUE),
        annual_deaths_change_Linear_high = quantile(annual_deaths_change_Linear, 0.975, na.rm = TRUE),
        
        n_boot_valid_Nonlinear = sum(!is.na(AF_Nonlinear)),
        AF_Nonlinear_mean = mean(AF_Nonlinear, na.rm = TRUE),
        AF_Nonlinear_low = quantile(AF_Nonlinear, 0.025, na.rm = TRUE),
        AF_Nonlinear_high = quantile(AF_Nonlinear, 0.975, na.rm = TRUE),
        death_change_Nonlinear_mean = mean(death_change_Nonlinear, na.rm = TRUE),
        death_change_Nonlinear_low = quantile(death_change_Nonlinear, 0.025, na.rm = TRUE),
        death_change_Nonlinear_high = quantile(death_change_Nonlinear, 0.975, na.rm = TRUE),
        annual_deaths_change_Nonlinear_mean = mean(annual_deaths_change_Nonlinear, na.rm = TRUE),
        annual_deaths_change_Nonlinear_low = quantile(annual_deaths_change_Nonlinear, 0.025, na.rm = TRUE),
        annual_deaths_change_Nonlinear_high = quantile(annual_deaths_change_Nonlinear, 0.975, na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  exposure_prediction_stats_all <- data.frame()
  gam_stats_all <- data.frame()
  threshold_table_all <- data.frame()
  policy_raw_indicator_snapshot_all <- data.frame()
  policy_score_snapshot_all <- data.frame()
  policy_exposure_prediction_snapshot_all <- data.frame()
  
  scenario_city_change_snapshot_all <- data.frame()
  policy_raw_indicator_city_snapshot_all <- data.frame()
  policy_score_city_snapshot_all <- data.frame()
  
  scenario_point_results_all <- data.frame()
  bootstrap_raw_all <- data.frame()
  reference_deaths_all <- data.frame()
  run_log_all <- data.frame()
  
  for (dname in names(data_list)) {
    if (isTRUE(show_message)) message("\nProcessing: ", dname, " [", time_scale, "]")
    
    common_df <- prepare_common_df(data_list[[dname]])
    ref <- compute_reference_deaths(common_df)
    
    reference_deaths_all <- bind_rows(reference_deaths_all, data.frame(
      disease = dname,
      time_scale = time_scale,
      base_year = base_year,
      death_reference_type = death_reference_type,
      base_year_observed_deaths_reference = ref$base_ref,
      avg_annual_observed_deaths_reference = ref$avg_ref,
      ReferenceDeaths = ref$ref_use,
      stringsAsFactors = FALSE
    ))
    
    vars_use <- intersect(exposure_vars, names(common_df))
    urban_use <- intersect(urban_exposure, names(common_df))
    
    if (length(vars_use) == 0 || length(urban_use) == 0) {
      run_log_all <- bind_rows(run_log_all, data.frame(
        disease = dname,
        time_scale = time_scale,
        message = "No exposure or urban variables found",
        stringsAsFactors = FALSE
      ))
      next
    }
    
    for (urban_var in urban_use) {
      prep <- prep_modifier_df(common_df, common_df[[urban_var]])
      df_u <- prep$data
      modifier_mean <- prep$modifier_mean
      modifier_sd <- prep$modifier_sd
      
      proj_future <- NULL
      if (!is.null(future_projection_table) && nrow(future_projection_table) > 0) {
        ft <- as.data.frame(future_projection_table)
        if (!(".city" %in% names(ft)) && "city" %in% names(ft)) ft$.city <- ft$city
        
        proj_future <- ft %>%
          filter(urban == urban_var, year %in% future_years)
      }
      
      for (exposure_var in vars_use) {
        if (isTRUE(show_message)) message("  Urban: ", urban_var, " | Exposure: ", exposure_var)
        
        exp_res <- tryCatch(fit_exposure_model(df_u, df_u, exposure_var), error = function(e) NULL)
        if (is.null(exp_res)) next
        
        df_u[[paste0(exposure_var, "_pred")]] <- exp_res$predicted[[paste0(exposure_var, "_pred")]]
        
        gam_list <- tryCatch(fit_mortality_models(df_u, exposure_var), error = function(e) NULL)
        if (is.null(gam_list)) next
        
        baseline_tot <- predict_totals(df_u, gam_list)
        
        exposure_prediction_stats_all <- bind_rows(exposure_prediction_stats_all, data.frame(
          disease = dname,
          time_scale = time_scale,
          analysis_type = "single_urban",
          urban = urban_var,
          exposure = exposure_var,
          model_type = exp_res$model_type,
          exposure_pred_r2 = exp_res$r2,
          exposure_pred_rmse = exp_res$rmse,
          stringsAsFactors = FALSE
        ))
        
        gam_stats_all <- bind_rows(gam_stats_all, data.frame(
          disease = dname,
          time_scale = time_scale,
          analysis_type = "single_urban",
          urban = urban_var,
          exposure = exposure_var,
          model = names(gam_list),
          dev_expl = vapply(gam_list, safe_dev_expl, numeric(1)),
          stringsAsFactors = FALSE
        ))
        
        if (isTRUE(run_future) &&
            (time_scale == "long" || isTRUE(run_future_short)) &&
            !is.null(proj_future) &&
            nrow(proj_future) > 0) {
          
          for (yy in future_years) {
            proj_y <- proj_future %>% filter(year == yy)
            if (nrow(proj_y) == 0) next
            
            base_same <- build_future_counterfactual(df_u, yy, NULL, modifier_mean, modifier_sd)
            base_same <- add_exposure_prediction(base_same, exp_res, exposure_var)
            base_tot <- predict_totals(base_same, gam_list)
            
            cf_full <- build_future_counterfactual(df_u, yy, proj_y, modifier_mean, modifier_sd)
            cf_full <- add_exposure_prediction(cf_full, exp_res, exposure_var)
            
            future_pathways <- make_pathway_counterfactuals(
              df_base = base_same,
              df_cf_full = cf_full,
              exposure_var = exposure_var
            )
            
            for (pw in names(future_pathways)) {
              cf_tot_pw <- predict_totals(future_pathways[[pw]], gam_list)
              
              scenario_point_results_all <- bind_rows(
                scenario_point_results_all,
                make_point_rows(
                  dname,
                  "future_projection",
                  urban_var,
                  exposure_var,
                  paste0("future_", yy),
                  base_tot,
                  cf_tot_pw,
                  extra = data.frame(
                    policy_framework = NA_character_,
                    delta = NA_real_,
                    pathway = pw,
                    future_year = yy,
                    stringsAsFactors = FALSE
                  )
                )
              )
              
              scenario_city_change_snapshot_all <- bind_rows(
                scenario_city_change_snapshot_all,
                build_scenario_city_change_snapshot(
                  df_base = base_same,
                  df_cf = future_pathways[[pw]],
                  gam_list = gam_list,
                  dname = dname,
                  analysis_type = "future_projection",
                  policy_framework = NA_character_,
                  scenario = paste0("future_", yy),
                  delta = NA_real_,
                  pathway = pw,
                  urban_var = urban_var,
                  exposure_var = exposure_var
                )
              )
            }
          }
        }
        
        thr_rule <- NULL
        
        if (isTRUE(run_threshold)) {
          curve_hit <- filter_curve_for_key(
            result_nonlinear_continuous_modification,
            dname,
            exposure_var,
            urban_var
          )
          
          if (nrow(curve_hit) > 0) {
            thr_rule <- classify_threshold_rule(curve_hit, df_u$modifier_raw)
            
            if (!is.null(thr_rule)) {
              threshold_table_all <- bind_rows(threshold_table_all, data.frame(
                disease = dname,
                time_scale = time_scale,
                urban = urban_var,
                exposure = exposure_var,
                thr_rule,
                stringsAsFactors = FALSE
              ))
              
              df_thr <- apply_threshold_to_modifier(df_u, thr_rule)
              df_thr <- refresh_modifier_scaled(df_thr, modifier_mean, modifier_sd)
              
              # full scenario 的 E'：用 U' 重新预测暴露
              df_thr_full <- add_exposure_prediction(df_thr, exp_res, exposure_var)
              
              # 拆成 full / modifier_only / exposure_only
              thr_pathways <- make_pathway_counterfactuals(
                df_base = df_u,
                df_cf_full = df_thr_full,
                exposure_var = exposure_var
              )
              
              for (pw in names(thr_pathways)) {
                thr_tot_pw <- predict_totals(thr_pathways[[pw]], gam_list)
                
                scenario_point_results_all <- bind_rows(
                  scenario_point_results_all,
                  make_point_rows(
                    dname,
                    "single_urban",
                    urban_var,
                    exposure_var,
                    paste0("threshold_", thr_rule$threshold_type[1]),
                    baseline_tot,
                    thr_tot_pw,
                    extra = bind_cols(
                      data.frame(
                        policy_framework = NA_character_,
                        delta = NA_real_,
                        pathway = pw,
                        stringsAsFactors = FALSE
                      ),
                      thr_rule
                    )
                  )
                )
                scenario_city_change_snapshot_all <- bind_rows(
                  scenario_city_change_snapshot_all,
                  build_scenario_city_change_snapshot(
                    df_base = df_u,
                    df_cf = thr_pathways[[pw]],
                    gam_list = gam_list,
                    dname = dname,
                    analysis_type = "single_urban",
                    policy_framework = NA_character_,
                    scenario = paste0("threshold_", thr_rule$threshold_type[1]),
                    delta = NA_real_,
                    pathway = pw,
                    urban_var = urban_var,
                    exposure_var = exposure_var
                  )
                )
              }
            }
          }
        }
        
        if (n_boot > 0) {
          boot_rows <- lapply(seq_len(n_boot), function(b) {
            set.seed(seed + b)
            df_b <- bootstrap_sample(df_u)
            
            if (nrow(df_b) < 10) return(NULL)
            
            m_mean_b <- mean(df_b$modifier_raw, na.rm = TRUE)
            m_sd_b <- safe_sd(df_b$modifier_raw)
            
            df_b <- refresh_modifier_scaled(df_b, m_mean_b, m_sd_b)
            
            exp_b <- tryCatch(
              fit_exposure_model(df_b, df_b, exposure_var),
              error = function(e) {
                message(
                  "[BOOT ERROR] exposure model failed | iter=", b,
                  " | disease=", dname,
                  " | urban=", urban_var,
                  " | exposure=", exposure_var,
                  " | error=", conditionMessage(e)
                )
                NULL
              }
            )
            
            if (is.null(exp_b)) {
              message(
                "[BOOT NULL] exposure model returned NULL | iter=", b,
                " | disease=", dname,
                " | urban=", urban_var,
                " | exposure=", exposure_var
              )
              return(NULL)
            }
            
            df_b[[paste0(exposure_var, "_pred")]] <- exp_b$predicted[[paste0(exposure_var, "_pred")]]
            
            gam_b <- tryCatch(
              fit_mortality_models(df_b, exposure_var),
              error = function(e) {
                message(
                  "[BOOT ERROR] mortality model failed | iter=", b,
                  " | disease=", dname,
                  " | urban=", urban_var,
                  " | exposure=", exposure_var,
                  " | error=", conditionMessage(e)
                )
                NULL
              }
            )
            
            if (is.null(gam_b)) {
              message(
                "[BOOT NULL] mortality model returned NULL | iter=", b,
                " | disease=", dname,
                " | urban=", urban_var,
                " | exposure=", exposure_var
              )
              return(NULL)
            }
            
            base_tot_b <- predict_totals(df_b, gam_b)
            out_b <- list()
            
            if (isTRUE(run_future) &&
                (time_scale == "long" || isTRUE(run_future_short)) &&
                !is.null(proj_future) &&
                nrow(proj_future) > 0) {
              
              for (yy in future_years) {
                proj_y <- proj_future %>% dplyr::filter(year == yy)
                if (nrow(proj_y) == 0) next
                
                base_same_b <- build_future_counterfactual(df_b, yy, NULL, m_mean_b, m_sd_b)
                base_same_b <- add_exposure_prediction(base_same_b, exp_b, exposure_var)
                base_tot_same_b <- predict_totals(base_same_b, gam_b)
                
                cf_b_full <- build_future_counterfactual(df_b, yy, proj_y, m_mean_b, m_sd_b)
                cf_b_full <- add_exposure_prediction(cf_b_full, exp_b, exposure_var)
                
                future_pathways_b <- make_pathway_counterfactuals(
                  df_base = base_same_b,
                  df_cf_full = cf_b_full,
                  exposure_var = exposure_var
                )
                
                for (pw in names(future_pathways_b)) {
                  cf_tot_b_pw <- predict_totals(future_pathways_b[[pw]], gam_b)
                  
                  out_b[[length(out_b) + 1]] <- data.frame(
                    disease = dname,
                    time_scale = time_scale,
                    analysis_type = "future_projection",
                    policy_framework = NA_character_,
                    pathway = pw,
                    urban = urban_var,
                    exposure = exposure_var,
                    scenario = paste0("future_", yy),
                    future_year = yy,
                    delta = NA_real_,
                    iter = b,
                    AF_Linear = if ("Linear" %in% names(cf_tot_b_pw)) {
                      (cf_tot_b_pw["Linear"] - base_tot_same_b["Linear"]) / base_tot_same_b["Linear"]
                    } else NA_real_,
                    AF_Nonlinear = if ("Nonlinear" %in% names(cf_tot_b_pw)) {
                      (cf_tot_b_pw["Nonlinear"] - base_tot_same_b["Nonlinear"]) / base_tot_same_b["Nonlinear"]
                    } else NA_real_,
                    death_change_Linear = if ("Linear" %in% names(cf_tot_b_pw)) {
                      cf_tot_b_pw["Linear"] - base_tot_same_b["Linear"]
                    } else NA_real_,
                    death_change_Nonlinear = if ("Nonlinear" %in% names(cf_tot_b_pw)) {
                      cf_tot_b_pw["Nonlinear"] - base_tot_same_b["Nonlinear"]
                    } else NA_real_,
                    stringsAsFactors = FALSE
                  )
                }
              }
            }
            
            if (isTRUE(run_threshold) && !is.null(thr_rule)) {
            
              df_thr_b <- apply_threshold_to_modifier(df_b, thr_rule)
              df_thr_b <- refresh_modifier_scaled(df_thr_b, m_mean_b, m_sd_b)
              
              df_thr_b_full <- add_exposure_prediction(df_thr_b, exp_b, exposure_var)
              
              thr_pathways_b <- make_pathway_counterfactuals(
                df_base = df_b,
                df_cf_full = df_thr_b_full,
                exposure_var = exposure_var
              )
              
              for (pw in names(thr_pathways_b)) {
                thr_tot_b_pw <- predict_totals(thr_pathways_b[[pw]], gam_b)
                
                out_b[[length(out_b) + 1]] <- data.frame(
                  disease = dname,
                  time_scale = time_scale,
                  analysis_type = "single_urban",
                  policy_framework = NA_character_,
                  pathway = pw,
                  urban = urban_var,
                  exposure = exposure_var,
                  scenario = paste0("threshold_", thr_rule$threshold_type[1]),
                  delta = NA_real_,
                  iter = b,
                  AF_Linear = if ("Linear" %in% names(thr_tot_b_pw)) {
                    (thr_tot_b_pw["Linear"] - base_tot_b["Linear"]) / base_tot_b["Linear"]
                  } else NA_real_,
                  AF_Nonlinear = if ("Nonlinear" %in% names(thr_tot_b_pw)) {
                    (thr_tot_b_pw["Nonlinear"] - base_tot_b["Nonlinear"]) / base_tot_b["Nonlinear"]
                  } else NA_real_,
                  death_change_Linear = if ("Linear" %in% names(thr_tot_b_pw)) {
                    thr_tot_b_pw["Linear"] - base_tot_b["Linear"]
                  } else NA_real_,
                  death_change_Nonlinear = if ("Nonlinear" %in% names(thr_tot_b_pw)) {
                    thr_tot_b_pw["Nonlinear"] - base_tot_b["Nonlinear"]
                  } else NA_real_,
                  stringsAsFactors = FALSE
                )
              }
            }
            
            if (length(out_b) == 0) return(NULL)
            bind_rows(out_b)
          })
          
          bootstrap_raw_all <- bind_rows(bootstrap_raw_all, bind_rows(boot_rows))
        }
      }
    }
    
    if (isTRUE(run_policy_integrated_be)) {
      if (is.null(pca_projection_params)) {
        stop("run_policy_integrated_be=TRUE 时必须提供 pca_projection_params。")
      }
      
      policy_registry <- get_policy_registry()
      delta <- as.numeric(policy_enhancement_delta[1])
      if (!is.finite(delta)) delta <- 0.30
      
      for (pp in seq_len(nrow(policy_registry))) {
        policy_framework <- policy_registry$policy_framework[pp]
        scenario_name <- paste0(policy_framework, "_raw_gap30")
        
        if (isTRUE(show_message)) message("  Policy scenario: ", scenario_name)
        
        cf_scores <- tryCatch(
          build_policy_counterfactual_raw(common_df, policy_framework, delta),
          error = function(e) NULL
        )
        
        if (is.null(cf_scores)) next
        
        policy_raw_indicator_snapshot_all <- bind_rows(
          policy_raw_indicator_snapshot_all,
          build_policy_raw_snapshot(common_df, cf_scores, dname, policy_framework, scenario_name, delta)
        )
        
        policy_score_snapshot_all <- bind_rows(
          policy_score_snapshot_all,
          build_policy_score_snapshot(common_df, cf_scores, dname, policy_framework, scenario_name, delta)
        )
        #
        policy_raw_indicator_city_snapshot_all <- bind_rows(
          policy_raw_indicator_city_snapshot_all,
          build_policy_raw_indicator_city_snapshot(
            df_base = common_df,
            df_cf = cf_scores,
            dname = dname,
            policy_framework = policy_framework,
            scenario = scenario_name,
            delta = delta
          )
        )
        
        policy_score_city_snapshot_all <- bind_rows(
          policy_score_city_snapshot_all,
          build_policy_score_city_snapshot(
            df_base = common_df,
            df_cf = cf_scores,
            dname = dname,
            policy_framework = policy_framework,
            scenario = scenario_name,
            delta = delta
          )
        )
        
        
        #
        for (urban_var in urban_use) {
          if (!(urban_var %in% names(cf_scores))) next
          
          prep <- prep_modifier_df(common_df, common_df[[urban_var]])
          df_u <- prep$data
          modifier_mean <- prep$modifier_mean
          modifier_sd <- prep$modifier_sd
          
          df_cf_base <- common_df
          df_cf_base$modifier_raw <- as.numeric(cf_scores[[urban_var]])
          df_cf_base <- refresh_modifier_scaled(df_cf_base, modifier_mean, modifier_sd)
          
          for (exposure_var in vars_use) {
            exp_res <- tryCatch(fit_exposure_model(df_u, df_u, exposure_var), error = function(e) NULL)
            if (is.null(exp_res)) next
            
            df_u[[paste0(exposure_var, "_pred")]] <- exp_res$predicted[[paste0(exposure_var, "_pred")]]
            
            gam_list <- tryCatch(fit_mortality_models(df_u, exposure_var), error = function(e) NULL)
            if (is.null(gam_list)) next
            
            baseline_tot <- predict_totals(df_u, gam_list)
            
            
            # full scenario 的 E'：用 U' 重新预测暴露
            df_cf_full <- add_exposure_prediction(df_cf_base, exp_res, exposure_var)
            
            policy_exposure_prediction_snapshot_all <- bind_rows(
              policy_exposure_prediction_snapshot_all,
              build_exposure_prediction_snapshot(
                df_base = df_u,
                df_cf = df_cf_full,
                dname = dname,
                analysis_type = "policy_integrated_be",
                policy_framework = policy_framework,
                scenario = scenario_name,
                delta = delta,
                urban_var = urban_var,
                exposure_var = exposure_var
              )
            )
            
            policy_pathways <- make_pathway_counterfactuals(
              df_base = df_u,
              df_cf_full = df_cf_full,
              exposure_var = exposure_var
            )
            
            for (pw in names(policy_pathways)) {
              cf_tot_pw <- predict_totals(policy_pathways[[pw]], gam_list)
              
              scenario_point_results_all <- bind_rows(
                scenario_point_results_all,
                make_point_rows(
                  dname,
                  "policy_integrated_be",
                  urban_var,
                  exposure_var,
                  scenario_name,
                  baseline_tot,
                  cf_tot_pw,
                  extra = data.frame(
                    policy_framework = policy_framework,
                    delta = delta,
                    pathway = pw,
                    stringsAsFactors = FALSE
                  )
                )
              )
              scenario_city_change_snapshot_all <- bind_rows(
                scenario_city_change_snapshot_all,
                build_scenario_city_change_snapshot(
                  df_base = df_u,
                  df_cf = policy_pathways[[pw]],
                  gam_list = gam_list,
                  dname = dname,
                  analysis_type = "policy_integrated_be",
                  policy_framework = policy_framework,
                  scenario = scenario_name,
                  delta = delta,
                  pathway = pw,
                  urban_var = urban_var,
                  exposure_var = exposure_var
                )
              )
            }
            
            exposure_prediction_stats_all <- bind_rows(exposure_prediction_stats_all, data.frame(
              disease = dname,
              time_scale = time_scale,
              analysis_type = "policy_integrated_be",
              policy_framework = policy_framework,
              scenario = scenario_name,
              delta = delta,
              urban = urban_var,
              exposure = exposure_var,
              model_type = exp_res$model_type,
              exposure_pred_r2 = exp_res$r2,
              exposure_pred_rmse = exp_res$rmse,
              stringsAsFactors = FALSE
            ))
            
            gam_stats_all <- bind_rows(gam_stats_all, data.frame(
              disease = dname,
              time_scale = time_scale,
              analysis_type = "policy_integrated_be",
              policy_framework = policy_framework,
              scenario = scenario_name,
              delta = delta,
              urban = urban_var,
              exposure = exposure_var,
              model = names(gam_list),
              dev_expl = vapply(gam_list, safe_dev_expl, numeric(1)),
              stringsAsFactors = FALSE
            ))
            
            
            if (n_boot > 0) {
              boot_rows <- lapply(seq_len(n_boot), function(b) {
                set.seed(seed + b)
                df_b0 <- bootstrap_sample(common_df)
                
                if (nrow(df_b0) < 10) return(NULL)
                
                cf_braw <- tryCatch(
                  build_policy_counterfactual_raw(df_b0, policy_framework, delta),
                  error = function(e) NULL
                )
                
                if (is.null(cf_braw) || !(urban_var %in% names(cf_braw))) return(NULL)
                
                prep_b <- prep_modifier_df(df_b0, df_b0[[urban_var]])
                df_b <- prep_b$data
                m_mean_b <- prep_b$modifier_mean
                m_sd_b <- prep_b$modifier_sd
                
                df_b_cf <- df_b0
                df_b_cf$modifier_raw <- as.numeric(cf_braw[[urban_var]])
                
                df_b_cf <- refresh_modifier_scaled(df_b_cf, m_mean_b, m_sd_b)
                
                exp_b <- tryCatch(fit_exposure_model(df_b, df_b, exposure_var), error = function(e) NULL)
                if (is.null(exp_b)) return(NULL)
                
                df_b[[paste0(exposure_var, "_pred")]] <- exp_b$predicted[[paste0(exposure_var, "_pred")]]
                
                gam_b <- tryCatch(fit_mortality_models(df_b, exposure_var), error = function(e) NULL)
                if (is.null(gam_b)) return(NULL)
                
                base_tot_b <- predict_totals(df_b, gam_b)
                
                df_b_cf_full <- add_exposure_prediction(df_b_cf, exp_b, exposure_var)
                
                policy_pathways_b <- make_pathway_counterfactuals(
                  df_base = df_b,
                  df_cf_full = df_b_cf_full,
                  exposure_var = exposure_var
                )
                
                out_policy_b <- lapply(names(policy_pathways_b), function(pw) {
                  cf_tot_b_pw <- predict_totals(policy_pathways_b[[pw]], gam_b)
                  
                  data.frame(
                    disease = dname,
                    time_scale = time_scale,
                    analysis_type = "policy_integrated_be",
                    policy_framework = policy_framework,
                    pathway = pw,
                    urban = urban_var,
                    exposure = exposure_var,
                    scenario = scenario_name,
                    delta = delta,
                    iter = b,
                    AF_Linear = if ("Linear" %in% names(cf_tot_b_pw)) {
                      (cf_tot_b_pw["Linear"] - base_tot_b["Linear"]) / base_tot_b["Linear"]
                    } else NA_real_,
                    AF_Nonlinear = if ("Nonlinear" %in% names(cf_tot_b_pw)) {
                      (cf_tot_b_pw["Nonlinear"] - base_tot_b["Nonlinear"]) / base_tot_b["Nonlinear"]
                    } else NA_real_,
                    death_change_Linear = if ("Linear" %in% names(cf_tot_b_pw)) {
                      cf_tot_b_pw["Linear"] - base_tot_b["Linear"]
                    } else NA_real_,
                    death_change_Nonlinear = if ("Nonlinear" %in% names(cf_tot_b_pw)) {
                      cf_tot_b_pw["Nonlinear"] - base_tot_b["Nonlinear"]
                    } else NA_real_,
                    stringsAsFactors = FALSE
                  )
                })
                
                bind_rows(out_policy_b)
              })
              
              bootstrap_raw_all <- bind_rows(bootstrap_raw_all, bind_rows(boot_rows))
            }
          }
        }
      }
    }
  }
  
  reference_key <- reference_deaths_all %>%
    dplyr::select(disease, death_reference_type, ReferenceDeaths) %>%
    distinct()
  
  if (nrow(scenario_point_results_all) > 0) {
    scenario_point_results_all <- scenario_point_results_all %>%
      left_join(reference_key, by = "disease") %>%
      mutate(
        annual_deaths_change = AF * ReferenceDeaths,
        annual_deaths_change_percent_based = AF_percent / 100 * ReferenceDeaths
      )
  }
  
  if (nrow(bootstrap_raw_all) > 0) {
    bootstrap_summary_all <- summarise_boot(bootstrap_raw_all, reference_key)
    
    bootstrap_raw_all <- bootstrap_raw_all %>%
      left_join(reference_key, by = "disease") %>%
      mutate(
        annual_deaths_change_Linear = AF_Linear * ReferenceDeaths,
        annual_deaths_change_Nonlinear = AF_Nonlinear * ReferenceDeaths
      )
  } else {
    bootstrap_summary_all <- data.frame()
  }
  
  write_if <- function(df, fn) {
    if (!is.null(df) && nrow(df) > 0) {
      write.csv(df, file.path(output_dir, fn), row.names = FALSE, fileEncoding = "UTF-8")
    }
  }
  
  write_if(reference_deaths_all, "reference_deaths_info.csv")
  write_if(exposure_prediction_stats_all, "exposure_prediction_stats_all.csv")
  write_if(gam_stats_all, "gam_stats_all.csv")
  write_if(threshold_table_all, "derived_threshold_table.csv")
  write_if(policy_raw_indicator_snapshot_all, "policy_raw_indicator_snapshot.csv")
  write_if(policy_score_snapshot_all, "policy_fixed_pca_score_snapshot.csv")
  write_if(policy_exposure_prediction_snapshot_all, "policy_exposure_prediction_snapshot.csv")
  
  write_if(scenario_city_change_snapshot_all, "scenario_city_change_snapshot.csv")
  write_if(policy_raw_indicator_city_snapshot_all, "policy_raw_indicator_city_snapshot.csv")
  write_if(policy_score_city_snapshot_all, "policy_score_city_snapshot.csv")
  
  write_if(scenario_point_results_all, "scenario_point_results_all.csv")
  write_if(bootstrap_raw_all, "bootstrap_raw_all.csv")
  write_if(bootstrap_summary_all, "bootstrap_summary_all.csv")
  write_if(run_log_all, "run_log.csv")
  
  invisible(list(
    reference_deaths_info = reference_deaths_all,
    exposure_prediction_stats = exposure_prediction_stats_all,
    gam_stats = gam_stats_all,
    threshold_table = threshold_table_all,
    policy_raw_indicator_snapshot = policy_raw_indicator_snapshot_all,
    policy_score_snapshot = policy_score_snapshot_all,
    policy_exposure_prediction_snapshot = policy_exposure_prediction_snapshot_all,
    scenario_city_change_snapshot = scenario_city_change_snapshot_all,
    policy_raw_indicator_city_snapshot = policy_raw_indicator_city_snapshot_all,
    policy_score_city_snapshot = policy_score_city_snapshot_all,
    scenario_point_results = scenario_point_results_all,
    bootstrap_raw = bootstrap_raw_all,
    bootstrap_summary = bootstrap_summary_all,
    run_log = run_log_all
  ))
}


# ============================================================
# File parser + parallel folder runner for scenario engine v4
# Replace the old run_af_from_rds_folder_v4() with this block
# ============================================================

parse_file_info_v4 <- function(fname, time_scale = c("long", "short", "auto")) {
  suppressPackageStartupMessages({
    library(stringr)
    library(tibble)
    library(dplyr)
  })
  
  time_scale <- match.arg(time_scale)
  
  nm0 <- tools::file_path_sans_ext(basename(fname))
  
  if (!stringr::str_starts(nm0, "diseases2_")) {
    return(tibble(
      disease = NA_character_,
      disease_base = NA_character_,
      disease_suffix = NA_character_,
      disease_version = NA_character_,
      strat_var = NA_character_,
      strat_level = NA_character_
    ))
  }
  
  rest <- stringr::str_remove(nm0, "^diseases2_")
  parts <- strsplit(rest, "_", fixed = TRUE)[[1]]
  
  if (length(parts) < 1) {
    return(tibble(
      disease = NA_character_,
      disease_base = NA_character_,
      disease_suffix = NA_character_,
      disease_version = NA_character_,
      strat_var = NA_character_,
      strat_level = NA_character_
    ))
  }
  
  is_short_file <- FALSE
  if (time_scale == "short") {
    is_short_file <- TRUE
  } else if (time_scale == "auto") {
    is_short_file <- length(parts) >= 2 && tolower(parts[2]) %in% c("short", "daily", "st")
  }
  
  # disease parsing
  if (is_short_file && length(parts) >= 2 && tolower(parts[2]) %in% c("short", "daily", "st")) {
    disease_part <- parts[1]
    disease_suffix <- parts[2]
    strata_parts <- if (length(parts) > 2) parts[3:length(parts)] else character(0)
    
    m <- stringr::str_match(disease_part, "^(total|cvd|resp|accident)(\\d*)$")
    base <- m[1, 2]
    ver <- m[1, 3]
    ver <- ifelse(is.na(ver) || ver == "", "", ver)
    
    disease_base <- if (!is.na(base)) paste0(base, ver) else disease_part
    disease <- paste0(disease_base, "_", disease_suffix)
    disease_version <- ifelse(ver == "", NA_character_, ver)
  } else {
    disease_part <- parts[1]
    disease_suffix <- NA_character_
    strata_parts <- if (length(parts) > 1) parts[2:length(parts)] else character(0)
    
    disease_base <- dplyr::case_when(
      grepl("^total2$", disease_part, ignore.case = TRUE) ~ "total2",
      grepl("^total$",  disease_part, ignore.case = TRUE) ~ "total2",
      grepl("^cvd2$", disease_part, ignore.case = TRUE) ~ "cvd2",
      grepl("^cvd$",  disease_part, ignore.case = TRUE) ~ "cvd2",
      grepl("^resp2$", disease_part, ignore.case = TRUE) ~ "resp2",
      grepl("^resp$",  disease_part, ignore.case = TRUE) ~ "resp2",
      grepl("^accident2$", disease_part, ignore.case = TRUE) ~ "accident2",
      grepl("^accident$",  disease_part, ignore.case = TRUE) ~ "accident2",
      TRUE ~ disease_part
    )
    
    disease <- disease_base
    disease_version <- stringr::str_match(disease_base, "([0-9]+)$")[1, 2]
  }
  
  # strata parsing
  strat_var <- NA_character_
  strat_level <- NA_character_
  
  if (length(strata_parts) > 0) {
    if (length(strata_parts) == 1) {
      s <- strata_parts[1]
      
      # sex2 / age1 / marital2 / education3 / urbanrate1
      m2 <- stringr::str_match(s, "^([A-Za-z]+)([0-9]+)$")
      
      if (!is.na(m2[1, 2])) {
        strat_var <- m2[1, 2]
        strat_level <- m2[1, 3]
      } else {
        strat_var <- s
        strat_level <- NA_character_
      }
    } else {
      # cluster_1 / Climate_Subtropical_Monsoon_Climate / GDP_low
      strat_var <- strata_parts[1]
      strat_level <- paste(strata_parts[-1], collapse = "_")
    }
    
    strat_var <- tolower(trimws(strat_var))
    strat_level <- trimws(strat_level)
    
    strat_var <- dplyr::case_when(
      strat_var %in% c("climate", "climatezone", "climate_zone") ~ "climate",
      strat_var %in% c("cluster", "becluster", "be_cluster") ~ "cluster",
      strat_var %in% c("urbanrate", "urbanization", "urbanizationrate", "urbanization_rate") ~ "urbanization",
      strat_var %in% c("gdp", "gdppc", "gdppercapita", "gdp_avg") ~ "gdp",
      strat_var %in% c("edu") ~ "education",
      TRUE ~ strat_var
    )
  }
  
  tibble(
    disease = disease,
    disease_base = disease_base,
    disease_suffix = disease_suffix,
    disease_version = disease_version,
    strat_var = strat_var,
    strat_level = strat_level
  )
}


make_file_info_v4 <- function(
    rds_dir,
    file_pattern = "^diseases2_.*\\.(rds|RDS)$",
    recursive = FALSE,
    time_scale = c("long", "short", "auto"),
    disease_keep = NULL,
    include_unstratified = TRUE,
    include_stratified = TRUE,
    exclude_strat_var = NULL,
    include_strat_var = NULL
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
    library(tibble)
  })
  
  time_scale <- match.arg(time_scale)
  
  files <- list.files(
    rds_dir,
    pattern = file_pattern,
    full.names = TRUE,
    recursive = recursive,
    ignore.case = TRUE
  )
  
  if (length(files) == 0) {
    stop("文件夹内没有找到 RDS 文件: ", rds_dir)
  }
  
  file_info <- purrr::map_dfr(files, function(f) {
    dplyr::bind_cols(
      tibble::tibble(
        file = f,
        source_file = basename(f),
        file_name = basename(f)
      ),
      parse_file_info_v4(
        fname = f,
        time_scale = time_scale
      )
    )
  }) %>%
    dplyr::filter(!is.na(disease_base)) %>%
    dplyr::mutate(
      is_stratified = !(is.na(strat_var) & is.na(strat_level))
    )
  
  # 疾病筛选
  if (!is.null(disease_keep)) {
    disease_keep <- trimws(as.character(disease_keep))
    
    file_info <- file_info %>%
      dplyr::filter(
        disease %in% disease_keep |
          disease_base %in% disease_keep
      )
  }
  
  # 是否保留未分层文件
  if (!isTRUE(include_unstratified)) {
    file_info <- file_info %>%
      dplyr::filter(is_stratified)
  }
  
  # 是否保留分层文件
  if (!isTRUE(include_stratified)) {
    file_info <- file_info %>%
      dplyr::filter(!is_stratified)
  }
  
  # 排除指定分层
  if (!is.null(exclude_strat_var)) {
    exclude_strat_var <- tolower(
      trimws(as.character(exclude_strat_var))
    )
    
    file_info <- file_info %>%
      dplyr::filter(
        is.na(strat_var) |
          !(strat_var %in% exclude_strat_var)
      )
  }
  
  # 仅保留指定分层
  # 注意：当 include_unstratified=TRUE 时，未分层文件仍然保留
  if (!is.null(include_strat_var)) {
    include_strat_var <- tolower(
      trimws(as.character(include_strat_var))
    )
    
    file_info <- file_info %>%
      dplyr::filter(
        is.na(strat_var) |
          strat_var %in% include_strat_var
      )
  }
  
  if (nrow(file_info) == 0) {
    stop(
      paste0(
        "筛选后没有可运行文件；请检查：",
        " disease_keep、include_unstratified、include_stratified、",
        "exclude_strat_var、include_strat_var 或文件命名。"
      )
    )
  }
  
  file_info
}


run_af_from_rds_folder_v4 <- function(
    rds_dir,
    time_scale = c("long", "short"),
    exposure_vars,
    urban_exposure,
    result_nonlinear_continuous_modification = NULL,
    future_projection_table = NULL,
    pca_projection_params = NULL,
    
    file_pattern = "^diseases2_.*\\.(rds|RDS)$",
    recursive = FALSE,
    
    disease_keep = NULL,
    include_unstratified = TRUE,
    include_stratified = TRUE,
    exclude_strat_var = NULL,
    include_strat_var = NULL,
    
    output_dir,
    n_workers = 1,
    ...
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
    library(tibble)
  })
  
  time_scale <- match.arg(time_scale)
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  file_info <- make_file_info_v4(
    rds_dir = rds_dir,
    file_pattern = file_pattern,
    recursive = recursive,
    time_scale = time_scale,
    disease_keep = disease_keep,
    include_unstratified = include_unstratified,
    include_stratified = include_stratified,
    exclude_strat_var = exclude_strat_var,
    include_strat_var = include_strat_var
  )
  
  write.csv(
    file_info,
    file.path(output_dir, "file_info_parsed.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  
  message("\nParsed files:")
  print(file_info %>% dplyr::select(source_file, disease, disease_base, disease_suffix, strat_var, strat_level, is_stratified))
  
  run_one_file <- function(i) {
    info_i <- file_info[i, ]
    
    message("\n========================================")
    message("Running file: ", info_i$source_file)
    message(
      "disease = ", info_i$disease,
      " | disease_base = ", info_i$disease_base,
      " | strat_var = ", info_i$strat_var,
      " | strat_level = ", info_i$strat_level
    )
    message("========================================")
    
    dat_i <- tryCatch(
      readRDS(info_i$file),
      error = function(e) {
        warning("读取失败: ", info_i$file, " | ", conditionMessage(e))
        NULL
      }
    )
    
    if (is.null(dat_i)) {
      return(list(
        reference_deaths_info = data.frame(),
        exposure_prediction_stats = data.frame(),
        gam_stats = data.frame(),
        threshold_table = data.frame(),
        policy_raw_indicator_snapshot = data.frame(),
        policy_score_snapshot = data.frame(),
        policy_exposure_prediction_snapshot = data.frame(),
        scenario_city_change_snapshot = data.frame(),
        policy_raw_indicator_city_snapshot = data.frame(),
        policy_score_city_snapshot = data.frame(),
        scenario_point_results = data.frame(),
        bootstrap_raw = data.frame(),
        bootstrap_summary = data.frame(),
        run_log = tibble::tibble(
          source_file = info_i$source_file,
          disease = info_i$disease,
          disease_base = info_i$disease_base,
          disease_suffix = info_i$disease_suffix,
          disease_version = info_i$disease_version,
          strat_var = info_i$strat_var,
          strat_level = info_i$strat_level,
          is_stratified_file = info_i$is_stratified,
          status = "read_failed"
        )
      ))
    }
    
    # 与 named-list runner 保持一致；
    # 使用标准疾病名称匹配 nonlinear result 中的 total2、cvd2 等
    dname_run <- info_i$disease_base
    
    out_i <- tryCatch(
      run_urban_resilience_scenarios_final_v4(
        data_list = setNames(list(dat_i), dname_run),
        time_scale = time_scale,
        exposure_vars = exposure_vars,
        urban_exposure = urban_exposure,
        result_nonlinear_continuous_modification = result_nonlinear_continuous_modification,
        future_projection_table = future_projection_table,
        pca_projection_params = pca_projection_params,
        output_dir = file.path(output_dir, tools::file_path_sans_ext(info_i$source_file)),
        ...
      ),
      error = function(e) {
        warning("运行失败: ", info_i$source_file, " | ", conditionMessage(e))
        NULL
      }
    )
    
    if (is.null(out_i)) {
      return(list(
        reference_deaths_info = data.frame(),
        exposure_prediction_stats = data.frame(),
        gam_stats = data.frame(),
        threshold_table = data.frame(),
        policy_raw_indicator_snapshot = data.frame(),
        policy_score_snapshot = data.frame(),
        policy_exposure_prediction_snapshot = data.frame(),
        
        scenario_city_change_snapshot = data.frame(),
        policy_raw_indicator_city_snapshot = data.frame(),
        policy_score_city_snapshot = data.frame(),
        
        scenario_point_results = data.frame(),
        bootstrap_raw = data.frame(),
        bootstrap_summary = data.frame(),
        
        run_log = tibble::tibble(
          source_file = info_i$source_file,
          disease = info_i$disease,
          disease_base = info_i$disease_base,
          disease_suffix = info_i$disease_suffix,
          disease_version = info_i$disease_version,
          strat_var = info_i$strat_var,
          strat_level = info_i$strat_level,
          is_stratified_file = info_i$is_stratified,
          status = "run_failed"
        )
      ))
    }
    
    tag_df <- function(df) {
      if (is.null(df) || nrow(df) == 0) return(data.frame())
      
      dplyr::bind_cols(
        tibble::tibble(
          source_file = rep(info_i$source_file, nrow(df)),
          disease_file = rep(info_i$disease, nrow(df)),
          disease_base = rep(info_i$disease_base, nrow(df)),
          disease_suffix = rep(info_i$disease_suffix, nrow(df)),
          disease_version = rep(info_i$disease_version, nrow(df)),
          strat_var = rep(info_i$strat_var, nrow(df)),
          strat_level = rep(info_i$strat_level, nrow(df)),
          is_stratified_file = rep(info_i$is_stratified, nrow(df))
        ),
        df
      )
    }
    
    out_i <- lapply(out_i, tag_df)
    
    out_i$run_log <- tibble::tibble(
      source_file = info_i$source_file,
      disease = info_i$disease,
      disease_base = info_i$disease_base,
      disease_suffix = info_i$disease_suffix,
      disease_version = info_i$disease_version,
      strat_var = info_i$strat_var,
      strat_level = info_i$strat_level,
      is_stratified_file = info_i$is_stratified,
      status = "success"
    )
    
    out_i
  }
  
  if (isTRUE(n_workers > 1)) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("请先安装 future.apply: install.packages('future.apply')")
    }
    if (!requireNamespace("future", quietly = TRUE)) {
      stop("请先安装 future: install.packages('future')")
    }
    
    old_plan <- future::plan()
    on.exit({
      try(future::plan(old_plan), silent = TRUE)
    }, add = TRUE)
    
    future::plan(future::multisession, workers = n_workers)
    
    all_out <- future.apply::future_lapply(
      seq_len(nrow(file_info)),
      run_one_file,
      future.seed = TRUE
    )
  } else {
    all_out <- lapply(seq_len(nrow(file_info)), run_one_file)
  }
  
  bind_component <- function(nm) {
    dplyr::bind_rows(lapply(all_out, function(x) {
      if (is.null(x[[nm]])) data.frame() else x[[nm]]
    }))
  }
  
  merged <- list(
    file_info = file_info,
    reference_deaths_info = bind_component("reference_deaths_info"),
    exposure_prediction_stats = bind_component("exposure_prediction_stats"),
    gam_stats = bind_component("gam_stats"),
    threshold_table = bind_component("threshold_table"),
    policy_raw_indicator_snapshot = bind_component("policy_raw_indicator_snapshot"),
    policy_score_snapshot = bind_component("policy_score_snapshot"),
    policy_exposure_prediction_snapshot = bind_component("policy_exposure_prediction_snapshot"),
    scenario_city_change_snapshot = bind_component("scenario_city_change_snapshot"),
    policy_raw_indicator_city_snapshot = bind_component("policy_raw_indicator_city_snapshot"),
    policy_score_city_snapshot = bind_component("policy_score_city_snapshot"),
    scenario_point_results = bind_component("scenario_point_results"),
    bootstrap_raw = bind_component("bootstrap_raw"),
    bootstrap_summary = bind_component("bootstrap_summary"),
    run_log = bind_component("run_log")
  )
  
  write_if <- function(df, fn) {
    if (!is.null(df) && nrow(df) > 0) {
      write.csv(
        df,
        file.path(output_dir, fn),
        row.names = FALSE,
        fileEncoding = "UTF-8"
      )
    }
  }
  
  write_if(merged$file_info, "ALL_file_info_parsed.csv")
  write_if(merged$reference_deaths_info, "ALL_reference_deaths_info.csv")
  write_if(merged$exposure_prediction_stats, "ALL_exposure_prediction_stats.csv")
  write_if(merged$gam_stats, "ALL_gam_stats.csv")
  write_if(merged$threshold_table, "ALL_threshold_table.csv")
  write_if(merged$policy_raw_indicator_snapshot, "ALL_policy_raw_indicator_snapshot.csv")
  write_if(merged$policy_score_snapshot, "ALL_policy_fixed_pca_score_snapshot.csv")
  write_if(merged$policy_exposure_prediction_snapshot, "ALL_policy_exposure_prediction_snapshot.csv")
  write_if(merged$scenario_city_change_snapshot, "ALL_scenario_city_change_snapshot.csv")
  write_if(merged$policy_raw_indicator_city_snapshot, "ALL_policy_raw_indicator_city_snapshot.csv")
  write_if(merged$policy_score_city_snapshot, "ALL_policy_score_city_snapshot.csv")
   write_if(merged$scenario_point_results, "ALL_scenario_point_results.csv")
  write_if(merged$bootstrap_raw, "ALL_bootstrap_raw.csv")
  write_if(merged$bootstrap_summary, "ALL_bootstrap_summary.csv")
  write_if(merged$run_log, "ALL_run_log.csv")
  
  invisible(merged)
}

# ============================================================
# Named-list runner for result_list_merged.rds
# 用于 result_list_merged[["result_total2_Climate_xxx"]] 这种结构
# ============================================================

parse_result_list_name_v4 <- function(x) {
  suppressPackageStartupMessages({
    library(stringr)
    library(tibble)
    library(dplyr)
  })
  
  nm0 <- as.character(x)
  nm <- nm0
  
  # 去掉 result_ 前缀
  nm <- stringr::str_remove(nm, "^result_")
  
  parts <- strsplit(nm, "_", fixed = TRUE)[[1]]
  
  disease_part <- parts[1]
  strata_parts <- if (length(parts) > 1) parts[2:length(parts)] else character(0)
  
  disease_base <- dplyr::case_when(
    grepl("^total2$", disease_part, ignore.case = TRUE) ~ "total2",
    grepl("^total$",  disease_part, ignore.case = TRUE) ~ "total2",
    grepl("^cvd2$", disease_part, ignore.case = TRUE) ~ "cvd2",
    grepl("^cvd$",  disease_part, ignore.case = TRUE) ~ "cvd2",
    grepl("^resp2$", disease_part, ignore.case = TRUE) ~ "resp2",
    grepl("^resp$",  disease_part, ignore.case = TRUE) ~ "resp2",
    grepl("^accident2$", disease_part, ignore.case = TRUE) ~ "accident2",
    grepl("^accident$",  disease_part, ignore.case = TRUE) ~ "accident2",
    TRUE ~ disease_part
  )
  
  disease_version <- stringr::str_match(disease_base, "([0-9]+)$")[1, 2]
  
  strat_var <- NA_character_
  strat_level <- NA_character_
  
  if (length(strata_parts) > 0) {
    if (length(strata_parts) == 1) {
      s <- strata_parts[1]
      
      # sex2 / age1 / marital2 / education3 / urbanrate1
      m2 <- stringr::str_match(s, "^([A-Za-z]+)([0-9]+)$")
      
      if (!is.na(m2[1, 2])) {
        strat_var <- m2[1, 2]
        strat_level <- m2[1, 3]
      } else {
        strat_var <- s
        strat_level <- NA_character_
      }
    } else {
      # Climate_Alpine Plateau Climate / Cluster_1 / GDP_low
      strat_var <- strata_parts[1]
      strat_level <- paste(strata_parts[-1], collapse = "_")
    }
    
    strat_var <- tolower(trimws(strat_var))
    strat_level <- trimws(strat_level)
    
    strat_var <- dplyr::case_when(
      strat_var %in% c("climate", "climatezone", "climate_zone") ~ "climate",
      strat_var %in% c("cluster", "becluster", "be_cluster") ~ "cluster",
      strat_var %in% c("urbanrate", "urbanization", "urbanizationrate", "urbanization_rate") ~ "urbanization",
      strat_var %in% c("gdp", "gdppc", "gdppercapita", "gdp_avg") ~ "gdp",
      strat_var %in% c("edu") ~ "education",
      TRUE ~ strat_var
    )
  }
  
  tibble::tibble(
    list_name = nm0,
    disease = disease_base,
    disease_base = disease_base,
    disease_version = disease_version,
    strat_var = strat_var,
    strat_level = strat_level,
    is_stratified = !(is.na(strat_var) & is.na(strat_level))
  )
}


make_list_info_v4 <- function(
    result_list,
    disease_keep = "total2",
    include_unstratified = TRUE,
    include_stratified = TRUE,
    exclude_strat_var = "climate",
    include_strat_var = NULL
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
  })
  
  info <- purrr::map_dfr(names(result_list), parse_result_list_name_v4)
  
  if (!is.null(disease_keep)) {
    info <- info %>%
      dplyr::filter(disease_base %in% disease_keep | disease %in% disease_keep)
  }
  
  if (!isTRUE(include_unstratified)) {
    info <- info %>% dplyr::filter(is_stratified)
  }
  
  if (!isTRUE(include_stratified)) {
    info <- info %>% dplyr::filter(!is_stratified)
  }
  
  if (!is.null(exclude_strat_var)) {
    info <- info %>%
      dplyr::filter(is.na(strat_var) | !(strat_var %in% exclude_strat_var))
  }
  
  if (!is.null(include_strat_var)) {
    info <- info %>%
      dplyr::filter(is.na(strat_var) | strat_var %in% include_strat_var)
  }
  
  info
}


run_af_from_named_list_v4 <- function(
    result_list,
    time_scale = c("long", "short"),
    exposure_vars,
    urban_exposure,
    result_nonlinear_continuous_modification = NULL,
    future_projection_table = NULL,
    pca_projection_params = NULL,
    
    disease_keep = "total2",
    include_unstratified = TRUE,
    include_stratified = TRUE,
    exclude_strat_var = "climate",
    include_strat_var = NULL,
    
    output_dir,
    n_workers = 1,
    ...
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
    library(tibble)
  })
  
  time_scale <- match.arg(time_scale)
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  list_info <- make_list_info_v4(
    result_list = result_list,
    disease_keep = disease_keep,
    include_unstratified = include_unstratified,
    include_stratified = include_stratified,
    exclude_strat_var = exclude_strat_var,
    include_strat_var = include_strat_var
  )
  
  if (nrow(list_info) == 0) {
    stop("筛选后没有可运行对象；检查 disease_keep / include_stratified / exclude_strat_var。")
  }
  
  write.csv(
    list_info,
    file.path(output_dir, "ALL_list_info_parsed.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  
  message("\n准备运行以下 list 元素：")
  print(list_info %>%
          dplyr::select(list_name, disease, disease_base, strat_var, strat_level, is_stratified))
  
  run_one_item <- function(i) {
    info_i <- list_info[i, ]
    
    message("\n========================================")
    message("Running list item: ", info_i$list_name)
    message(
      "disease = ", info_i$disease,
      " | strat_var = ", info_i$strat_var,
      " | strat_level = ", info_i$strat_level
    )
    message("========================================")
    
    dat_i <- result_list[[info_i$list_name]]
    
    if (is.null(dat_i)) {
      return(list(
        reference_deaths_info = data.frame(),
        exposure_prediction_stats = data.frame(),
        gam_stats = data.frame(),
        threshold_table = data.frame(),
        policy_raw_indicator_snapshot = data.frame(),
        policy_score_snapshot = data.frame(),
        policy_exposure_prediction_snapshot = data.frame(),
        scenario_point_results = data.frame(),
        bootstrap_raw = data.frame(),
        bootstrap_summary = data.frame(),
        run_log = tibble::tibble(
          source_file = info_i$list_name,
          disease = info_i$disease,
          disease_base = info_i$disease_base,
          strat_var = info_i$strat_var,
          strat_level = info_i$strat_level,
          is_stratified_file = info_i$is_stratified,
          status = "object_null"
        )
      ))
    }
    
    # 注意：这里用 disease_base 作为 run_urban 的 disease 名称
    # 这样可以匹配 nonlinear csv 里的 total2 / accident2 等 disease 字段
    dname_run <- info_i$disease_base
    
    out_i <- tryCatch(
      run_urban_resilience_scenarios_final_v4(
        data_list = setNames(list(dat_i), dname_run),
        time_scale = time_scale,
        exposure_vars = exposure_vars,
        urban_exposure = urban_exposure,
        result_nonlinear_continuous_modification = result_nonlinear_continuous_modification,
        future_projection_table = future_projection_table,
        pca_projection_params = pca_projection_params,
        output_dir = file.path(output_dir, make.names(info_i$list_name)),
        ...
      ),
      error = function(e) {
        warning("运行失败: ", info_i$list_name, " | ", conditionMessage(e))
        NULL
      }
    )
    
    if (is.null(out_i)) {
      return(list(
        reference_deaths_info = data.frame(),
        exposure_prediction_stats = data.frame(),
        gam_stats = data.frame(),
        threshold_table = data.frame(),
        policy_raw_indicator_snapshot = data.frame(),
        policy_score_snapshot = data.frame(),
        policy_exposure_prediction_snapshot = data.frame(),
        scenario_point_results = data.frame(),
        bootstrap_raw = data.frame(),
        bootstrap_summary = data.frame(),
        run_log = tibble::tibble(
          source_file = info_i$list_name,
          disease = info_i$disease,
          disease_base = info_i$disease_base,
          strat_var = info_i$strat_var,
          strat_level = info_i$strat_level,
          is_stratified_file = info_i$is_stratified,
          status = "run_failed"
        )
      ))
    }
    
    tag_df <- function(df) {
      if (is.null(df) || nrow(df) == 0) return(data.frame())
      
      dplyr::bind_cols(
        tibble::tibble(
          source_file = rep(info_i$list_name, nrow(df)),
          disease_file = rep(info_i$disease, nrow(df)),
          disease_base = rep(info_i$disease_base, nrow(df)),
          disease_version = rep(info_i$disease_version, nrow(df)),
          strat_var = rep(info_i$strat_var, nrow(df)),
          strat_level = rep(info_i$strat_level, nrow(df)),
          is_stratified_file = rep(info_i$is_stratified, nrow(df))
        ),
        df
      )
    }
    
    out_i <- lapply(out_i, tag_df)
    
    out_i$run_log <- tibble::tibble(
      source_file = info_i$list_name,
      disease = info_i$disease,
      disease_base = info_i$disease_base,
      disease_version = info_i$disease_version,
      strat_var = info_i$strat_var,
      strat_level = info_i$strat_level,
      is_stratified_file = info_i$is_stratified,
      status = "success"
    )
    
    out_i
  }
  
  if (isTRUE(n_workers > 1)) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("请先安装 future.apply: install.packages('future.apply')")
    }
    if (!requireNamespace("future", quietly = TRUE)) {
      stop("请先安装 future: install.packages('future')")
    }
    
    old_plan <- future::plan()
    on.exit({
      try(future::plan(old_plan), silent = TRUE)
    }, add = TRUE)
    
    future::plan(future::multisession, workers = n_workers)
    
    all_out <- future.apply::future_lapply(
      seq_len(nrow(list_info)),
      run_one_item,
      future.seed = TRUE
    )
  } else {
    all_out <- lapply(seq_len(nrow(list_info)), run_one_item)
  }
  
  bind_component <- function(nm) {
    dplyr::bind_rows(lapply(all_out, function(x) {
      if (is.null(x[[nm]])) data.frame() else x[[nm]]
    }))
  }
  
  merged <- list(
    list_info = list_info,
    reference_deaths_info = bind_component("reference_deaths_info"),
    exposure_prediction_stats = bind_component("exposure_prediction_stats"),
    gam_stats = bind_component("gam_stats"),
    threshold_table = bind_component("threshold_table"),
    policy_raw_indicator_snapshot = bind_component("policy_raw_indicator_snapshot"),
    policy_score_snapshot = bind_component("policy_score_snapshot"),
    policy_exposure_prediction_snapshot = bind_component("policy_exposure_prediction_snapshot"),
    scenario_city_change_snapshot = bind_component("scenario_city_change_snapshot"),
    policy_raw_indicator_city_snapshot = bind_component("policy_raw_indicator_city_snapshot"),
    policy_score_city_snapshot = bind_component("policy_score_city_snapshot"),
    scenario_point_results = bind_component("scenario_point_results"),
    bootstrap_raw = bind_component("bootstrap_raw"),
    bootstrap_summary = bind_component("bootstrap_summary"),
    run_log = bind_component("run_log")
  )
  
  write_if <- function(df, fn) {
    if (!is.null(df) && nrow(df) > 0) {
      write.csv(
        df,
        file.path(output_dir, fn),
        row.names = FALSE,
        fileEncoding = "UTF-8"
      )
    }
  }
  
  write_if(merged$list_info, "ALL_list_info_parsed.csv")
  write_if(merged$reference_deaths_info, "ALL_reference_deaths_info.csv")
  write_if(merged$exposure_prediction_stats, "ALL_exposure_prediction_stats.csv")
  write_if(merged$gam_stats, "ALL_gam_stats.csv")
  write_if(merged$threshold_table, "ALL_threshold_table.csv")
  write_if(merged$policy_raw_indicator_snapshot, "ALL_policy_raw_indicator_snapshot.csv")
  write_if(merged$policy_score_snapshot, "ALL_policy_fixed_pca_score_snapshot.csv")
  write_if(merged$scenario_point_results, "ALL_scenario_point_results.csv")
  write_if(merged$policy_exposure_prediction_snapshot, "ALL_policy_exposure_prediction_snapshot.csv")
  write_if(merged$scenario_city_change_snapshot, "ALL_scenario_city_change_snapshot.csv")
  write_if(merged$policy_raw_indicator_city_snapshot, "ALL_policy_raw_indicator_city_snapshot.csv")
  write_if(merged$policy_score_city_snapshot, "ALL_policy_score_city_snapshot.csv")
  
  
  write_if(merged$bootstrap_raw, "ALL_bootstrap_raw.csv")
  write_if(merged$bootstrap_summary, "ALL_bootstrap_summary.csv")
  write_if(merged$run_log, "ALL_run_log.csv")
  
  
  invisible(merged)
}
library(tibble)
library(stringr)
library(dplyr)
library(tidyverse)
library(data.table)
library(Hmisc)
library(lubridate)
library(stringr)
library(sf)
library(purrr)
library(splines)
library(dplyr)
library(mgcv)
library(dlnm)
library(ggplot2)
library(interactionR)
library(metafor)
library(ggh4x)
library(ggsci)
library(dlnm)
library(splines)
library(xtable)
library(mvmeta)
library(ggplot2)
library(patchwork)
library(compareGroups)
library(terra)
#install.packages("magrittr")
#install.packages("exactextractr")
#library(exactextractr)
#install.packages("furrr")
library(furrr)
#install.packages("humidity")
library(sp)
library("humidity")
library(magrittr)
library(haven)
library(doParallel)
library(parallel)
library(geepack)
library(zoo)
library(gnm)
library(MASS)
library(future.apply)
library(haven)
library(tableone)
library(cowplot)
library(boot)
library(future)
#install.packages('doRNG')
library(doRNG)
# library(geosphere)
# data<-within(data,{
#   wgs84_coords<-bd2wgs(bd_lon,bd_lat)
#   wgs84_lon<-wgs84_coords[,1]
#   wgs84_lat<-wgs84_coords[,2]
# })

##############################################################读取下面的数据##############
########################################################################################
#######################################################################################
#####################################################################################
#读取数据####
rm(list=ls())
diseases2_total2<-readRDS("E:/PKU_fanyi/PM2.5/data_city2/diseases2_total2.rds")
diseases2_cvd2<-readRDS("E:/PKU_fanyi/PM2.5/data_city2/diseases2_cvd2.rds")
diseases2_resp2<-readRDS("E:/PKU_fanyi/PM2.5/data_city2/diseases2_resp2.rds")
diseases2_accident2<-readRDS("E:/PKU_fanyi/PM2.5/data_city2/diseases2_accident2.rds")
#
diseases2_total2<-readRDS("E:/TMP/data/diseases2_total2.rds")
diseases2_cvd2<-readRDS("E:/TMP/data/diseases2_cvd2.rds")
diseases2_resp2<-readRDS("E:/TMP/data/diseases2_resp2.rds")
diseases2_accident2<-readRDS("E:/TMP/data/diseases2_accident2.rds")

#做城市因素与结局直接关联#####

# 1) 暴露变量清单 & 数据集列表
datasets <- list(
  total    = diseases2_total2,
  cvd      = diseases2_cvd2,
  resp     = diseases2_resp2,
  accident = diseases2_accident2
)

ns_value <- 4


# 2) 5维度映射
#
urban_environment_vars <- c(
  #"urbanization_rate","popolation_age_65_Ratio",
  #"Male_Ratio",
  #"Cropland",
  "forest_shrub_grassland","water_wetland","impervious",
  paste0(c(
    #"total.population",
    "Residential","Workplace","Industrial_production","Education_culture",
    "Healthcare","Government","Commerce_services","Leisure_tourism","Logistics_delivery",
    "Urban_public_transit","Intercity_passenger_hub","Motorized_transport_support",
    "total_len_km","major_len_km","minor_len_km"
  ), "_per_totalArea")
)

# PDS_vars <- c("urbanization_rate","popolation_age_65_Ratio",
#               #"Male_Ratio",
#               "total.population_per_totalArea")
LandCover_vars <- c(
  #"Cropland",
  "forest_shrub_grassland","water_wetland","impervious")
Functional_Service_Density_vars <- c(
  "Residential_per_totalArea","Workplace_per_totalArea","Industrial_production_per_totalArea",
  "Education_culture_per_totalArea","Healthcare_per_totalArea","Government_per_totalArea",
  "Commerce_services_per_totalArea","Leisure_tourism_per_totalArea","Logistics_delivery_per_totalArea"
)
Transport_Facility_Density_vars <- c(
  "Urban_public_transit_per_totalArea","Intercity_passenger_hub_per_totalArea","Motorized_transport_support_per_totalArea"
)
Road_Network_Density_vars <- c(
  "total_len_km_per_totalArea","major_len_km_per_totalArea","minor_len_km_per_totalArea"
)
domain_levels <- c(
                   #"PDS",
                   "LandCover","Functional_Service_Density",
                   "Transport_Facility_Density","Road_Network_Density")
# --- 1) 原始指标 domain_map（你已有逻辑）---
domain_map_raw <- tibble(exposure = urban_environment_vars) %>%
  mutate(
    domain = case_when(
     # exposure %in% PDS_vars ~ "PDS",
      exposure %in% LandCover_vars ~ "LandCover",
      exposure %in% Functional_Service_Density_vars ~ "Functional_Service_Density",
      exposure %in% Transport_Facility_Density_vars ~ "Transport_Facility_Density",
      exposure %in% Road_Network_Density_vars ~ "Road_Network_Density",
      TRUE ~ NA_character_
    ),
    feature_type = "raw",
    score_scope = NA_character_
  ) %>%
  filter(!is.na(domain))

# --- 2) PCA score 列（维度得分 + 总得分）---
# 维度得分列名：score_<domain>_pca
pca_domain_scores <- tibble(
  exposure = paste0("score_", domain_levels, "_pca"),
  domain = domain_levels,
  feature_type = "score_pca",
  score_scope = "domain"
)

# 总分列名（按你现在生成的列名）
pca_total_scores <- tibble(
  exposure = c("score_total_pca_allvars", "score_total_pca_2stage"),
  domain = "Total",
  feature_type = "score_pca",
  score_scope = "total"
)

domain_map_pca <- bind_rows(pca_domain_scores, pca_total_scores)

# --- 3) FA score 列（维度得分 + 总得分）---
fa_domain_scores <- tibble(
  exposure = paste0("score_", domain_levels, "_fa"),
  domain = domain_levels,
  feature_type = "score_fa",
  score_scope = "domain"
)

fa_total_scores <- tibble(
  exposure = c("score_total_fa_allvars", "score_total_fa_2stage"),
  domain = "Total",
  feature_type = "score_fa",
  score_scope = "total"
)

domain_map_fa <- bind_rows(fa_domain_scores, fa_total_scores)

# --- 4) 合并成一个总 domain_map ---
domain_map <- bind_rows(domain_map_raw, domain_map_pca, domain_map_fa)
domain_map<- domain_map %>%
  mutate(
    # 1) domain顺序
    domain_order = match(domain, c(domain_levels, "Total")),
    
    # 2) 类型顺序：raw -> score_pca -> score_fa
    type_order = case_when(
      feature_type == "raw" ~ 1L,
      feature_type == "score_pca" ~ 2L,
      feature_type == "score_fa" ~ 3L,
      TRUE ~ 99L
    ),
    
    # 3) score内部顺序：domain 得分在前，total 在后；raw 用 1 占位
    scope_order = case_when(
      feature_type == "raw" ~ 1L,
      score_scope == "domain" ~ 1L,
      score_scope == "total" ~ 2L,
      TRUE ~ 99L
    ),
    
    # 4) raw 内部顺序：严格按你最开始 urban_environment_vars 的顺序
    raw_order = if_else(feature_type == "raw",
                        match(exposure, urban_environment_vars),
                        NA_integer_),
    
    # 5) score 内部：按 domain_levels 排序；Total 放最后；total里再规定 allvars 在 2stage 前
    score_order = case_when(
      feature_type != "raw" & domain %in% domain_levels ~ match(domain, domain_levels),
      feature_type != "raw" & domain == "Total" ~ 999L,
      TRUE ~ NA_integer_
    ),
    total_score_order = case_when(
      exposure %in% c("score_total_pca_allvars","score_total_fa_allvars") ~ 1L,
      exposure %in% c("score_total_pca_2stage","score_total_fa_2stage") ~ 2L,
      TRUE ~ NA_integer_
    )
  ) %>%
  arrange(domain_order, type_order, scope_order,
          raw_order, score_order, total_score_order, exposure) %>%
  dplyr::select(exposure, domain, feature_type, score_scope)
urban_environment_vars<-domain_map$exposure

# 3) 工具函数：给定 beta,se 与 delta，输出 RR/PC + 95%CI（正确缩放）
effect_with_delta <- function(beta, se, delta, z = qnorm(0.975)) {
  # logRR = beta * delta
  logrr  <- beta * delta
  selog  <- se   * delta
  rr     <- exp(logrr)
  rr_l   <- exp(logrr - z * selog)
  rr_u   <- exp(logrr + z * selog)
  
  pc   <- (rr   - 1) * 100
  pc_l <- (rr_l - 1) * 100
  pc_u <- (rr_u - 1) * 100
  
  list(rr = rr, rr_l = rr_l, rr_u = rr_u,
       pc = pc, pc_l = pc_l, pc_u = pc_u)
}

# 4) 单暴露 GAM：返回 beta/se/p + Δ=0.1/IQR/SD 的 RR/PC/CI
fit_one_exposure <- function(data, exposure, ns_value) {
  
  # 暴露不在数据里就跳过
  if (!exposure %in% names(data)) return(NULL)
  
  # ---- 强制要求这些列必须存在（否则直接停或返回NULL）
  required <- c(
    "n","dead_year","total.population","GDP_avg_quintile",
    "urbanization_rate",
    "popolation_age_65_Ratio","Male_Ratio",
    "t2m_mean_summer","t2m_mean_winter","t2m_sd_summer","t2m_sd_winter",
    "climate_zone_name_code", exposure
  )
  miss <- setdiff(required, names(data))
  if (length(miss) > 0) {
    message("Skip exposure=", exposure, " because missing columns: ", paste(miss, collapse = ", "))
    return(NULL)
  }
  
  # 只取模型需要的列
  dat <- data[, required, drop = FALSE]
  dat<-dat%>%
    #group_by(dead_year) %>%
    mutate(urbanization_rate_quintile = ntile(urbanization_rate, 5))
  # 确保关键列类型合理
  dat$dead_year <- as.numeric(dat$dead_year)
  dat$GDP_avg_quintile <- as.factor(dat$GDP_avg_quintile)
  dat$urbanization_rate_quintile <- as.factor(dat$urbanization_rate_quintile)
  dat$climate_zone_name_code <- as.factor(dat$climate_zone_name_code)
  
  # offset
  dat$lpopulation <- log(dat$total.population)
  
  
  # GDP固定调整；暴露是GDP类变量则不调整GDP（避免共线/重复）
  gdp_like <- c("GDP_avg_quintile", "gdp_per_capita", "GDP_per_capita", "GDP_pc", "gdp_pc")
  adjust_gdp <- !(exposure %in% gdp_like)
  
  rhs <- paste0(
    exposure,
    " + ns(dead_year, ", ns_value, ")",
    if (adjust_gdp) " + as.factor(GDP_avg_quintile)" else "",
    " + offset(lpopulation)",
    "+as.factor(urbanization_rate_quintile)",
    "+popolation_age_65_Ratio","+Male_Ratio",
    " + t2m_mean_summer + t2m_mean_winter + t2m_sd_summer + t2m_sd_winter",
    " + s(climate_zone_name_code, bs='re')"
  )
  
  fml <- as.formula(paste0("n ~ ", rhs))
  
  # ---- 用 tryCatch 捕获 mgcv 的内部怪错，避免整个循环崩
  m <- tryCatch(
    mgcv::gam(
      fml,
      family = quasipoisson,
      data = dat,
      method = "REML"
    ),
    error = function(e) {
      message("GAM failed for exposure=", exposure, " | ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(m)) return(NULL)
  
  pt <- as.data.frame(summary(m)$p.table)
  if (!exposure %in% rownames(pt)) return(NULL)
  
  beta <- unname(pt[exposure, "Estimate"])
  se   <- unname(pt[exposure, "Std. Error"])
  pval <- unname(pt[exposure, "Pr(>|t|)"])
  
  # --- 在“用于该模型”的有效样本中计算 IQR/SD（避免CI不对齐） ---
  mf <- model.frame(m)
  x  <- mf[[exposure]]
  
  delta_10 <- 10
  delta_iqr <- as.numeric(IQR(x, na.rm = TRUE))
  delta_sd  <- as.numeric(stats::sd(x, na.rm = TRUE))
  if (!is.finite(delta_iqr) || delta_iqr <= 0) delta_iqr <- NA_real_
  if (!is.finite(delta_sd)  || delta_sd  <= 0) delta_sd  <- NA_real_
  
  e1   <- effect_with_delta(beta, se, 1)
  e10  <- effect_with_delta(beta, se, delta_10)
  eiqr <- if (is.na(delta_iqr)) NULL else effect_with_delta(beta, se, delta_iqr)
  esd  <- if (is.na(delta_sd))  NULL else effect_with_delta(beta, se, delta_sd)
  
  out <- tibble::tibble(
    exposure = exposure,
    beta = beta, se = se, p = pval,
    ns = ns_value,
    family = m$family$family,
    gcv_REML = summary(m)$sp.criterion[["REML"]],
    r2 = summary(m)$r.sq,
    dev_expl = summary(m)$dev.expl,
    x_IQR = delta_iqr,
    x_SD  = delta_sd,
    RR_1 = e1$rr, RR_1_l = e1$rr_l, RR_1_u = e1$rr_u,
    PC_1 = e1$pc, PC_1_l = e1$pc_l, PC_1_u = e1$pc_u,
    RR_10 = e10$rr, RR_10_l = e10$rr_l, RR_10_u = e10$rr_u,
    PC_10 = e10$pc, PC_10_l = e10$pc_l, PC_10_u = e10$pc_u
  )
  
  if (is.null(eiqr)) {
    out <- out %>% mutate(
      RR_IQR = NA_real_, RR_IQR_l = NA_real_, RR_IQR_u = NA_real_,
      PC_IQR = NA_real_, PC_IQR_l = NA_real_, PC_IQR_u = NA_real_
    )
  } else {
    out <- out %>% mutate(
      RR_IQR = eiqr$rr, RR_IQR_l = eiqr$rr_l, RR_IQR_u = eiqr$rr_u,
      PC_IQR = eiqr$pc, PC_IQR_l = eiqr$pc_l, PC_IQR_u = eiqr$pc_u
    )
  }
  
  if (is.null(esd)) {
    out <- out %>% mutate(
      RR_SD = NA_real_, RR_SD_l = NA_real_, RR_SD_u = NA_real_,
      PC_SD = NA_real_, PC_SD_l = NA_real_, PC_SD_u = NA_real_
    )
  } else {
    out <- out %>% mutate(
      RR_SD = esd$rr, RR_SD_l = esd$rr_l, RR_SD_u = esd$rr_u,
      PC_SD = esd$pc, PC_SD_l = esd$pc_l, PC_SD_u = esd$pc_u
    )
  }
  
  out <- out %>% mutate(
    PCCI_1 = paste0(round(PC_1,2), " (", round(PC_1_l,2), ", ", round(PC_1_u,2), ")"),
    PCCI_10 = paste0(round(PC_10,2), " (", round(PC_10_l,2), ", ", round(PC_10_u,2), ")"),
    PCCI_IQR = ifelse(is.na(PC_IQR), NA_character_,
                      paste0(round(PC_IQR,2), " (", round(PC_IQR_l,2), ", ", round(PC_IQR_u,2), ")")),
    PCCI_SD  = ifelse(is.na(PC_SD), NA_character_,
                      paste0(round(PC_SD,2), " (", round(PC_SD_l,2), ", ", round(PC_SD_u,2), ")"))
  )
  
  out
}

# 5) 批量跑：4结局 × 全部暴露（直接关联 + FDR）

all_res <- list()

for (dname in names(datasets)) {
  dat <- datasets[[dname]]
  dat<-dat %>% as.data.frame()
  res_d <- lapply(urban_environment_vars, function(x) fit_one_exposure(dat, x, ns_value))
  res_d <- bind_rows(res_d)
  
  if (nrow(res_d) == 0) next
  
  # 每个结局内部做 FDR（BH）
  res_d$q_fdr <- p.adjust(res_d$p, method = "BH")
  res_d$disease <- dname
  
  all_res[[dname]] <- res_d
}

result_direct_urban_environment <- bind_rows(all_res)
# write.csv(result_direct_urban_environment,"E:/PKU_fanyi/PM2.5/result/result_direct/result_direct_urban_environment.csv",
#           fileEncoding = 'GBK',row.names = F)
write.csv(result_direct_urban_environment,"E:/Deskbook/Result_generate/1_Urban_environment/result_direct/result_direct_urban_environment.csv",
          fileEncoding = 'GBK',row.names = F)

####做二分类环境因素与结局二分类（按中位数）单暴露 GAM：High vs Low####

fit_one_exposure_median2group <- function(data, exposure, ns_value,
                                          min_n = 30, min_each_group = 10) {
  
  if (!exposure %in% names(data)) return(NULL)
  
  required <- c(
    "n","dead_year","total.population","GDP_avg_quintile",
    "urbanization_rate",
    "popolation_age_65_Ratio","Male_Ratio",
    "t2m_mean_summer","t2m_mean_winter","t2m_sd_summer","t2m_sd_winter",
    "climate_zone_name_code", exposure
  )
  miss <- setdiff(required, names(data))
  if (length(miss) > 0) {
    message("Skip exposure=", exposure, " because missing columns: ", paste(miss, collapse = ", "))
    return(NULL)
  }
  
  dat <- data[, required, drop = FALSE]
  dat<- dat %>%
    mutate(
      urbanization_rate_quintile = as.factor(dplyr::ntile(urbanization_rate, 5)))
  dat$dead_year <- as.numeric(dat$dead_year)
  dat$GDP_avg_quintile <- as.factor(dat$GDP_avg_quintile)
  dat$climate_zone_name_code <- as.factor(dat$climate_zone_name_code)
  dat$lpopulation <- log(dat$total.population)
  
  # complete-case，避免 na.omit 后变 0 行造成 mgcv 怪错
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]
  if (nrow(dat) < min_n) return(NULL)
  
  x <- dat[[exposure]]
  if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) return(NULL)
  
  # 中位数（用有效样本计算）
  med <- as.numeric(stats::median(x, na.rm = TRUE))
  
  # 二分：High(>median) vs Low(<=median)
  # 若很多值正好等于median，会导致 High 太少；这里用 >median 更常见
  dat$x_bin <- ifelse(x > med, 1L, 0L)
  
  n0 <- sum(dat$x_bin == 0, na.rm = TRUE)
  n1 <- sum(dat$x_bin == 1, na.rm = TRUE)
  if (n0 < min_each_group || n1 < min_each_group) {
    message("Skip exposure=", exposure, " because groups too imbalanced: n0=", n0, ", n1=", n1)
    return(NULL)
  }
  
  # GDP重复调整保护
  gdp_like <- c("GDP_avg_quintile", "gdp_per_capita", "GDP_per_capita", "GDP_pc", "gdp_pc")
  adjust_gdp <- !(exposure %in% gdp_like)
  
  rhs <- paste0(
    "x_bin",
    " + ns(dead_year, ", ns_value, ")",
    if (adjust_gdp) " + as.factor(GDP_avg_quintile)" else "",
    " + offset(lpopulation)",
    "+as.factor(urbanization_rate_quintile)",
    "+popolation_age_65_Ratio","+Male_Ratio",
    " + t2m_mean_summer + t2m_mean_winter + t2m_sd_summer + t2m_sd_winter",
    " + s(climate_zone_name_code, bs='re')"
  )
  
  fml <- as.formula(paste0("n ~ ", rhs))
  
  m <- tryCatch(
    mgcv::gam(
      fml,
      family = quasipoisson,
      data = dat,
      method = "REML"
    ),
    error = function(e) {
      message("GAM(median2group) failed for exposure=", exposure, " | ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(m)) return(NULL)
  
  pt <- as.data.frame(summary(m)$p.table)
  if (!"x_bin" %in% rownames(pt)) return(NULL)
  
  beta <- unname(pt["x_bin", "Estimate"])
  se   <- unname(pt["x_bin", "Std. Error"])
  pval <- unname(pt["x_bin", "Pr(>|t|)"])
  z <- qnorm(0.975)
  
  RR   <- exp(beta)
  RR_l <- exp(beta - z * se)
  RR_u <- exp(beta + z * se)
  
  tibble(
    exposure = exposure,
    cut_type = "median",
    cut_value = med,
    group_def = "High(>median) vs Low(<=median)",
    n_used = nrow(dat),
    n_low = n0,
    n_high = n1,
    beta = beta,
    se = se,
    p = pval,
    RR = RR, RR_l = RR_l, RR_u = RR_u,
    PC = (RR - 1) * 100,
    PC_l = (RR_l - 1) * 100,
    PC_u = (RR_u - 1) * 100
  )
}

# 8) 批量跑：raw + score 的 median 二分关联

res_2group_list <- list()

for (dname in names(datasets)) {
  cat("Processing disease:", dname, "\n")
  cat("Available datasets:", paste(names(datasets), collapse = ", "), "\n")
  dat <- datasets[[dname]]
  dat<-dat %>% as.data.frame()
  res_d <- lapply(urban_environment_vars, function(x){
    fit_one_exposure_median2group(dat, x, ns_value)
  }) %>% bind_rows()
  
  if (nrow(res_d) == 0) next
  
  res_d$q_fdr <- p.adjust(res_d$p, method = "BH")
  res_d$disease <- dname
  res_d$feature_type <- "binary_median"
  
  res_2group_list[[dname]] <- res_d
}

result_direct_urban_environment_2group <- bind_rows(res_2group_list)
# write.csv(result_direct_urban_environment_2group,"E:/PKU_fanyi/PM2.5/result/result_direct/result_direct_urban_environment_2group.csv",
#           fileEncoding = 'GBK',row.names = F)
write.csv(result_direct_urban_environment_2group,"E:/Deskbook/Result_generate/1_Urban_environment/result_direct/result_direct_urban_environment_2group.csv",
          fileEncoding = 'GBK',row.names = F)
#整理二分类的结果
result_direct_urban_environment_2group <-fread("E:/Deskbook/Result_generate/1_Urban_environment/result_direct/result_direct_urban_environment_2group.csv")
names(result_direct_urban_environment_2group)
table(result_direct_urban_environment_2group$exposure)
table(result_direct_urban_environment_2group$disease)
#
# 普通指标顺序
indicator_order <- c(
  "water_wetland",
  "forest_shrub_grassland",
  "impervious",
  "Residential_per_totalArea",
  "Workplace_per_totalArea",
  "Industrial_production_per_totalArea",
  "Education_culture_per_totalArea",
  "Healthcare_per_totalArea",
  "Government_per_totalArea",
  "Commerce_services_per_totalArea",
  "Leisure_tourism_per_totalArea",
  "Logistics_delivery_per_totalArea",
  "Urban_public_transit_per_totalArea",
  "Intercity_passenger_hub_per_totalArea",
  "Motorized_transport_support_per_totalArea",
  "total_len_km_per_totalArea",
  "major_len_km_per_totalArea",
  "minor_len_km_per_totalArea"
)

# PCA 指标顺序
pca_order <- c(
  "score_LandCover_pca",
  "score_Functional_Service_Density_pca",
  "score_Transport_Facility_Density_pca",
  "score_Road_Network_Density_pca",
  "score_total_pca_allvars"
)

exposure_order <- c(indicator_order, pca_order)

# P 值格式：P < 0.001 用 ×10^-4 这种格式；其他保留原始数值
format_p <- function(p) {
  sapply(p, function(x) {
    if (is.na(x)) {
      return(NA_character_)
    }
    
    if (x < 0.001) {
      sci <- formatC(x, format = "e", digits = 2)
      sci_split <- strsplit(sci, "e")[[1]]
      
      base <- as.numeric(sci_split[1])
      exp_num <- as.integer(sci_split[2])
      
      return(paste0(sprintf("%.2f", base), "×10^", exp_num))
    } else {
      return(as.character(x))
    }
  })
}

# 筛选 total 和 accident
dt <- copy(result_direct_urban_environment_2group)

dt2 <- dt[
  disease %in% c("total", "accident") &
    exposure %in% exposure_order
]

# 检查是否存在同一 exposure-disease 多行
dup_check <- dt2[, .N, by = .(exposure, disease)][N > 1]

if (nrow(dup_check) > 0) {
  print(dup_check)
  stop("存在同一 exposure 和 disease 对应多行结果，请检查是否需要进一步按 cut_type 或 group_def 筛选。")
}

# 重新整理方向：如果 RR > 1，则 beta 取反，使 RR_new < 1
dt2[, beta_new := beta]

dt2[exp(beta_new) > 1, beta_new := -beta_new]

# 根据 beta_new 和 se 重新计算 RR、PC 和 95% CI
dt2[, RR_new := exp(beta_new)]
dt2[, RR_l_new := exp(beta_new - 1.96 * se)]
dt2[, RR_u_new := exp(beta_new + 1.96 * se)]

dt2[, PC_new := (RR_new - 1) * 100]
dt2[, PC_l_new := (RR_l_new - 1) * 100]
dt2[, PC_u_new := (RR_u_new - 1) * 100]

# 格式化 PC 和 95%CI，保留两位小数
dt2[, PC_95CI := sprintf("%.2f (%.2f, %.2f)", PC_new, PC_l_new, PC_u_new)]

# 格式化 P 值
dt2[, P_value := format_p(p)]

# 转宽表
table_2group <- dcast(
  dt2,
  exposure ~ disease,
  value.var = c("PC_95CI", "P_value")
)

# 按指定 exposure 顺序排列
table_2group[, exposure := factor(exposure, levels = exposure_order)]
table_2group <- table_2group[order(exposure)]

# 最终列顺序：total 的 PC95%CI 和 P，然后 accident 的 PC95%CI 和 P
table_2group_final <- table_2group[
  ,
  .(
    exposure = as.character(exposure),
    total_PC_95CI = PC_95CI_total,
    total_P = P_value_total,
    accident_PC_95CI = PC_95CI_accident,
    accident_P = P_value_accident
  )
]

# 导出结果
outfile <- "E:/Deskbook/Result_generate/1_Urban_environment/result_direct/result_direct_urban_environment_2group_2.csv"

write.csv(
  table_2group_final,
  outfile,fileEncoding = 'GBK',
  row.names = F
)

# 查看结果
print(table_2group_final)
####直接非线性关联####


# 1) 参数
ns_value <- 4         # 与线性模型一致
k_smooth <- 3         # 比 k=3 更合理；你可做敏感性：6/8/10
grid_n   <- 100       # 每条曲线输出点数
q_lo     <- 0.01      # x网格下界分位数，避免极端值
q_hi     <- 0.99
z        <- qnorm(0.975)


# 2) 工具函数

mode_fun <- function(v) {
  v <- v[!is.na(v)]
  if (length(v) == 0) return(NA)
  ux <- unique(v)
  ux[which.max(tabulate(match(v, ux)))]
}

# 构造dat_cc：包含原始 dead_year、lpopulation（而不是 ns(...) / offset(...)）
prep_dat_cc <- function(data, urban_environment_var_) {
  required <- c(
    "n","dead_year","total.population","GDP_avg_quintile",
    "urbanization_rate",
    "popolation_age_65_Ratio","Male_Ratio",
    "t2m_mean_summer","t2m_mean_winter","t2m_sd_summer","t2m_sd_winter",
    "climate_zone_name_code", urban_environment_var_
  )
  miss <- setdiff(required, names(data))
  if (length(miss) > 0) return(NULL)
  
  dat_cc <- data %>%
    mutate(
      urbanization_rate_quintile = as.factor(dplyr::ntile(urbanization_rate, 5)),
      dead_year = as.numeric(dead_year),
      GDP_avg_quintile = as.factor(GDP_avg_quintile),
      climate_zone_name_code = droplevels(as.factor(climate_zone_name_code)),
      lpopulation = log(total.population)
    ) %>%
    dplyr::select(
      n, dead_year, total.population, lpopulation,
      GDP_avg_quintile,urbanization_rate_quintile,
      popolation_age_65_Ratio,Male_Ratio,
      t2m_mean_summer, t2m_mean_winter, t2m_sd_summer, t2m_sd_winter,
      climate_zone_name_code,
      all_of(urban_environment_var_)
    ) %>%
    na.omit()
  
  if (nrow(dat_cc) < 30) return(NULL)
  if (nlevels(dat_cc$climate_zone_name_code) < 2) return(NULL)
  
  # urban_environment_var_ 变常数也跳过
  if (sd(dat_cc[[urban_environment_var_]], na.rm = TRUE) == 0) return(NULL)
  
  dat_cc
}


# 3) 主循环：导出非线性曲线结果

result_direct_nonlinear <- tibble()

for (dname in names(datasets)) {
  cat("Processing disease:", dname, "\n")
  cat("Available datasets:", paste(names(datasets), collapse = ", "), "\n")
  dat <- datasets[[dname]]
  dat<-dat %>% as.data.frame()
  
  for (urban_environment_var_ in urban_environment_vars) {
    
    if (!urban_environment_var_ %in% names(dat)) next
    
    dat_cc <- prep_dat_cc(dat, urban_environment_var_)
    if (is.null(dat_cc)) next
    
    # GDP重复调整保护（与线性一致）
    gdp_like <- c("GDP_avg_quintile","gdp_per_capita","GDP_per_capita","GDP_pc","gdp_pc")
    adjust_gdp <- !(urban_environment_var_ %in% gdp_like)
    
    # 公式：非线性 s(urban_environment_var_) + ns(year) + offset + met + random effect
    fml <- as.formula(paste0(
      "n ~ s(", urban_environment_var_, ", k=", k_smooth, ")",
      " + ns(dead_year, ", ns_value, ")",
      if (adjust_gdp) " + as.factor(GDP_avg_quintile)" else "",
      " + offset(lpopulation)",
      "+as.factor(urbanization_rate_quintile)",
      "+popolation_age_65_Ratio","+Male_Ratio",
      " + t2m_mean_summer + t2m_mean_winter + t2m_sd_summer + t2m_sd_winter",
      " + s(climate_zone_name_code, bs='re')"
    ))
    
    model <- tryCatch(
      gam(fml, family = quasipoisson, data = dat_cc, method = "REML"),
      error = function(e) {
        message("GAM failed: ", dname, " | ", urban_environment_var_, " | ", conditionMessage(e))
        return(NULL)
      }
    )
    if (is.null(model)) next
    
    # 预测网格（从 dat_cc 的 urban_environment_var_ 分布来）
    x_vec <- dat_cc[[urban_environment_var_]]
    x_grid <- seq(quantile(x_vec, q_lo, na.rm = TRUE),
                  quantile(x_vec, q_hi, na.rm = TRUE),
                  length.out = grid_n)
    
    # 参考值：中位数（相对中位数RR=1更容易解释）
    x_ref <- as.numeric(median(x_vec, na.rm = TRUE))
    
    # newdata：从 dat_cc 复制一行（包含原始 dead_year/lpopulation 等）
    newdata <- dat_cc[rep(1, grid_n), , drop = FALSE]
    newdata[[urban_environment_var_]] <- x_grid
    
    # 固定协变量（代表值）
    newdata$dead_year <- median(dat_cc$dead_year, na.rm = TRUE)
    newdata$total.population <- median(dat_cc$total.population, na.rm = TRUE)
    newdata$lpopulation <- log(newdata$total.population)
    
    newdata$popolation_age_65_Ratio <- mean(dat_cc$popolation_age_65_Ratio, na.rm = TRUE)
    newdata$Male_Ratio <- mean(dat_cc$Male_Ratio, na.rm = TRUE)
    
    
    newdata$t2m_mean_summer <- mean(dat_cc$t2m_mean_summer, na.rm = TRUE)
    newdata$t2m_mean_winter <- mean(dat_cc$t2m_mean_winter, na.rm = TRUE)
    newdata$t2m_sd_summer   <- mean(dat_cc$t2m_sd_summer, na.rm = TRUE)
    newdata$t2m_sd_winter   <- mean(dat_cc$t2m_sd_winter, na.rm = TRUE)
    
    if (adjust_gdp) newdata$GDP_avg_quintile <- mode_fun(dat_cc$GDP_avg_quintile)
    newdata$urbanization_rate_quintile <- mode_fun(dat_cc$urbanization_rate_quintile)
    newdata$climate_zone_name_code <- mode_fun(dat_cc$climate_zone_name_code)
    newdata$climate_zone_name_code <- droplevels(as.factor(newdata$climate_zone_name_code))
    
    # 取平滑项的 term 与 se
    pr <- predict(model, newdata = newdata, type = "terms", se.fit = TRUE)
    
    # 找到该 urban_environment_var_ 的平滑项列名（通常是 "s(urban_environment_var_)"）
    term_name <- grep(paste0("^s\\(", urban_environment_var_, "\\)$"), colnames(pr$fit), value = TRUE)
    if (length(term_name) != 1) {
      # 有时会出现带冒号/不同命名，备用匹配
      term_name <- grep(paste0("s\\(", urban_environment_var_, "\\)"), colnames(pr$fit), value = TRUE)
      if (length(term_name) < 1) next
      term_name <- term_name[1]
    }
    
    fit_term <- as.numeric(pr$fit[, term_name])
    se_term  <- as.numeric(pr$se.fit[, term_name])
    
    # 参考点预测（单点）
    new_ref <- newdata[1, , drop = FALSE]
    new_ref[[urban_environment_var_]] <- x_ref
    pr_ref <- predict(model, newdata = new_ref, type = "terms", se.fit = TRUE)
    fit_ref <- as.numeric(pr_ref$fit[, term_name])
    
    # 相对参考值的 logRR 曲线
    logRR   <- fit_term - fit_ref
    logRR_l <- logRR - z * se_term
    logRR_u <- logRR + z * se_term
    
    RR   <- exp(logRR)
    RR_l <- exp(logRR_l)
    RR_u <- exp(logRR_u)
    
    result_df <- tibble(
      disease = dname,
      urban_environment_var = urban_environment_var_,
      x = x_grid,
      x_ref = x_ref,
      term = fit_term,
      se_term = se_term,
      logRR = logRR,
      RR = RR, RR_l = RR_l, RR_u = RR_u,
      PC = (RR - 1) * 100,
      PC_l = (RR_l - 1) * 100,
      PC_u = (RR_u - 1) * 100,
      ns = ns_value,
      k = k_smooth,
      n_used = nrow(dat_cc)
    )
    
    result_direct_nonlinear <- bind_rows(result_direct_nonlinear, result_df)
  }
}
# write.csv(result_direct_nonlinear,"E:/PKU_fanyi/PM2.5/result/result_direct/result_direct_urban_environment_nonlinear.csv",
#           fileEncoding = 'GBK',row.names = F)
write.csv(result_direct_nonlinear,"E:/Deskbook/Result_generate/1_Urban_environment/result_direct/result_direct_urban_environment_nonlinear.csv",
          fileEncoding = 'GBK',row.names = F)
#长期修饰作用####
####先做污染物与结局的直接关联####
#统一进行标化
#
ranges <- list(
  c("sum_heat_intensity0__Tmean_90_10", "mean_BC")
)
#
template_dt <- datasets$total
vars_to_std <- unique(unlist(lapply(ranges, \(rg) col_range(template_dt, rg[1], rg[2]))))
#
# ---- 标化函数：z-score（生成新列 z_*** 或覆盖原列）
standardize_one_df <- function(df, vars, prefix = "z_", overwrite = FALSE,
                               skip_binary = TRUE, min_n = 30){
  
  df <- data.table::as.data.table(df)
  
  # 只保留本 df 中存在的列
  v_ok <- intersect(vars, names(df))
  if (length(v_ok) == 0) return(df)
  
  # 只保留 numeric（兼容老 data.table：用 with=FALSE）
  is_num <- vapply(df[, v_ok, with = FALSE], is.numeric, logical(1))
  v_ok <- v_ok[is_num]
  if (length(v_ok) == 0) return(df)
  
  # 可选：跳过 0/1 变量
  if (skip_binary) {
    is_bin <- vapply(v_ok, function(v){
      x <- df[[v]]
      x <- x[is.finite(x)]
      length(x) > 0 && all(x %in% c(0,1))
    }, logical(1))
    v_ok <- v_ok[!is_bin]
    if (length(v_ok) == 0) return(df)
  }
  
  for (v in v_ok) {
    x  <- df[[v]]
    x2 <- x[is.finite(x)]
    if (length(x2) < min_n) next
    
    m <- mean(x2)
    s <- stats::sd(x2)
    if (!is.finite(s) || s == 0) next
    
    z <- (x - m) / s
    
    if (overwrite) {
      df[[v]] <- z
    } else {
      df[[paste0(prefix, v)]] <- z
    }
  }
  
  df
}

# ---- 批量处理 datasets
datasets <- purrr::imap(datasets, ~ standardize_one_df(.x, vars_to_std,
                                                       prefix = "z_", overwrite = FALSE,
                                                       skip_binary = TRUE))

#选暴露因素
ranges <- list(
  c("z_sum_heat_intensity0__Tmean_90_10", "z_mean_BC")
)
#
template_dt <- datasets$total
exposure_vars <- unique(unlist(lapply(ranges, \(rg) col_range(template_dt, rg[1], rg[2]))))

fit_one_exposure <- function(data, exposure,
                             ns_value = 4,
                             delta_1 = 1,
                             min_n = 30,
                             treat_z_prefix_as_scaled = TRUE) {
  
  if (!exposure %in% names(data)) return(NULL)
  
  required <- c(
    "n","dead_year","total.population","GDP_avg_quintile",
    "urbanization_rate",
    "popolation_age_65_Ratio","Male_Ratio",
    "t2m_mean_summer","t2m_mean_winter","t2m_sd_summer","t2m_sd_winter",
    "climate_zone_name_code", exposure
  )
  miss <- setdiff(required, names(data))
  if (length(miss) > 0) return(NULL)
  
  dat <- data[, required, drop = FALSE] %>%
    mutate(
      dead_year = as.numeric(dead_year),
      GDP_avg_quintile = as.factor(GDP_avg_quintile),
      urbanization_rate_quintile = as.factor(ntile(urbanization_rate, 5)),
      climate_zone_name_code = droplevels(as.factor(climate_zone_name_code)),
      lpopulation = log(total.population)
    )
  
  dat_cc <- dat[complete.cases(dat), , drop = FALSE]
  if (nrow(dat_cc) < min_n) return(NULL)
  if (nlevels(dat_cc$climate_zone_name_code) < 2) return(NULL)
  
  x0 <- dat_cc[[exposure]]
  x1 <- x0[is.finite(x0)]
  if (length(x1) == 0) return(NULL)
  
  is_binary <- all(x1 %in% c(0, 1)) && length(unique(x1)) <= 2
  is_scaled <- (treat_z_prefix_as_scaled && grepl("^z_", exposure))
  if (!is_binary && sd(x1, na.rm = TRUE) == 0) return(NULL)
  
  fml <- as.formula(paste0(
    "n ~ ", exposure,
    " + splines::ns(dead_year, ", ns_value, ")",
    " + as.factor(GDP_avg_quintile)",
    " + offset(lpopulation)",
    "+as.factor(urbanization_rate_quintile)",
    "+popolation_age_65_Ratio","+Male_Ratio",
    " + t2m_mean_summer + t2m_mean_winter + t2m_sd_summer + t2m_sd_winter",
    " + s(climate_zone_name_code, bs='re')"
  ))
  
  m <- tryCatch(
    mgcv::gam(fml, family = quasipoisson, data = dat_cc, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(m)) return(NULL)
  
  pt <- as.data.frame(summary(m)$p.table)
  if (!exposure %in% rownames(pt)) return(NULL)
  
  beta <- unname(pt[exposure, "Estimate"])
  se   <- unname(pt[exposure, "Std. Error"])
  pval <- unname(pt[exposure, "Pr(>|t|)"])
  
  e1 <- effect_with_delta(beta, se, delta_1)
  
  tibble(
    exposure = exposure,
    exposure_type = case_when(
      is_binary ~ "binary_0_1",
      is_scaled ~ "z_scaled",
      TRUE      ~ "continuous_raw"
    ),
    beta = beta, se = se, p = pval,
    ns = ns_value,
    r2 = summary(m)$r.sq,
    dev_expl = summary(m)$dev.expl,
    RR_1   = e1$rr,   RR_1_l = e1$rr_l, RR_1_u = e1$rr_u,
    PC_1   = e1$pc,   PC_1_l = e1$pc_l, PC_1_u = e1$pc_u
  )
}
# 2) 批量跑所有暴露：feature_type 不再 match.arg 限制
classify_exposure_auto <- function(x){
  case_when(
    grepl("^z_composite_event_", x) ~ "composite_event_scaled",
    grepl("^composite_event_", x)   ~ "composite_event_raw",
    grepl("^z_composite_intensity_", x) ~ "composite_intensity_scaled",
    grepl("^composite_intensity_", x)   ~ "composite_intensity_raw",
    grepl("^z_score_", x) ~ "score_scaled",
    grepl("^score_", x)   ~ "score_raw",
    grepl("TAP_O3|O3", x, ignore.case = TRUE) ~ "o3",
    grepl("heat_intensity", x, fixed = TRUE) ~ "heat_intensity",
    grepl("heat_event", x, fixed = TRUE) | grepl("heat_wave", x, fixed = TRUE) ~ "heat_event_days",
    grepl("PM2.5", x, fixed = TRUE) & grepl("event_intensity0", x, fixed = TRUE) ~ "pm_event_intensity",
    grepl("PM2.5", x, fixed = TRUE) & grepl("event0", x, fixed = TRUE) ~ "pm_event_days",
    grepl("mean_PM2.5", x, fixed = TRUE) | grepl("mean_SO4_|mean_NO3_|mean_NH4_|mean_OM_|mean_BC_", paste0("_", x, "_")) ~ "pm25_component",
    TRUE ~ "other"
  )
}

run_exposure_direct <- function(datasets, exposure_vars,
                                ns_value = 4,
                                delta_1 = 1,
                                feature_type = "auto",
                                verbose = TRUE) {
  
  res_list <- list()
  
  for (dname in names(datasets)) {
    if (verbose) cat("Processing disease:", dname, "\n")
    
    dat <- as.data.frame(datasets[[dname]])
    vars_use <- intersect(exposure_vars, names(dat))
    if (length(vars_use) == 0) next
    
    res_d <- lapply(vars_use, function(v) {
      fit_one_exposure(dat, v, ns_value = ns_value, delta_1 = delta_1)
    }) %>% bind_rows()
    
    if (is.null(res_d) || nrow(res_d) == 0) next
    
    res_d$q_fdr <- p.adjust(res_d$p, method = "BH")
    res_d$disease <- dname
    
    res_d$feature_type <- if (identical(feature_type, "auto")) {
      classify_exposure_auto(res_d$exposure)
    } else as.character(feature_type)
    
    res_list[[dname]] <- res_d
  }
  
  bind_rows(res_list)
}

# 4) 用 beta 加权生成 score/composite，并写回 datasets（同时生成 z_ 版本）
add_beta_weighted_scores <- function(datasets,
                                     beta_table,
                                     score_specs,
                                     make_z = TRUE,
                                     z_prefix = "z_",
                                     min_n_z = 30,
                                     normalize_by_absw = FALSE,
                                     na_to_zero = FALSE,
                                     weight_filter = NULL) {
  
  score_specs <- score_specs[vapply(score_specs, function(v) length(v) > 0, logical(1))]
  score_names <- names(score_specs)
  
  weights_used <- vector("list", length(datasets)); names(weights_used) <- names(datasets)
  
  for (dn in names(datasets)) {
    dat <- as.data.frame(datasets[[dn]])
    weights_used[[dn]] <- list()
    
    for (sc in score_names) {
      
      vars0 <- unique(score_specs[[sc]])
      vars  <- intersect(vars0, names(dat))
      
      if (length(vars) == 0) {
        dat[[sc]] <- NA_real_
        if (make_z) dat[[paste0(z_prefix, sc)]] <- NA_real_
        weights_used[[dn]][[sc]] <- tibble(exposure = character(), beta = numeric())
        next
      }
      
      w <- beta_table %>%
        filter(disease == dn, exposure %in% vars) %>%
        dplyr::select(exposure, beta, p, q_fdr)
      
      if (!is.null(weight_filter)) {
        w <- w %>% filter(!!weight_filter)
      }
      
      if (nrow(w) == 0) {
        dat[[sc]] <- NA_real_
        if (make_z) dat[[paste0(z_prefix, sc)]] <- NA_real_
        weights_used[[dn]][[sc]] <- w %>% select(exposure, beta)
        next
      }
      
      X <- as.matrix(dat[, w$exposure, drop = FALSE])
      if (na_to_zero) X[is.na(X)] <- 0
      
      b <- w$beta
      if (normalize_by_absw) {
        denom <- sum(abs(b), na.rm = TRUE)
        if (is.finite(denom) && denom > 0) b <- b / denom
      }
      
      s <- as.numeric(X %*% b)
      dat[[sc]] <- s
      
      if (make_z) {
        s2 <- s[is.finite(s)]
        if (length(s2) >= min_n_z && stats::sd(s2) > 0) {
          dat[[paste0(z_prefix, sc)]] <- (s - mean(s2)) / stats::sd(s2)
        } else {
          dat[[paste0(z_prefix, sc)]] <- NA_real_
        }
      }
      
      weights_used[[dn]][[sc]] <- w %>% dplyr::select(exposure, beta)
    }
    
    datasets[[dn]] <- dat
  }
  
  list(datasets = datasets, weights_used = weights_used)
}

# 5) 你的“基础暴露变量”（把你列的都放进来）
base_exposure_vars <- c(
  # O3
  "z_mean_TAP_O3",
  
  # PM2.5 + 组分
  "z_mean_PM2.5_TAP","z_mean_SO4","z_mean_NO3","z_mean_NH4","z_mean_OM","z_mean_BC",
  
  # 热强度（90_10, 95_5, 92.5_7.5）
  "z_sum_heat_intensity0__Tmean_90_10","z_sum_heat_intensity_mean2__Tmean_90_10","z_sum_heat_intensity_sum2__Tmean_90_10",
  "z_sum_heat_intensity0__Tmean_95_5","z_sum_heat_intensity_mean2__Tmean_95_5","z_sum_heat_intensity_sum2__Tmean_95_5",
  "z_sum_heat_intensity0__Tmean_92.5_7.5","z_sum_heat_intensity_mean2__Tmean_92.5_7.5","z_sum_heat_intensity_sum2__Tmean_92.5_7.5",
  
  # PM2.5 超污染事件强度
  "z_sum_PM2.5_event_intensity0","z_sum_PM2.5_TAP_event_intensity0",
  "z_sum_lag01PM2.5_event_intensity0","z_sum_lag01PM2.5_TAP_event_intensity0",
  
  # 热事件天数（95_5, 92.5_7.5；只用 event0 + wave01）
  "z_sum_heat_event0__Tmean_95_5","z_sum_heat_wave01__Tmean_95_5",
  "z_sum_heat_event0__Tmean_92.5_7.5","z_sum_heat_wave01__Tmean_92.5_7.5",
  
  # PM2.5 超污染事件天数
  "z_sum_PM2.5_event0","z_sum_PM2.5_TAP_event0",
  "z_sum_lag01PM2.5_event0","z_sum_lag01PM2.5_TAP_event0"
)

# 6) 按你要求构建：PM25mix + composite_event_* + composite_intensity_*
#    - composite_event：2(阈值) × 4(PM事件天数) = 8
#    - composite_intensity：3(阈值) × 4(PM事件强度) = 12

# (1) PM2.5 mix
pm25mix_vars <- c("z_mean_PM2.5_TAP","z_mean_PM2.5_TAP",
                  "z_mean_SO4","z_mean_NO3","z_mean_NH4","z_mean_OM","z_mean_BC")

# (2) 复合事件天数：热事件（按阈值分组） + PM事件天数（逐个）
heat_event_vars <- c(
  "z_sum_heat_event0__Tmean_95_5",
  "z_sum_heat_wave01__Tmean_95_5",
  "z_sum_heat_event0__Tmean_92.5_7.5",
  "z_sum_heat_wave01__Tmean_92.5_7.5"
)

pm_event_days_vars <- c(
  "z_sum_PM2.5_event0",
  "z_sum_PM2.5_TAP_event0",
  "z_sum_lag01PM2.5_event0",
  "z_sum_lag01PM2.5_TAP_event0"
)

composite_event_specs <- list()
for (hv in heat_event_vars) {
  for (pmv in pm_event_days_vars) {
    heat_tag <- sub("^z_", "", hv)
    pm_tag   <- sub("^z_", "", pmv)
    nm <- paste0("composite_event_", heat_tag, "_", pm_tag)
    composite_event_specs[[nm]] <- c(hv, pmv)   # 注意：这里是“各自”，只放 2 个变量
  }
}

# (3) 复合强度：热强度（按阈值分组） + PM事件强度（逐个）
heat_intensity_vars <- c(
  "z_sum_heat_intensity0__Tmean_90_10",
  "z_sum_heat_intensity_mean2__Tmean_90_10",
  "z_sum_heat_intensity_sum2__Tmean_90_10",
  
  "z_sum_heat_intensity0__Tmean_95_5",
  "z_sum_heat_intensity_mean2__Tmean_95_5",
  "z_sum_heat_intensity_sum2__Tmean_95_5",
  
  "z_sum_heat_intensity0__Tmean_92.5_7.5",
  "z_sum_heat_intensity_mean2__Tmean_92.5_7.5",
  "z_sum_heat_intensity_sum2__Tmean_92.5_7.5"
)

pm_event_intensity_vars <- c(
  "z_sum_PM2.5_event_intensity0",
  "z_sum_PM2.5_TAP_event_intensity0",
  "z_sum_lag01PM2.5_event_intensity0",
  "z_sum_lag01PM2.5_TAP_event_intensity0"
)

composite_intensity_specs <- list()
for (hv in heat_intensity_vars) {
  for (pmv in pm_event_intensity_vars) {
    heat_tag <- sub("^z_", "", hv)
    pm_tag   <- sub("^z_", "", pmv)
    nm <- paste0("composite_intensity_", heat_tag, "_", pm_tag)
    composite_intensity_specs[[nm]] <- c(hv, pmv)  # “各自”，只放 2 个变量
  }
}
# 汇总 score_specs（你要的全部指标）
score_specs <- c(
  list(score_PM25mix = pm25mix_vars),
  composite_event_specs,
  composite_intensity_specs
)

# 7) 一次性流程：先跑基础暴露 -> 用 beta 生成 composite -> 写回 datasets -> 再把新变量一起跑

# A) 先跑基础暴露，拿到每个 disease 的 beta（用来加权 composite）
result_base <- run_exposure_direct(
  datasets = datasets,
  exposure_vars = base_exposure_vars,
  ns_value = 4,
  delta_1 = 1,
  feature_type = "auto",
  verbose = TRUE
)

# B) 生成 composite/score，并写回 datasets（同时生成 z_ 版本：z_score_* / z_composite_*）
out_scores <- add_beta_weighted_scores(
  datasets = datasets,
  beta_table = result_base,
  score_specs = score_specs,
  make_z = TRUE,
  z_prefix = "z_",
  normalize_by_absw = FALSE,  # 你如果希望不同 composite 可比性更强，可改 TRUE
  na_to_zero = FALSE,
  weight_filter = NULL        # 如果只想用显著权重：quote(q_fdr <= 0.05)
)

datasets <- out_scores$datasets
weights_used <- out_scores$weights_used

# C) 把新生成的 composite/score 也加入要跑的暴露集合
new_vars_raw <- names(score_specs)
new_vars_z   <- paste0("z_", names(score_specs))
exposure_vars_all <- unique(c(base_exposure_vars, new_vars_raw, new_vars_z))

# D) 再跑一次：基础暴露 + 所有 composite/score（含 z_ 版本）
result_exposure_all <- run_exposure_direct(
  datasets = datasets,
  exposure_vars = exposure_vars_all,
  ns_value = 4,
  delta_1 = 1,
  feature_type = "auto",
  verbose = TRUE
)
# write.csv(result_direct_pollution,"E:/PKU_fanyi/PM2.5/result/result_long_modification/result_pollution_direct.csv",
#           fileEncoding = 'GBK',row.names = F)
write.csv(result_exposure_all,"E:/Deskbook/Result_generate/1_Urban_environment/result_long_modification/result_exposure_all.csv",
          fileEncoding = 'GBK',row.names = F)
####线性修饰作用####
exposure_vars_all <-c("z_mean_PM2.5_TAP","z_mean_SO4","z_mean_NO3",
                      "z_mean_NH4","z_mean_OM","z_mean_BC","score_PM25mix",
                      "z_sum_heat_intensity0__Tmean_92.5_7.5",
                      # "z_sum_heat_intensity_mean2__Tmean_92.5_7.5",
                      # "z_sum_heat_intensity_sum2__Tmean_92.5_7.5",
                      "z_sum_PM2.5_TAP_event_intensity0",
                      "z_sum_heat_event0__Tmean_92.5_7.5",
                      "z_sum_PM2.5_TAP_event0",
                      "composite_event_sum_heat_event0__Tmean_92.5_7.5_sum_PM2.5_TAP_event0",
                      "composite_intensity_sum_heat_intensity0__Tmean_92.5_7.5_sum_PM2.5_TAP_event_intensity0"
                      #"composite_intensity_sum_heat_intensity_mean2__Tmean_92.5_7.5_sum_PM2.5_TAP_event_intensity0",
                      #"composite_intensity_sum_heat_intensity_sum2__Tmean_92.5_7.5_sum_PM2.5_TAP_event_intensity0"
)
z <- qnorm(0.975)

# 生成城市因素网格：用有效样本分位数，避免score不在0-1时出错
make_urban_grid <- function(x, n = 40, q_lo = 0.01, q_hi = 0.99) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NULL)
  lo <- as.numeric(quantile(x, q_lo, na.rm = TRUE))
  hi <- as.numeric(quantile(x, q_hi, na.rm = TRUE))
  if (!is.finite(lo) || !is.finite(hi) || lo == hi) return(NULL)
  if (lo >= -1e-6 && hi <= 1 + 1e-6) return(seq(0, 1, length.out = n))
  seq(lo, hi, length.out = n)
}

# ====== 关键：effect_with_delta 的“向量版” 
effect_with_delta_vec <- function(beta, se, delta = 1, z = qnorm(0.975)) {
  b <- beta * delta
  s <- se * abs(delta)
  l <- b - z * s
  u <- b + z * s
  rr   <- exp(b); rr_l <- exp(l); rr_u <- exp(u)
  pc   <- (rr   - 1) * 100
  pc_l <- (rr_l - 1) * 100
  pc_u <- (rr_u - 1) * 100
  tibble::tibble(RR = rr, RR_l = rr_l, RR_u = rr_u,
                 PC = pc, PC_l = pc_l, PC_u = pc_u)
}

# ====== 选择 delta 的规则（核心逻辑）
# effect_scale 取值： "auto" / "1" / "SD" / "IQR"
effect_scale <- "1"

pick_delta <- function(exposure_name, pol_x, effect_scale = c("auto","1","SD","IQR")) {
  effect_scale <- match.arg(effect_scale)
  
  pol_x <- pol_x[is.finite(pol_x)]
  if (!length(pol_x)) return(NULL)
  
  delta_iqr <- as.numeric(IQR(pol_x, na.rm = TRUE))
  delta_sd  <- as.numeric(sd(pol_x, na.rm = TRUE))
  if (!is.finite(delta_iqr) || delta_iqr <= 0) delta_iqr <- NA_real_
  if (!is.finite(delta_sd)  || delta_sd  <= 0) delta_sd  <- NA_real_
  
  # 如果是 z_ 标化变量：1 unit = 1 SD
  is_scaled <- grepl("^z_", exposure_name)
  
  if (is_scaled) {
    return(list(delta = 1, scale_used = "SD", pol_IQR = delta_iqr, pol_SD = delta_sd))
  }
  
  # 非 z_：按 effect_scale 决定 delta
  if (effect_scale == "1") {
    return(list(delta = 1, scale_used = "1", pol_IQR = delta_iqr, pol_SD = delta_sd))
  }
  if (effect_scale == "SD") {
    if (is.na(delta_sd)) return(NULL)
    return(list(delta = delta_sd, scale_used = "SD", pol_IQR = delta_iqr, pol_SD = delta_sd))
  }
  if (effect_scale == "IQR") {
    if (is.na(delta_iqr)) return(NULL)
    return(list(delta = delta_iqr, scale_used = "IQR", pol_IQR = delta_iqr, pol_SD = delta_sd))
  }
  
  # auto：优先 IQR，IQR 不可用则用 SD
  if (!is.na(delta_iqr)) {
    return(list(delta = delta_iqr, scale_used = "IQR", pol_IQR = delta_iqr, pol_SD = delta_sd))
  } else if (!is.na(delta_sd)) {
    return(list(delta = delta_sd, scale_used = "SD", pol_IQR = delta_iqr, pol_SD = delta_sd))
  } else {
    return(NULL)
  }
}

result_long_continuous_modification <- tibble::tibble()

for (dataset_name in names(datasets)) {
  
  data <- datasets[[dataset_name]] %>%
    dplyr::mutate(
      dead_year = as.numeric(dead_year),
      urbanization_rate_quintile = as.factor(dplyr::ntile(urbanization_rate, 5)),
      GDP_avg_quintile = as.factor(GDP_avg_quintile),
      climate_zone_name_code = droplevels(as.factor(climate_zone_name_code)),
      lpopulation = log(total.population)
    )
  
  for (exposure_ in exposure_vars_all) {
    
    if (!exposure_ %in% names(data)) next
    
    for (urban_ in urban_environment_vars) {
      
      if (!urban_ %in% names(data)) next
      
      fml <- as.formula(paste0(
        "n ~ ", exposure_, " * ", urban_,
        " + splines::ns(dead_year, 4)",
        " + as.factor(GDP_avg_quintile)",
        " + offset(lpopulation)",
        "+as.factor(urbanization_rate_quintile)",
        "+popolation_age_65_Ratio","+Male_Ratio",
        " + t2m_mean_summer + t2m_mean_winter + t2m_sd_summer + t2m_sd_winter",
        " + s(climate_zone_name_code, bs='re')"
      ))
      
      m <- tryCatch(
        mgcv::gam(fml, family = quasipoisson, data = data, na.action = na.omit, method = "REML"),
        error = function(e) NULL
      )
      if (is.null(m)) next
      
      mf <- model.frame(m)
      if (!exposure_ %in% names(mf) || !urban_ %in% names(mf)) next
      
      # 只处理连续变量（否则交互项不止一个系数）
      if (!is.numeric(mf[[exposure_]]) || !is.numeric(mf[[urban_]])) next
      
      pol_x <- mf[[exposure_]]
      urb_x <- mf[[urban_]]
      
      # 选择 delta：z_ -> delta=1 (per SD)；非 z_ -> 按 effect_scale
      dd <- pick_delta(exposure_, pol_x, effect_scale = effect_scale)
      if (is.null(dd)) next
      delta_use <- dd$delta
      effect_scale_used <- dd$scale_used
      
      grid <- make_urban_grid(urb_x)
      if (is.null(grid)) next
      
      cf <- coef(m)
      V  <- vcov(m)
      
      main_term <- exposure_
      int1 <- paste0(exposure_, ":", urban_)
      int2 <- paste0(urban_, ":", exposure_)
      int_term <- if (int1 %in% names(cf)) int1 else if (int2 %in% names(cf)) int2 else NA_character_
      if (is.na(int_term)) next
      if (!main_term %in% names(cf)) next
      
      b1 <- unname(cf[main_term])
      b2 <- unname(cf[int_term])
      
      v11 <- V[main_term, main_term]
      v22 <- V[int_term, int_term]
      v12 <- V[main_term, int_term]
      
      beta_per1 <- b1 + b2 * grid
      var_per1  <- v11 + (grid^2) * v22 + 2 * grid * v12
      se_per1   <- sqrt(pmax(var_per1, 0))
      
      # 统一用 effect_with_delta 计算：per delta_use（z_ 时 delta_use=1 即 per 1SD）
      eff <- effect_with_delta_vec(beta_per1, se_per1, delta = delta_use, z = z)
      
      out <- tibble::tibble(
        disease = dataset_name,
        exposure = exposure_,
        urban = urban_,
        x = grid,
        beta_per1 = beta_per1,
        se_per1 = se_per1,
        effect_scale = effect_scale_used,
        delta = delta_use,
      ) %>% dplyr::bind_cols(eff)
      
      result_long_continuous_modification <- dplyr::bind_rows(result_long_continuous_modification, out)
    }
  }
}
# write.csv(result_long_continuous_modification,"E:/PKU_fanyi/PM2.5/result/result_long_modification/result_pollution_continuous_modification.csv",
#           fileEncoding = 'GBK',row.names = F)
fwrite(result_long_continuous_modification,"E:/Deskbook/Result_generate/1_Urban_environment/result_long_modification/result_exposure_linear_continuous_modification.csv",
          row.names = F)
####二分类修饰作用####
# 结果表
result_long_cat_modification <- tibble()

for (dataset_name in names(datasets)) {
  
  data0 <- datasets[[dataset_name]] %>%
    mutate(
      dead_year = as.numeric(dead_year),
      urbanization_rate_quintile = as.factor(dplyr::ntile(urbanization_rate, 5)),
      GDP_avg_quintile = as.factor(GDP_avg_quintile),
      climate_zone_name_code = as.factor(climate_zone_name_code),
      lpopulation = log(total.population)
    )
  
  # 城市因素：二分类分层（按中位数）
  for (strate_factor in urban_environment_vars) {
    
    if (!strate_factor %in% names(data0)) next
    
    med <- median(data0[[strate_factor]], na.rm = TRUE)
    if (!is.finite(med)) next
    
    # 每个 strate_factor 单独生成 strata
    data <- data0 %>%
      mutate(.strata = ifelse(.data[[strate_factor]] > med, 1L, 0L))
    
    for (fenzu in 0:1) {
      
      data1 <- data %>% filter(.strata == fenzu)
      
      if (nrow(data1) < 50) next
      if (nlevels(droplevels(data1$climate_zone_name_code)) < 2) next
      
      for (exposure_ in exposure_vars_all) {
        
        if (!exposure_ %in% names(data1)) next
        
        needed <- c(
          "n","dead_year","GDP_avg_quintile","lpopulation",
          "urbanization_rate_quintile",
          "popolation_age_65_Ratio","Male_Ratio",
          "t2m_mean_summer","t2m_mean_winter","t2m_sd_summer","t2m_sd_winter",
          "climate_zone_name_code", exposure_
        )
        
        dat_mod <- data1[, needed, drop = FALSE] %>% na.omit()
        
        if (nrow(dat_mod) < 30) next
        if (!is.numeric(dat_mod[[exposure_]])) next
        if (sd(dat_mod[[exposure_]], na.rm = TRUE) == 0) next
        if (nlevels(droplevels(dat_mod$climate_zone_name_code)) < 2) next
        
        fml <- as.formula(paste0(
          "n ~ ", exposure_,
          " + splines::ns(dead_year, ", ns_value, ")",
          " + as.factor(GDP_avg_quintile)",
          " + offset(lpopulation)",
          "+as.factor(urbanization_rate_quintile)",
          "+popolation_age_65_Ratio","+Male_Ratio",
          " + t2m_mean_summer + t2m_mean_winter + t2m_sd_summer + t2m_sd_winter",
          " + s(climate_zone_name_code, bs='re')"
        ))
        
        model <- tryCatch(
          mgcv::gam(fml, family = quasipoisson, data = dat_mod, method = "REML"),
          error = function(e) NULL
        )
        if (is.null(model)) next
        
        pt <- as.data.frame(summary(model)$p.table)
        if (!exposure_ %in% rownames(pt)) next
        
        beta <- unname(pt[exposure_, "Estimate"])
        se   <- unname(pt[exposure_, "Std. Error"])
        pval <- unname(pt[exposure_, "Pr(>|t|)"])
        
        # 用模型有效样本来选 delta（z_ -> 1SD；非z_按 effect_scale）
        mf <- model.frame(model)
        x  <- mf[[exposure_]]
        
        dd <- pick_delta(exposure_, x, effect_scale = effect_scale)
        if (is.null(dd)) next
        
        delta_use <- dd$delta
        eff_scale_used <- dd$scale_used
        
        # 统一用 effect_with_delta（向量版）算 RR/PC（这里 beta/se 是标量也能用）
        eff <- effect_with_delta_vec(beta = beta, se = se, delta = delta_use, z = z)
        
        out <- tibble(
          disease = dataset_name,
          strate_factor = strate_factor,
          cut_type = "median",
          cut_value = med,
          fenzu = fenzu,  # 0=Low,1=High
          group_label = ifelse(fenzu == 1, "High(>median)", "Low(<=median)"),
          exposure = exposure_,
          beta_per1 = beta,     # 每1单位（模型系数口径）
          se_per1   = se,
          p = pval,
          effect_scale = eff_scale_used,
          delta = delta_use
        ) %>% bind_cols(eff)
        
        result_long_cat_modification <- bind_rows(result_long_cat_modification, out)
      }
    }
  }
}
# write.csv(result_long_cat_modification,"E:/PKU_fanyi/PM2.5/result/result_long_modification/result_pollution_cat_modification.csv",
#           fileEncoding = 'GBK',row.names = F)
fwrite(result_long_cat_modification,"E:/Deskbook/Result_generate/1_Urban_environment/result_long_modification/result_exposure_cat_modification.csv",
          row.names = F)



####非线性修饰作用#####
mode_fun <- function(v) {
  v <- v[!is.na(v)]
  if (!length(v)) return(NA)
  ux <- unique(v)
  ux[which.max(tabulate(match(v, ux)))]
}

# ---- 核心：用 lpmatrix 计算 “对 exposure 的偏导”
# 返回 list(derivative, se)
derivative_by_lpmatrix <- function(model, newdata, var, eps) {
  nd1 <- newdata; nd2 <- newdata
  nd1[[var]] <- nd1[[var]] + eps
  nd2[[var]] <- nd2[[var]] - eps
  
  X1 <- predict(model, newdata = nd1, type = "lpmatrix")
  X2 <- predict(model, newdata = nd2, type = "lpmatrix")
  dX <- (X1 - X2) / (2 * eps)
  
  beta <- coef(model)
  V <- vcov(model)
  
  deriv <- as.numeric(dX %*% beta)
  se    <- sqrt(pmax(0, rowSums((dX %*% V) * dX)))
  list(derivative = deriv, se = se)
}

# ---- 非线性修饰：支持 exposure_vars_all + 多个 pol_ref 分位点
run_nonlinear_effectmod <- function(data, disease_name,
                                    exposure_, urban_,
                                    ns_value = 4,
                                    k_smooth = 3,
                                    k_ti = 3,
                                    n_grid = 40,
                                    min_n = 100,
                                    q_lo = 0.01,
                                    q_hi = 0.99,
                                    ref_probs = c(0.10, 0.25, 0.50, 0.75, 0.90),
                                    z = qnorm(0.975)) {
  
  need <- c("n","dead_year","total.population","GDP_avg_quintile",
            "urbanization_rate",
            "popolation_age_65_Ratio","Male_Ratio",
            "t2m_mean_summer","t2m_mean_winter","t2m_sd_summer","t2m_sd_winter",
            "climate_zone_name_code", exposure_, urban_)
  if (any(!need %in% names(data))) return(NULL)
  
  dat_cc <- data %>%
    mutate(
      dead_year = as.numeric(dead_year),
      urbanization_rate_quintile = as.factor(dplyr::ntile(urbanization_rate, 5)),
      GDP_avg_quintile = as.factor(GDP_avg_quintile),
      climate_zone_name_code = droplevels(as.factor(climate_zone_name_code)),
      lpopulation = log(total.population)
    ) %>%
    dplyr::select(all_of(need), urbanization_rate_quintile, lpopulation) %>%
    na.omit()
  
  if (nrow(dat_cc) < min_n) return(NULL)
  if (!is.numeric(dat_cc[[exposure_]]) || !is.numeric(dat_cc[[urban_]])) return(NULL)
  if (sd(dat_cc[[exposure_]], na.rm = TRUE) == 0) return(NULL)
  if (sd(dat_cc[[urban_]], na.rm = TRUE) == 0) return(NULL)
  if (nlevels(dat_cc$climate_zone_name_code) < 2) return(NULL)
  
  fml <- as.formula(paste0(
    "n ~ s(", exposure_, ", k=", k_smooth, ")",
    " + s(", urban_, ", k=", k_smooth, ")",
    " + ti(", exposure_, ", ", urban_, ", k=c(", k_ti, ",", k_ti, "))",
    " + splines::ns(dead_year, ", ns_value, ")",
    " + as.factor(GDP_avg_quintile)",
    " + offset(lpopulation)",
    #"+as.factor(urbanization_rate_quintile)",
    #"+popolation_age_65_Ratio","+Male_Ratio",
    " + t2m_mean_summer + t2m_mean_winter + t2m_sd_summer + t2m_sd_winter",
    " + s(climate_zone_name_code, bs='re')"
  ))
  
  m <- tryCatch(
    mgcv::gam(fml, family = quasipoisson, data = dat_cc, method = "REML", na.action = NULL),
    error = function(e) NULL
  )
  if (is.null(m)) return(NULL)
  
  # delta 固定 = 1（你要求“标尺用1”）
  # 若 exposure_ 是 z_ 标化变量，则 +1 = +1 SD（但 delta 仍为 1）
  delta_use <- 1
  effect_scale_used <- ifelse(grepl("^z_", exposure_), "SD(z_: +1=+1SD)", "1(unit)")
  
  # # 记录一下 IQR/SD（不用于 delta，只用于信息输出）
  ex <- dat_cc[[exposure_]]
  delta_iqr <- as.numeric(IQR(ex, na.rm = TRUE))
  delta_sd  <- as.numeric(sd(ex, na.rm = TRUE))
  if (!is.finite(delta_iqr) || delta_iqr <= 0) delta_iqr <- NA_real_
  if (!is.finite(delta_sd)  || delta_sd  <= 0) delta_sd  <- NA_real_
  
  # urban 网格（1%-99%）
  urb_vec <- dat_cc[[urban_]]
  lo <- as.numeric(quantile(urb_vec, q_lo, na.rm = TRUE))
  hi <- as.numeric(quantile(urb_vec, q_hi, na.rm = TRUE))
  if (!is.finite(lo) || !is.finite(hi) || lo == hi) return(NULL)
  urb_grid <- seq(lo, hi, length.out = n_grid)
  
  # newdata 的“模板”：先固定协变量
  newdata0 <- dat_cc[rep(1, length(urb_grid)), , drop = FALSE]
  newdata0[[urban_]] <- urb_grid
  
  newdata0$dead_year <- median(dat_cc$dead_year, na.rm = TRUE)
  newdata0$total.population <- median(dat_cc$total.population, na.rm = TRUE)
  newdata0$lpopulation <- log(newdata0$total.population)
  
  newdata$popolation_age_65_Ratio <- mean(dat_cc$popolation_age_65_Ratio, na.rm = TRUE) 
  newdata$Male_Ratio <- mean(dat_cc$Male_Ratio, na.rm = TRUE)
  
  
  newdata0$t2m_mean_summer <- mean(dat_cc$t2m_mean_summer, na.rm = TRUE)
  newdata0$t2m_mean_winter <- mean(dat_cc$t2m_mean_winter, na.rm = TRUE)
  newdata0$t2m_sd_summer   <- mean(dat_cc$t2m_sd_summer, na.rm = TRUE)
  newdata0$t2m_sd_winter   <- mean(dat_cc$t2m_sd_winter, na.rm = TRUE)
  
  newdata0$GDP_avg_quintile <- factor(mode_fun(dat_cc$GDP_avg_quintile),
                                      levels = levels(dat_cc$GDP_avg_quintile))
   newdata$urbanization_rate_quintile <-factor(mode_fun(dat_cc$urbanization_rate_quintile), 
                                levels = levels(dat_cc$urbanization_rate_quintile))
  newdata0$climate_zone_name_code <- factor(mode_fun(dat_cc$climate_zone_name_code),
                                            levels = levels(dat_cc$climate_zone_name_code))
  
  # eps：用 exposure 的尺度
  sdx <- sd(ex, na.rm = TRUE)
  eps <- if (is.finite(sdx) && sdx > 0) 1e-4 * sdx else 1e-4
  
  # ===== 多个 pol_ref 分位点一起跑并合并输出
  out_list <- list()
  
  for (pp in ref_probs) {
    
    pol_ref <- as.numeric(quantile(ex, pp, na.rm = TRUE))
    if (!is.finite(pol_ref)) next
    
    newdata <- newdata0
    newdata[[exposure_]] <- pol_ref
    
    der <- derivative_by_lpmatrix(m, newdata, var = exposure_, eps = eps)
    d_beta <- der$derivative
    d_se   <- der$se
    
    # 用 effect_with_delta_vec：这里 delta=1
    eff <- effect_with_delta_vec(d_beta, d_se, delta = delta_use, z = z)
    
    out_list[[as.character(pp)]] <- tibble::tibble(
      disease = disease_name,
      exposure = exposure_,
      urban = urban_,
      x = urb_grid,
      pol_ref = pol_ref,
      pol_ref_p = pp,
      beta_per1 = d_beta,
      se_per1 = d_se,
      effect_scale = effect_scale_used,
      delta = delta_use
    ) %>% dplyr::bind_cols(eff)
  }
  
  dplyr::bind_rows(out_list)
}

# ========== 主循环：exposure_vars_all 全部跑 
result_long_modification_nonlinear <- tibble()

for (dataset_name in names(datasets)) {
  dat0 <- datasets[[dataset_name]]
  
  for (exposure_ in exposure_vars_all) {
    if (!exposure_ %in% names(dat0)) next
    
    for (urban_ in urban_environment_vars) {
      if (!urban_ %in% names(dat0)) next
      
      tmp <- run_nonlinear_effectmod(
        data = dat0,
        disease_name = dataset_name,
        exposure_ = exposure_,
        urban_ = urban_,
        ns_value = 4,
        k_smooth = 3,
        k_ti = 3,
        n_grid = 40,
        min_n = 100,
        ref_probs = c(0.10, 0.25, 0.50, 0.75, 0.90),
        z = qnorm(0.975)
      )
      
      if (!is.null(tmp) && nrow(tmp) > 0) {
        result_long_modification_nonlinear <- dplyr::bind_rows(result_long_modification_nonlinear, tmp)
      }
    }
  }
}
# write.csv(result_long_modification_nonlinear,"E:/PKU_fanyi/PM2.5/result/result_long_modification/result_pollution_continous_modification_nonlinear.csv",
#           fileEncoding = 'GBK',row.names = F)
fwrite(result_long_modification_nonlinear,"E:/Deskbook/Result_generate/1_Urban_environment/result_long_modification/result_exposure_nonlinear_continous_modification.csv",
          row.names = F)
#中介作用#####
## 1) 参数设置
n_boot <- 1000
ns_value <- 4
# cluster bootstrap 按城市抽样
cluster_col <- "city_code"

# 面板控制开关
use_city_fe   <- TRUE   # 推荐 TRUE：城市固定效应
use_year_fe   <- FALSE  # FALSE = 用 ns(dead_year, 4)；TRUE = 年份固定效应
use_climate_re <- TRUE  # 保留你原来的 climate zone random effect
future::plan(future::multisession, workers = 6)
#根据气候区水平决定是否加随机效应
make_re_term <- function(df, col = "climate_zone_name_code",use_re=TRUE) {
  if (!use_re) return("")
  if (!col %in% names(df)) return("")
  v <- droplevels(as.factor(df[[col]]))
  if (nlevels(v) >= 2) paste0(" + s(", col, ", bs='re')") else ""
}
#按照城市抽样
bootstrap_sample <- function(data, cluster_col = NULL) {
  n <- nrow(data)
  if (is.null(cluster_col) || !(cluster_col %in% names(data))) {
    idx <- sample.int(n, n, replace = TRUE)
    return(data[idx, , drop = FALSE])
  }
  
  # cluster bootstrap: sample cities, then bind all rows within sampled cities
  clus <- unique(data[[cluster_col]])
  clus_boot <- sample(clus, length(clus), replace = TRUE)
  dplyr::bind_rows(lapply(clus_boot, function(cc) data[data[[cluster_col]] == cc, , drop = FALSE]))
}
#
# 城市固定效应项
make_city_fe_term <- function(df, city_col = "city_code", use_city_fe = TRUE) {
  if (!use_city_fe) return("")
  if (!city_col %in% names(df)) return("")
  v <- df[[city_col]]
  if (all(is.na(v))) return("")
  paste0(" + as.factor(", city_col, ")")
}
#
# 年份控制项：默认延续你原来的 ns(dead_year, 4)，也可切换 year FE
make_year_term <- function(df, year_col = "dead_year", ns_value = 4, use_year_fe = FALSE) {
  if (!year_col %in% names(df)) return("")
  if (use_year_fe) {
    return(paste0(" + as.factor(", year_col, ")"))
  } else {
    return(paste0(" + splines::ns(", year_col, ", ", ns_value, ")"))
  }
}
#
prep_model_data <- function(data, urban_, exposure_,
                            cluster_col = "city_code",
                            climate_col = "climate_zone_name_code") {
  
  base_need <- c(
    "n", "dead_year", "total.population", "GDP_avg_quintile",
    "urbanization_rate",
    "popolation_age_65_Ratio", "Male_Ratio",
    "t2m_mean_summer", "t2m_mean_winter", "t2m_sd_summer", "t2m_sd_winter",
    urban_, exposure_
  )
  
  extra_need <- c()
  if (cluster_col %in% names(data)) extra_need <- c(extra_need, cluster_col)
  if (climate_col %in% names(data)) extra_need <- c(extra_need, climate_col)
  
  need <- unique(c(base_need, extra_need))
  if (any(!need %in% names(data))) return(NULL)
  
  out <- data %>%
    mutate(
      dead_year = as.numeric(dead_year),
      GDP_avg_quintile = as.factor(GDP_avg_quintile),
      urbanization_rate_quintile = as.factor(dplyr::ntile(urbanization_rate, 5)),
      lpopulation = ifelse(is.finite(total.population) & total.population > 0,
                           log(total.population), NA_real_)
    )
  
  if (climate_col %in% names(out)) {
    out[[climate_col]] <- droplevels(as.factor(out[[climate_col]]))
  }
  if (cluster_col %in% names(out)) {
    out[[cluster_col]] <- as.factor(out[[cluster_col]])
  }
  
  out <- out %>%
    dplyr::select(all_of(need), lpopulation) %>%
    na.omit()
  
  out
}
#
#
# 统一生成 delta
get_delta_info <- function(dat_scale, urban_) {
  delta_per1   <- 1
  delta_per10  <- 10
  delta_perIQR <- as.numeric(IQR(dat_scale[[urban_]], na.rm = TRUE))
  delta_perSD  <- as.numeric(sd(dat_scale[[urban_]],  na.rm = TRUE))
  
  if (!is.finite(delta_perIQR) || delta_perIQR <= 0) delta_perIQR <- NA_real_
  if (!is.finite(delta_perSD)  || delta_perSD  <= 0) delta_perSD  <- NA_real_
  
  tibble(
    delta_per1   = delta_per1,
    delta_per10  = delta_per10,
    delta_perIQR = delta_perIQR,
    delta_perSD  = delta_perSD
  )
}

# 提取系数
safe_coef <- function(fit, term) {
  if (is.null(fit)) return(NA_real_)
  cf <- tryCatch(stats::coef(fit), error = function(e) NULL)
  if (is.null(cf)) return(NA_real_)
  if (!term %in% names(cf)) return(NA_real_)
  unname(cf[term])
}

# 结果汇总函数
summarise_boot <- function(df, cols_keep) {
  means <- df %>%
    summarise(across(all_of(cols_keep), ~ mean(.x, na.rm = TRUE)))
  
  cis <- df %>%
    summarise(across(
      all_of(cols_keep),
      list(
        l = ~ as.numeric(quantile(.x, 0.025, na.rm = TRUE)),
        u = ~ as.numeric(quantile(.x, 0.975, na.rm = TRUE))
      )
    ))
  
  list(means = means, cis = cis)
}
### 3) 单次原始样本拟合：方法1 = 近似路径分解
fit_path_approx_once <- function(data, urban_, exposure_,
                                 ns_value = 4,
                                 cluster_col = "city_code",
                                 climate_col = "climate_zone_name_code",
                                 use_city_fe = TRUE,
                                 use_year_fe = FALSE,
                                 use_climate_re = TRUE) {
  
  dat <- prep_model_data(
    data = data,
    urban_ = urban_,
    exposure_ = exposure_,
    cluster_col = cluster_col,
    climate_col = climate_col
  )
  
  if (is.null(dat) || nrow(dat) < 80) {
    return(tibble(
      beta_a = NA_real_,
      beta_b = NA_real_,
      beta_c2 = NA_real_,
      indirect = NA_real_,
      direct = NA_real_,
      total = NA_real_
    ))
  }
  
  if (!is.numeric(dat[[urban_]]) || !is.numeric(dat[[exposure_]])) {
    return(tibble(
      beta_a = NA_real_,
      beta_b = NA_real_,
      beta_c2 = NA_real_,
      indirect = NA_real_,
      direct = NA_real_,
      total = NA_real_
    ))
  }
  
  if (sd(dat[[urban_]], na.rm = TRUE) == 0 || sd(dat[[exposure_]], na.rm = TRUE) == 0) {
    return(tibble(
      beta_a = NA_real_,
      beta_b = NA_real_,
      beta_c2 = NA_real_,
      indirect = NA_real_,
      direct = NA_real_,
      total = NA_real_
    ))
  }
  
  year_term <- make_year_term(dat, "dead_year", ns_value = ns_value, use_year_fe = use_year_fe)
  city_fe_term <- make_city_fe_term(dat, city_col = cluster_col, use_city_fe = use_city_fe)
  re_term <- make_re_term(dat, col = climate_col, use_re = use_climate_re)
  
  # 1) mediator model
  f_med <- as.formula(paste0(
    exposure_, " ~ ", urban_,
    year_term,
    " + as.factor(GDP_avg_quintile)",
     " + as.factor(urbanization_rate_quintile)",
     " + popolation_age_65_Ratio + Male_Ratio",
    " + t2m_mean_summer + t2m_mean_winter + t2m_sd_summer + t2m_sd_winter",
    city_fe_term,
    re_term
  ))
  
  med_fit <- tryCatch(
    mgcv::gam(f_med, family = gaussian(), data = dat, method = "REML"),
    error = function(e) NULL
  )
  
  beta_a <- safe_coef(med_fit, urban_)
  if (!is.finite(beta_a)) {
    return(tibble(
      beta_a = NA_real_,
      beta_b = NA_real_,
      beta_c2 = NA_real_,
      indirect = NA_real_,
      direct = NA_real_,
      total = NA_real_
    ))
  }
  
  # 2) outcome model
  f_out <- as.formula(paste0(
    "n ~ ", urban_, " + ", exposure_,
    year_term,
    " + as.factor(GDP_avg_quintile)",
     " + as.factor(urbanization_rate_quintile)",
     " + popolation_age_65_Ratio + Male_Ratio",
    " + offset(lpopulation)",
    " + t2m_mean_summer + t2m_mean_winter + t2m_sd_summer + t2m_sd_winter",
    city_fe_term,
    re_term
  ))
  
  out_fit <- tryCatch(
    mgcv::gam(f_out, family = quasipoisson(), data = dat, method = "REML"),
    error = function(e) NULL
  )
  
  beta_b  <- safe_coef(out_fit, exposure_)
  beta_c2 <- safe_coef(out_fit, urban_)
  
  if (!is.finite(beta_b) || !is.finite(beta_c2)) {
    return(tibble(
      beta_a = beta_a,
      beta_b = NA_real_,
      beta_c2 = NA_real_,
      indirect = NA_real_,
      direct = NA_real_,
      total = NA_real_
    ))
  }
  
  indirect <- beta_a * beta_b
  direct   <- beta_c2
  total    <- indirect + direct
  
  tibble(
    beta_a = beta_a,
    beta_b = beta_b,
    beta_c2 = beta_c2,
    indirect = indirect,
    direct = direct,
    total = total
  )
}
###单次原始样本拟合：方法2 = difference method
fit_difference_once <- function(data, urban_, exposure_,
                                ns_value = 4,
                                cluster_col = "city_code",
                                climate_col = "climate_zone_name_code",
                                use_city_fe = TRUE,
                                use_year_fe = FALSE,
                                use_climate_re = TRUE) {
  
  dat <- prep_model_data(
    data = data,
    urban_ = urban_,
    exposure_ = exposure_,
    cluster_col = cluster_col,
    climate_col = climate_col
  )
  
  if (is.null(dat) || nrow(dat) < 80) {
    return(tibble(
      total_reduced = NA_real_,
      direct_full = NA_real_,
      indirect_diff = NA_real_
    ))
  }
  
  if (!is.numeric(dat[[urban_]]) || !is.numeric(dat[[exposure_]])) {
    return(tibble(
      total_reduced = NA_real_,
      direct_full = NA_real_,
      indirect_diff = NA_real_
    ))
  }
  
  if (sd(dat[[urban_]], na.rm = TRUE) == 0) {
    return(tibble(
      total_reduced = NA_real_,
      direct_full = NA_real_,
      indirect_diff = NA_real_
    ))
  }
  
  year_term <- make_year_term(dat, "dead_year", ns_value = ns_value, use_year_fe = use_year_fe)
  city_fe_term <- make_city_fe_term(dat, city_col = cluster_col, use_city_fe = use_city_fe)
  re_term <- make_re_term(dat, col = climate_col, use_re = use_climate_re)
  
  # reduced model：不调 exposure，urban 系数视为 total effect 近似
  f_red <- as.formula(paste0(
    "n ~ ", urban_,
    year_term,
    " + as.factor(GDP_avg_quintile)",
     " + as.factor(urbanization_rate_quintile)",
     " + popolation_age_65_Ratio + Male_Ratio",
    " + offset(lpopulation)",
    " + t2m_mean_summer + t2m_mean_winter + t2m_sd_summer + t2m_sd_winter",
    city_fe_term,
    re_term
  ))
  
  red_fit <- tryCatch(
    mgcv::gam(f_red, family = quasipoisson(), data = dat, method = "REML"),
    error = function(e) NULL
  )
  
  total_reduced <- safe_coef(red_fit, urban_)
  
  # full model：调 exposure 后 urban 系数视为 direct effect 近似
  f_full <- as.formula(paste0(
    "n ~ ", urban_, " + ", exposure_,
    year_term,
    " + as.factor(GDP_avg_quintile)",
     " + as.factor(urbanization_rate_quintile)",
     " + popolation_age_65_Ratio + Male_Ratio",
    " + offset(lpopulation)",
    " + t2m_mean_summer + t2m_mean_winter + t2m_sd_summer + t2m_sd_winter",
    city_fe_term,
    re_term
  ))
  
  full_fit <- tryCatch(
    mgcv::gam(f_full, family = quasipoisson(), data = dat, method = "REML"),
    error = function(e) NULL
  )
  
  direct_full <- safe_coef(full_fit, urban_)
  
  if (!is.finite(total_reduced) || !is.finite(direct_full)) {
    return(tibble(
      total_reduced = total_reduced,
      direct_full = direct_full,
      indirect_diff = NA_real_
    ))
  }
  
  indirect_diff <- total_reduced - direct_full
  
  tibble(
    total_reduced = total_reduced,
    direct_full = direct_full,
    indirect_diff = indirect_diff
  )
}

# 2) 组合 + 并行跑

## 5) bootstrap：方法1 = 近似路径分解

bootstrap_path_approx <- function(data, urban_, exposure_,
                                  ns_value = 4,
                                  cluster_col = "city_code",
                                  climate_col = "climate_zone_name_code",
                                  use_city_fe = TRUE,
                                  use_year_fe = FALSE,
                                  use_climate_re = TRUE) {
  
  boot_raw <- bootstrap_sample(data, cluster_col = cluster_col)
  
  fit_path_approx_once(
    data = boot_raw,
    urban_ = urban_,
    exposure_ = exposure_,
    ns_value = ns_value,
    cluster_col = cluster_col,
    climate_col = climate_col,
    use_city_fe = use_city_fe,
    use_year_fe = use_year_fe,
    use_climate_re = use_climate_re
  ) %>%
    as.list() %>%
    unlist()
}


## 6) bootstrap：方法2 = difference method

bootstrap_difference <- function(data, urban_, exposure_,
                                 ns_value = 4,
                                 cluster_col = "city_code",
                                 climate_col = "climate_zone_name_code",
                                 use_city_fe = TRUE,
                                 use_year_fe = FALSE,
                                 use_climate_re = TRUE) {
  
  boot_raw <- bootstrap_sample(data, cluster_col = cluster_col)
  
  fit_difference_once(
    data = boot_raw,
    urban_ = urban_,
    exposure_ = exposure_,
    ns_value = ns_value,
    cluster_col = cluster_col,
    climate_col = climate_col,
    use_city_fe = use_city_fe,
    use_year_fe = use_year_fe,
    use_climate_re = use_climate_re
  ) %>%
    as.list() %>%
    unlist()
}


## 7) 组合待分析 pair
exposure_vars_all <-c(
                     "z_mean_PM2.5_TAP",
                     #"z_mean_SO4","z_mean_NO3",
                     # "z_mean_NH4","z_mean_OM","z_mean_BC","score_PM25mix",
                      "z_sum_heat_intensity0__Tmean_92.5_7.5",
                      # "z_sum_heat_intensity_mean2__Tmean_92.5_7.5",
                      # "z_sum_heat_intensity_sum2__Tmean_92.5_7.5",
                      #"z_sum_PM2.5_TAP_event_intensity0",
                      "z_sum_heat_event0__Tmean_92.5_7.5",
                      "z_sum_PM2.5_TAP_event0"
                     # "composite_event_sum_heat_event0__Tmean_92.5_7.5_sum_PM2.5_TAP_event0",
                     # "composite_intensity_sum_heat_intensity0__Tmean_92.5_7.5_sum_PM2.5_TAP_event_intensity0"
                      #"composite_intensity_sum_heat_intensity_mean2__Tmean_92.5_7.5_sum_PM2.5_TAP_event_intensity0",
                      #"composite_intensity_sum_heat_intensity_sum2__Tmean_92.5_7.5_sum_PM2.5_TAP_event_intensity0"
)
all_pairs <- tidyr::expand_grid(
  disease  = names(datasets),
  urban    = c(domain_map_pca$exposure, "forest_shrub_grassland", "water_wetland"),
  exposure = exposure_vars_all
)
all_pairs <-all_pairs %>% 
  filter(disease %in% c("accident","total"))

## 8) 主并行循环

res_list <- furrr::future_map(
  seq_len(nrow(all_pairs)),
  function(i) {
    dn <- all_pairs$disease[i]
    urban_ <- all_pairs$urban[i]
    exposure_ <- all_pairs$exposure[i]
    
    message(sprintf("[Pair %d/%d] %s | U=%s | M=%s",
                    i, nrow(all_pairs), dn, urban_, exposure_))
    
    dat <- datasets[[dn]]
    if (is.null(dat)) return(NULL)
    if (!urban_ %in% names(dat) || !exposure_ %in% names(dat)) return(NULL)
    
    dat_scale <- prep_model_data(
      data = dat,
      urban_ = urban_,
      exposure_ = exposure_,
      cluster_col = cluster_col,
      climate_col = "climate_zone_name_code"
    )
    if (is.null(dat_scale) || nrow(dat_scale) < 80) return(NULL)
    
    delta_tbl <- get_delta_info(dat_scale, urban_)
    
    ## ---------- 原始样本结果 
    raw_path <- fit_path_approx_once(
      data = dat,
      urban_ = urban_,
      exposure_ = exposure_,
      ns_value = ns_value,
      cluster_col = cluster_col,
      climate_col = "climate_zone_name_code",
      use_city_fe = use_city_fe,
      use_year_fe = use_year_fe,
      use_climate_re = use_climate_re
    ) %>%
      mutate(
        disease = dn,
        urban = urban_,
        exposure = exposure_,
        method = "path_approx"
      ) %>%
      bind_cols(delta_tbl, .)
    
    raw_diff <- fit_difference_once(
      data = dat,
      urban_ = urban_,
      exposure_ = exposure_,
      ns_value = ns_value,
      cluster_col = cluster_col,
      climate_col = "climate_zone_name_code",
      use_city_fe = use_city_fe,
      use_year_fe = use_year_fe,
      use_climate_re = use_climate_re
    ) %>%
      mutate(
        disease = dn,
        urban = urban_,
        exposure = exposure_,
        method = "difference_method"
      ) %>%
      bind_cols(delta_tbl, .)
    
    raw_result <- bind_rows(raw_path, raw_diff)
    
    ## ---------- bootstrap：方法1 
    boots_path <- replicate(
      n_boot,
      bootstrap_path_approx(
        data = dat,
        urban_ = urban_,
        exposure_ = exposure_,
        ns_value = ns_value,
        cluster_col = cluster_col,
        climate_col = "climate_zone_name_code",
        use_city_fe = use_city_fe,
        use_year_fe = use_year_fe,
        use_climate_re = use_climate_re
      ),
      simplify = "matrix"
    )
    boots_path <- t(boots_path)
    
    boots_path_df <- as_tibble(boots_path) %>%
      mutate(
        disease  = dn,
        urban    = urban_,
        exposure = exposure_,
        method   = "path_approx",
        boot_id  = row_number()
      )
    
    boots_path_df <- bind_cols(
      delta_tbl[rep(1, nrow(boots_path_df)), ],
      boots_path_df
    ) %>%
      dplyr::select(
        disease, urban, exposure, method, boot_id,
        delta_per1, delta_per10, delta_perIQR, delta_perSD,
        beta_a, beta_b, beta_c2, indirect, direct, total
      ) %>%
      filter(if_all(c(beta_a, beta_b, beta_c2, indirect, direct, total), is.finite))
    
    n_used_path <- nrow(boots_path_df)
    
    boot_path_summary <- NULL
    if (n_used_path >= 20) {
      tmp_sum <- summarise_boot(
        boots_path_df,
        cols_keep = c("beta_a", "beta_b", "beta_c2", "indirect", "direct", "total")
      )
      
      boot_path_summary <- tibble(
        disease = dn,
        urban = urban_,
        exposure = exposure_,
        method = "path_approx",
        n_boot_used = n_used_path
      ) %>%
        bind_cols(delta_tbl, .) %>%
        bind_cols(
          setNames(as.list(tmp_sum$means[1, ]), paste0(names(tmp_sum$means), "_mean")),
          as.list(tmp_sum$cis[1, ])
        )
    }
    
    ## ---------- bootstrap：方法2 
    boots_diff <- replicate(
      n_boot,
      bootstrap_difference(
        data = dat,
        urban_ = urban_,
        exposure_ = exposure_,
        ns_value = ns_value,
        cluster_col = cluster_col,
        climate_col = "climate_zone_name_code",
        use_city_fe = use_city_fe,
        use_year_fe = use_year_fe,
        use_climate_re = use_climate_re
      ),
      simplify = "matrix"
    )
    boots_diff <- t(boots_diff)
    
    boots_diff_df <- as_tibble(boots_diff) %>%
      mutate(
        disease  = dn,
        urban    = urban_,
        exposure = exposure_,
        method   = "difference_method",
        boot_id  = row_number()
      )
    
    boots_diff_df <- bind_cols(
      delta_tbl[rep(1, nrow(boots_diff_df)), ],
      boots_diff_df
    ) %>%
      dplyr::select(
        disease, urban, exposure, method, boot_id,
        delta_per1, delta_per10, delta_perIQR, delta_perSD,
        total_reduced, direct_full, indirect_diff
      ) %>%
      filter(if_all(c(total_reduced, direct_full, indirect_diff), is.finite))
    
    n_used_diff <- nrow(boots_diff_df)
    
    boot_diff_summary <- NULL
    if (n_used_diff >= 20) {
      tmp_sum2 <- summarise_boot(
        boots_diff_df,
        cols_keep = c("total_reduced", "direct_full", "indirect_diff")
      )
      
      boot_diff_summary <- tibble(
        disease = dn,
        urban = urban_,
        exposure = exposure_,
        method = "difference_method",
        n_boot_used = n_used_diff
      ) %>%
        bind_cols(delta_tbl, .) %>%
        bind_cols(
          setNames(as.list(tmp_sum2$means[1, ]), paste0(names(tmp_sum2$means), "_mean")),
          as.list(tmp_sum2$cis[1, ])
        )
    }
    
    message(sprintf("  -> done | path n=%d | diff n=%d", n_used_path, n_used_diff))
    
    list(
      raw_result = raw_result,
      boot_draws_path = boots_path_df,
      boot_summary_path = boot_path_summary,
      boot_draws_diff = boots_diff_df,
      boot_summary_diff = boot_diff_summary
    )
  },
  .options = furrr::furrr_options(
    seed = TRUE,
    packages = c("dplyr","tibble","mgcv","splines","purrr","tidyr","data.table")
  )
)


## 9) 汇总结果
res_list <- purrr::compact(res_list)

result_mediation_raw <- bind_rows(map(res_list, "raw_result"))

result_mediation_draws_path <- bind_rows(map(res_list, "boot_draws_path"))
result_mediation_summary_path <- bind_rows(map(res_list, "boot_summary_path"))

result_mediation_draws_diff <- bind_rows(map(res_list, "boot_draws_diff"))
result_mediation_summary_diff <- bind_rows(map(res_list, "boot_summary_diff"))

# 合并成总表（便于你后续统一处理）
result_mediation_draws_all <- bind_rows(
  result_mediation_draws_path,
  result_mediation_draws_diff
)

result_mediation_summary_all <- bind_rows(
  result_mediation_summary_path,
  result_mediation_summary_diff
)


## 10) 保存结果


# 原始样本结果
write.csv(
  result_mediation_raw,
  file.path(out_dir, "result_mediation_raw_all_methods.csv"),
  fileEncoding = "GBK",
  row.names = FALSE
)

# 方法1：近似路径分解 bootstrap draws / summary
write.csv(
  result_mediation_draws_path,
  file.path(out_dir, "result_mediation_draws_path_approx.csv"),
  fileEncoding = "GBK",
  row.names = FALSE
)

write.csv(
  result_mediation_summary_path,
  file.path(out_dir, "result_mediation_summary_path_approx.csv"),
  fileEncoding = "GBK",
  row.names = FALSE
)

# 方法2：difference method bootstrap draws / summary
write.csv(
  result_mediation_draws_diff,
  file.path(out_dir, "result_mediation_draws_difference_method.csv"),
  fileEncoding = "GBK",
  row.names = FALSE
)

write.csv(
  result_mediation_summary_diff,
  file.path(out_dir, "result_mediation_summary_difference_method.csv"),
  fileEncoding = "GBK",
  row.names = FALSE
)

# 合并总表
write.csv(
  result_mediation_draws_all,
  file.path(out_dir, "result_mediation_draws_all_methods.csv"),
  fileEncoding = "GBK",
  row.names = FALSE
)

write.csv(
  result_mediation_summary_all,
  file.path(out_dir, "result_mediation_summary_all_methods.csv"),
  fileEncoding = "GBK",
  row.names = FALSE
)
#输出文件
diseases2_total2    <- datasets[["total"]] %>% as.data.frame()
diseases2_cvd2      <- datasets[["cvd"]]%>% as.data.frame()
diseases2_resp2    <- datasets[["resp"]]%>% as.data.frame()
diseases2_accident2 <- datasets[["accident"]]%>% as.data.frame()
#
write_rds(diseases2_total2,"E:/PKU_fanyi/PM2.5/data_city2/diseases2_total2.rds")
write_rds(diseases2_cvd2,"E:/PKU_fanyi/PM2.5/data_city2/diseases2_cvd2.rds")
write_rds(diseases2_resp2,"E:/PKU_fanyi/PM2.5/data_city2/diseases2_resp2.rds")
write_rds(diseases2_accident2,"E:/PKU_fanyi/PM2.5/data_city2/diseases2_accident2.rds")

###处理中介的表格####
result_exposure_mediation_summary <-read.csv("E:/Deskbook/Result_generate/1_Urban_environment/result_mediation/result_exposure_mediation_summary.csv")
# ---- 1) 工具函数：beta(log RR) -> RR, PC(%) 
beta_to_rr <- function(beta) exp(beta)
# （可选）若你的 beta 是“每单位增加”，需要用 IQR 修正：
beta_to_pc_with_iqr <- function(beta, iqr) 100 * (exp(beta * iqr) - 1)
names(result_exposure_mediation_summary)
result_exposure_mediation_summary <-result_exposure_mediation_summary %>% 
  filter(disease %in% c("total","accident")&
           urban %nin% c("forest_shrub_grassland","water_wetland")&
           exposure %in% c("z_mean_PM2.5_TAP","z_sum_heat_intensity0__Tmean_92.5_7.5"))
#

# ---- 2) 主整理函数 
tidy_mediation_summary <- function(df,
                                   use_iqr = T,
                                   iqr_col = "delta_perSD") {
  # use_iqr = TRUE 代表 beta 需要乘以 IQR 后再换算（beta 是 per-unit）
  # iqr_col  给出 IQR 数值所在列名（你目前这列看起来是 PC，不一定是 IQR 数值；只有在 use_iqr=TRUE 时才用到）
  
  df %>%
    mutate(
      # ---- 显著性：间接效应 beta 的95%CI是否跨0 ----
      indirect_sig = ifelse(indirect_l > 0 & indirect_u > 0, TRUE,
                            ifelse(indirect_l < 0 & indirect_u < 0, TRUE, FALSE)),
      
      # ---- RR（NDE/NIE/TE） ----
      RR_NDE_mean = beta_to_rr(direct_mean),
      RR_NDE_l    = beta_to_rr(direct_l),
      RR_NDE_u    = beta_to_rr(direct_u),
      
      RR_NIE_mean = beta_to_rr(indirect_mean),
      RR_NIE_l    = beta_to_rr(indirect_l),
      RR_NIE_u    = beta_to_rr(indirect_u),
      
      RR_TE_mean  = beta_to_rr(total_mean),
      RR_TE_l     = beta_to_rr(total_l),
      RR_TE_u     = beta_to_rr(total_u)
    ) %>%
    mutate(
      # ---- PC（% change）及95%CI ----
      PC_NDE_mean = if (!use_iqr) beta_to_pc(direct_mean) else beta_to_pc_with_iqr(direct_mean, .data[[iqr_col]]),
      PC_NDE_l    = if (!use_iqr) beta_to_pc(direct_l)    else beta_to_pc_with_iqr(direct_l,    .data[[iqr_col]]),
      PC_NDE_u    = if (!use_iqr) beta_to_pc(direct_u)    else beta_to_pc_with_iqr(direct_u,    .data[[iqr_col]]),
      
      PC_NIE_mean = if (!use_iqr) beta_to_pc(indirect_mean) else beta_to_pc_with_iqr(indirect_mean, .data[[iqr_col]]),
      PC_NIE_l    = if (!use_iqr) beta_to_pc(indirect_l)    else beta_to_pc_with_iqr(indirect_l,    .data[[iqr_col]]),
      PC_NIE_u    = if (!use_iqr) beta_to_pc(indirect_u)    else beta_to_pc_with_iqr(indirect_u,    .data[[iqr_col]]),
      
      PC_TE_mean  = if (!use_iqr) beta_to_pc(total_mean) else beta_to_pc_with_iqr(total_mean, .data[[iqr_col]]),
      PC_TE_l     = if (!use_iqr) beta_to_pc(total_l)    else beta_to_pc_with_iqr(total_l,    .data[[iqr_col]]),
      PC_TE_u     = if (!use_iqr) beta_to_pc(total_u)    else beta_to_pc_with_iqr(total_u,    .data[[iqr_col]])
    ) %>%
    mutate(
      # ---- 介导比例（RR尺度）: RR_NDE*(RR_NIE-1)/(RR_NDE*RR_NIE-1) ----
      # 这里的 l/u 用端点代入做近似（严格应基于bootstrap逐次计算再取分位数）
      med_prop_mean = (RR_NDE_mean * (RR_NIE_mean - 1)) / (RR_NDE_mean * RR_NIE_mean - 1),
      med_prop_l    = (RR_NDE_l    * (RR_NIE_l    - 1)) / (RR_NDE_l    * RR_NIE_l    - 1),
      med_prop_u    = (RR_NDE_u    * (RR_NIE_u    - 1)) / (RR_NDE_u    * RR_NIE_u    - 1),
      
      # ---- 显示用字符串：PC(95%CI) ----
      direct_PC_95CI   = sprintf("%.2f (%.2f, %.2f)", PC_NDE_mean, PC_NDE_l, PC_NDE_u),
      indirect_PC_95CI = sprintf("%.2f (%.2f, %.2f)", PC_NIE_mean, PC_NIE_l, PC_NIE_u),
      total_PC_95CI    = sprintf("%.2f (%.2f, %.2f)", PC_TE_mean,  PC_TE_l,  PC_TE_u),
      
      # 介导比例显示（通常用%更直观）
      med_prop_pct_mean = 100 * med_prop_mean,
      med_prop_pct_mean= paste0(sprintf("%.2f", med_prop_pct_mean),"%")
    ) %>%
    dplyr::select(
      disease, urban, exposure, n_boot_used, delta_perSD,
      direct_mean,direct_l,direct_u,
      indirect_mean,indirect_l,indirect_u,
      total_mean,total_l,total_u,
      PC_NDE_mean,PC_NDE_l,PC_NDE_u,
      PC_NIE_mean,PC_NIE_l,PC_NIE_u,
      PC_TE_mean,PC_TE_l,PC_TE_u,
      direct_PC_95CI, indirect_PC_95CI, total_PC_95CI,
      indirect_sig,
      med_prop_pct_mean
      
    )
}

# ---- 3) 运行 
mediation_tidy <- tidy_mediation_summary(result_exposure_mediation_summary,
                                         use_iqr = T)  # 通常你这种 summary 直接 FALSE
mediation_tidy <- mediation_tidy %>%
  mutate(
    PC_NDE_se_l = (PC_NDE_mean - PC_NDE_l) / 1.96,
    PC_NDE_se_u = (PC_NDE_u - PC_NDE_mean) / 1.96,
    
    PC_NIE_se_l = (PC_NIE_mean - PC_NIE_l) / 1.96,
    PC_NIE_se_u = (PC_NIE_u - PC_NIE_mean) / 1.96,
    
    PC_TE_se_l  = (PC_TE_mean  - PC_TE_l)  / 1.96,
    PC_TE_se_u  = (PC_TE_u     - PC_TE_mean) / 1.96
  )

write.csv(mediation_tidy,
          "E:/Deskbook/Result_generate/1_Urban_environment/result_mediation/mediation_tidy.csv",
          row.names = F)

#COV分层修饰作用#####
set.seed(123)


# 0) PATHS & SETTINGS
in_dir  <- "E:/PKU_fanyi/PM2.5/data_city"
out_dir <- "E:/PKU_fanyi/PM2.5/result/result_long_strate"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ns_value <- 4
z975 <- qnorm(0.975)

n_workers <- 2
future::plan(future::multisession, workers = n_workers)

# 你已经跑好的 exposure direct 结果（关键！不再重跑）
result_exposure_all_file <- "E:/PKU_fanyi/PM2.5/result/result_long_modification/result_exposure_direct_all.csv"
# ↑ 改成你真实那份“之前跑好的” result_exposure_all.csv / result_exposure_all_strata.csv 路径


# 0.1) REQUIRED OBJECTS CHECK
need_objs <- c("base_exposure_vars", "score_specs", "urban_environment_vars")
miss_objs <- need_objs[!vapply(need_objs, exists, logical(1), envir = .GlobalEnv)]
if (length(miss_objs) > 0) {
  stop(
    "脚本运行前请先在环境中定义/加载以下对象：\n- ",
    paste(miss_objs, collapse = "\n- "),
    "\n\n例如：source('your_define_vars.R') 或 load('xxx.RData')"
  )
}


# 1) Helper: parse disease + strata info from filename
parse_strata_info <- function(fname) {
  nm <- tools::file_path_sans_ext(basename(fname))
  rest <- sub("^diseases2_", "", nm)
  rest2 <- sub("^(total|cvd|resp|accident)\\d*", "", rest)
  rest2 <- sub("^_", "", rest2)
  
  if (rest2 == "") {
    return(tibble(strat_var = NA_character_, strat_level = NA_character_))
  }
  
  # sex1/age2
  if (grepl("^[A-Za-z]+\\d+$", rest2)) {
    return(tibble(
      strat_var = sub("\\d+$", "", rest2),
      strat_level = sub("^[A-Za-z]+", "", rest2)
    ))
  }
  # cluster_1
  if (grepl("^[A-Za-z]+_\\d+$", rest2)) {
    return(tibble(
      strat_var = sub("_\\d+$", "", rest2),
      strat_level = sub("^.*_", "", rest2)
    ))
  }
  # Climate_Alpine Plateau Climate / 其它：按第一个 "_" 切
  tibble(
    strat_var = sub("_.*$", "", rest2),
    strat_level = sub("^[^_]+_", "", rest2)
  )
}

parse_disease_info <- function(fname) {
  nm <- tools::file_path_sans_ext(basename(fname))
  rest <- sub("^diseases2_", "", nm)
  m <- stringr::str_match(rest, "^(total|cvd|resp|accident)(\\d*)")
  if (is.na(m[1,2])) {
    return(tibble(disease_base = NA_character_, disease_version = 0L))
  }
  dv <- m[1,3]
  dv <- ifelse(is.na(dv) || dv == "", "0", dv)
  tibble(disease_base = m[1,2], disease_version = as.integer(dv))
}

# 2) Read all rds files => datasets_all + meta
rds_files <- list.files(in_dir, pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(rds_files) > 0)

meta <- tibble(file = rds_files) %>%
  mutate(dataset_name = tools::file_path_sans_ext(basename(file))) %>%
  bind_cols(
    purrr::map_dfr(rds_files, parse_disease_info),
    purrr::map_dfr(rds_files, parse_strata_info)
  ) %>% filter(!is.na(strat_var))

datasets_all <- setNames(lapply(meta$file, function(f) {
  x <- readRDS(f)
  as.data.frame(x)
}), meta$dataset_name)

# 3) Utilities
safe_log <- function(x){
  ifelse(is.finite(x) & x > 0, log(x), NA_real_)
}

ensure_z_vars <- function(df, z_vars, prefix = "z_", min_n = 30){
  df <- as.data.frame(df)
  for (zv in z_vars) {
    if (zv %in% names(df)) next
    raw <- sub(paste0("^", prefix), "", zv)
    if (!raw %in% names(df)) next
    if (!is.numeric(df[[raw]])) next
    x <- df[[raw]]
    x2 <- x[is.finite(x)]
    if (length(x2) < min_n) next
    s <- sd(x2); if (!is.finite(s) || s == 0) next
    df[[zv]] <- (x - mean(x2)) / s
  }
  df
}

prep_common <- function(df) {
  df %>%
    mutate(
      dead_year = as.numeric(dead_year),
      GDP_avg_quintile = as.factor(GDP_avg_quintile),
      urbanization_rate_quintile = as.factor(ntile(urbanization_rate, 5)),
      climate_zone_name_code = droplevels(as.factor(climate_zone_name_code)),
      lpopulation = safe_log(total.population)
    )
}

effect_with_delta <- function(beta, se, delta, z = qnorm(0.975)) {
  logrr  <- beta * delta
  selog  <- se   * abs(delta)
  rr     <- exp(logrr)
  rr_l   <- exp(logrr - z * selog)
  rr_u   <- exp(logrr + z * selog)
  
  pc   <- (rr   - 1) * 100
  pc_l <- (rr_l - 1) * 100
  pc_u <- (rr_u - 1) * 100
  
  list(rr = rr, rr_l = rr_l, rr_u = rr_u,
       pc = pc, pc_l = pc_l, pc_u = pc_u)
}

effect_with_delta_vec <- function(beta, se, delta = 1, z = qnorm(0.975)) {
  b <- beta * delta
  s <- se * abs(delta)
  l <- b - z * s
  u <- b + z * s
  rr   <- exp(b); rr_l <- exp(l); rr_u <- exp(u)
  pc   <- (rr   - 1) * 100
  pc_l <- (rr_l - 1) * 100
  pc_u <- (rr_u - 1) * 100
  tibble(RR = rr, RR_l = rr_l, RR_u = rr_u,
         PC = pc, PC_l = pc_l, PC_u = pc_u)
}

pick_delta <- function(exposure_name, x, effect_scale = c("auto","1","SD","IQR")) {
  effect_scale <- match.arg(effect_scale)
  x <- x[is.finite(x)]
  if (!length(x)) return(NULL)
  
  delta_iqr <- as.numeric(IQR(x, na.rm = TRUE))
  delta_sd  <- as.numeric(sd(x, na.rm = TRUE))
  if (!is.finite(delta_iqr) || delta_iqr <= 0) delta_iqr <- NA_real_
  if (!is.finite(delta_sd)  || delta_sd  <= 0) delta_sd  <- NA_real_
  
  is_scaled <- grepl("^z_", exposure_name)
  
  if (is_scaled) {
    return(list(delta = 1, scale_used = "SD", x_IQR = delta_iqr, x_SD = delta_sd))
  }
  
  if (effect_scale == "1")  return(list(delta = 1, scale_used = "1",  x_IQR = delta_iqr, x_SD = delta_sd))
  if (effect_scale == "SD") return(if (is.na(delta_sd)) NULL else list(delta = delta_sd, scale_used = "SD", x_IQR = delta_iqr, x_SD = delta_sd))
  if (effect_scale == "IQR")return(if (is.na(delta_iqr)) NULL else list(delta = delta_iqr, scale_used = "IQR", x_IQR = delta_iqr, x_SD = delta_sd))
  
  if (!is.na(delta_iqr)) return(list(delta = delta_iqr, scale_used = "IQR", x_IQR = delta_iqr, x_SD = delta_sd))
  if (!is.na(delta_sd))  return(list(delta = delta_sd,  scale_used = "SD",  x_IQR = delta_iqr, x_SD = delta_sd))
  NULL
}

make_urban_grid <- function(x, n = 40, q_lo = 0.01, q_hi = 0.99) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NULL)
  lo <- as.numeric(quantile(x, q_lo, na.rm = TRUE))
  hi <- as.numeric(quantile(x, q_hi, na.rm = TRUE))
  if (!is.finite(lo) || !is.finite(hi) || lo == hi) return(NULL)
  if (lo >= -1e-6 && hi <= 1 + 1e-6) return(seq(0, 1, length.out = n))
  seq(lo, hi, length.out = n)
}

# 4) Read precomputed exposure-direct results (YOUR SCHEMA)
#    and extract overall weights for composite score

read_result_exposure_all <- function(path) {
  stopifnot(file.exists(path))
  
  x <- tryCatch(
    readr::read_csv(path, show_col_types = FALSE,
                    locale = readr::locale(encoding = "GB18030")),
    error = function(e) readr::read_csv(path, show_col_types = FALSE)
  ) %>% tibble::as_tibble()
  
  # ---- 标准化列名：disease -> disease_base
  if ("disease" %in% names(x) && !"disease_base" %in% names(x)) {
    x <- dplyr::rename(x, disease_base = disease)
  }
  
  # ---- 必要列检查
  need <- c("disease_base", "exposure", "beta", "se", "p")
  miss <- setdiff(need, names(x))
  if (length(miss) > 0) {
    stop("result_exposure_all 文件缺少必要列：", paste(miss, collapse = ", "))
  }
  
  # ---- 类型修正（防止读成字符）
  x <- x %>%
    mutate(
      disease_base = as.character(disease_base),
      exposure     = as.character(exposure),
      beta         = as.numeric(beta),
      se           = as.numeric(se),
      p            = as.numeric(p)
    )
  
  # ---- 若没有 q_fdr 就补一个（你现在已有 q_fdr，这里只是兜底）
  if (!"q_fdr" %in% names(x)) {
    x <- x %>%
      group_by(disease_base) %>%
      mutate(q_fdr = p.adjust(p, method = "BH")) %>%
      ungroup()
  } else {
    x <- x %>% mutate(q_fdr = as.numeric(q_fdr))
  }
  
  # ---- 给后续统一接口用：补齐缺失列（全 NA）
  if (!"dataset_name" %in% names(x)) x$dataset_name <- NA_character_
  if (!"strat_var"    %in% names(x)) x$strat_var    <- NA_character_
  if (!"strat_level"  %in% names(x)) x$strat_level  <- NA_character_
  
  x
}

extract_beta_base_overall <- function(result_exposure_all_df, base_exposure_vars) {
  
  # 只取你要用于加权构建 score 的 exposure（通常就是 base_exposure_vars）
  out <- result_exposure_all_df %>%
    filter(exposure %in% base_exposure_vars) %>%
    # 如果同一 disease+exposure 不小心有多条，优先保留 p 最小的一条
    arrange(disease_base, exposure, p) %>%
    group_by(disease_base, exposure) %>%
    slice(1) %>%
    ungroup()
  
  # 给你留足信息（后续调试/追溯方便）
  keep_cols <- intersect(
    c("disease_base","exposure","exposure_type","feature_type",
      "beta","se","p","q_fdr","ns","r2","dev_expl",
      "RR_1","RR_1_l","RR_1_u","PC_1","PC_1_l","PC_1_u"),
    names(out)
  )
  
  out %>% dplyr::select(all_of(keep_cols))
}

# ---- 读取你已有的 overall 结果 
result_exposure_all_precomputed <- read_result_exposure_all(result_exposure_all_file)

# ---- 抽取 overall 权重 
beta_base_overall <- extract_beta_base_overall(
  result_exposure_all_df = result_exposure_all_precomputed,
  base_exposure_vars = base_exposure_vars
)


# 5) Build composite scores (weights from OVERALL disease_base)
add_beta_weighted_scores_strata <- function(datasets_all, meta, beta_table_overall,
                                            score_specs,
                                            make_z = TRUE, z_prefix = "z_",
                                            min_n_z = 30,
                                            normalize_by_absw = FALSE,
                                            na_to_zero = FALSE,
                                            weight_filter = NULL) {
  
  score_specs <- score_specs[vapply(score_specs, function(v) length(v) > 0, logical(1))]
  score_names <- names(score_specs)
  
  weights_used <- vector("list", length(datasets_all))
  names(weights_used) <- names(datasets_all)
  
  for (nm in names(datasets_all)) {
    dat <- as.data.frame(datasets_all[[nm]])
    mi <- meta %>% filter(dataset_name == nm) %>% slice(1)
    db <- mi$disease_base
    
    weights_used[[nm]] <- list()
    
    for (sc in score_names) {
      vars0 <- unique(score_specs[[sc]])
      vars  <- intersect(vars0, names(dat))
      
      if (length(vars) == 0) {
        dat[[sc]] <- NA_real_
        if (make_z) dat[[paste0(z_prefix, sc)]] <- NA_real_
        weights_used[[nm]][[sc]] <- tibble(exposure = character(), beta = numeric())
        next
      }
      
      w <- beta_table_overall %>%
        filter(disease_base == db, exposure %in% vars) %>%
        dplyr::select(exposure, beta, p, q_fdr)
      
      if (!is.null(weight_filter)) {
        w <- w %>% filter(!!weight_filter)
      }
      
      if (nrow(w) == 0) {
        dat[[sc]] <- NA_real_
        if (make_z) dat[[paste0(z_prefix, sc)]] <- NA_real_
        weights_used[[nm]][[sc]] <- w %>% dplyr::select(exposure, beta)
        next
      }
      
      X <- as.matrix(dat[, w$exposure, drop = FALSE])
      if (na_to_zero) X[is.na(X)] <- 0
      
      b <- w$beta
      if (normalize_by_absw) {
        denom <- sum(abs(b), na.rm = TRUE)
        if (is.finite(denom) && denom > 0) b <- b / denom
      }
      
      s <- as.numeric(X %*% b)
      dat[[sc]] <- s
      
      if (make_z) {
        s2 <- s[is.finite(s)]
        if (length(s2) >= min_n_z && sd(s2) > 0) {
          dat[[paste0(z_prefix, sc)]] <- (s - mean(s2)) / sd(s2)
        } else {
          dat[[paste0(z_prefix, sc)]] <- NA_real_
        }
      }
      
      weights_used[[nm]][[sc]] <- w %>% dplyr::select(exposure, beta)
    }
    
    datasets_all[[nm]] <- dat
  }
  
  list(datasets_all = datasets_all, weights_used = weights_used)
}
#统一进行标化
#
ranges <- list(
  c("sum_heat_intensity0__Tmean_90_10", "mean_BC")
)
#
template_dt <- datasets_all[["diseases2_accident2_age1"]]
vars_to_std <- unique(unlist(lapply(ranges, \(rg) col_range(template_dt, rg[1], rg[2]))))
#
# ---- 标化函数：z-score（生成新列 z_*** 或覆盖原列）
standardize_one_df <- function(df, vars, prefix = "z_", overwrite = FALSE,
                               skip_binary = TRUE, min_n = 30){
  
  df <- data.table::as.data.table(df)
  
  # 只保留本 df 中存在的列
  v_ok <- intersect(vars, names(df))
  if (length(v_ok) == 0) return(df)
  
  # 只保留 numeric（兼容老 data.table：用 with=FALSE）
  is_num <- vapply(df[, v_ok, with = FALSE], is.numeric, logical(1))
  v_ok <- v_ok[is_num]
  if (length(v_ok) == 0) return(df)
  
  # 可选：跳过 0/1 变量
  if (skip_binary) {
    is_bin <- vapply(v_ok, function(v){
      x <- df[[v]]
      x <- x[is.finite(x)]
      length(x) > 0 && all(x %in% c(0,1))
    }, logical(1))
    v_ok <- v_ok[!is_bin]
    if (length(v_ok) == 0) return(df)
  }
  
  for (v in v_ok) {
    x  <- df[[v]]
    x2 <- x[is.finite(x)]
    if (length(x2) < min_n) next
    
    m <- mean(x2)
    s <- stats::sd(x2)
    if (!is.finite(s) || s == 0) next
    
    z <- (x - m) / s
    
    if (overwrite) {
      df[[v]] <- z
    } else {
      df[[paste0(prefix, v)]] <- z
    }
  }
  
  df
}

# ---- 批量处理 datasets_all
datasets_all <- purrr::imap(datasets_all, ~ standardize_one_df(.x, vars_to_std,
                                                               prefix = "z_", overwrite = FALSE,
                                                               skip_binary = TRUE))

#
scores_out <- add_beta_weighted_scores_strata(
  datasets_all = datasets_all,
  meta = meta,
  beta_table_overall = beta_base_overall,
  score_specs = score_specs,
  make_z = TRUE,
  z_prefix = "z_",
  normalize_by_absw = FALSE,
  na_to_zero = FALSE,
  weight_filter = NULL
)

datasets_all <- scores_out$datasets_all
weights_used <- scores_out$weights_used
#saveRDS(weights_used, file.path(out_dir, "weights_used_by_dataset.rds"))

new_vars_raw <- names(score_specs)
new_vars_z   <- paste0("z_", names(score_specs))
exposure_vars_all <- unique(c(base_exposure_vars, new_vars_raw))

# 6) STEP C: Exposure direct outputs (NO RE-FIT)
result_exposure_all <- result_exposure_all_precomputed



#### 7) Urban direct association (linear / median2group / nonlinear)####
#    注意：这里“标尺 = 1单位”，就是你问的“在哪改？”
#    -> 在 fit_urban_direct_linear_one() 里 pick_delta(..., effect_scale="1")
fit_urban_direct_linear_one <- function(data, urban_var, ns_value = 4, min_n = 30) {
  if (!urban_var %in% names(data)) return(NULL)
  
  required <- c(
    "n","dead_year","total.population","GDP_avg_quintile",
    "urbanization_rate",
    "popolation_age_65_Ratio","Male_Ratio",
    "t2m_mean_summer","t2m_mean_winter","t2m_sd_summer","t2m_sd_winter",
    "climate_zone_name_code", urban_var
  )
  miss <- setdiff(required, names(data))
  if (length(miss) > 0) return(NULL)
  
  dat <- data[, required, drop = FALSE] %>%
    prep_common() %>%
    filter(complete.cases(.))
  
  if (nrow(dat) < min_n) return(NULL)
  if (nlevels(dat$climate_zone_name_code) < 2) return(NULL)
  if (!is.numeric(dat[[urban_var]])) return(NULL)
  if (sd(dat[[urban_var]], na.rm = TRUE) == 0) return(NULL)
  dat<- dat %>%
    mutate(
      urbanization_rate_quintile = as.factor(dplyr::ntile(urbanization_rate, 5)))
  dat$dead_year <- as.numeric(dat$dead_year)
  dat$GDP_avg_quintile <- as.factor(dat$GDP_avg_quintile)
  dat$climate_zone_name_code <- as.factor(dat$climate_zone_name_code)
  dat$lpopulation <- log(dat$total.population)
  fml <- as.formula(paste0(
    "n ~ ", urban_var,
    " + splines::ns(dead_year, ", ns_value, ")",
    " + as.factor(GDP_avg_quintile)",
    " + offset(lpopulation)",
    #"+as.factor(urbanization_rate_quintile)",
    #"+popolation_age_65_Ratio","+Male_Ratio",
    " + t2m_mean_summer + t2m_mean_winter + t2m_sd_summer + t2m_sd_winter",
    " + s(climate_zone_name_code, bs='re')"
  ))
  
  m <- tryCatch(mgcv::gam(fml, family = quasipoisson, data = dat, method = "REML"),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)
  
  pt <- as.data.frame(summary(m)$p.table)
  if (!urban_var %in% rownames(pt)) return(NULL)
  
  beta <- unname(pt[urban_var, "Estimate"])
  se   <- unname(pt[urban_var, "Std. Error"])
  pval <- unname(pt[urban_var, "Pr(>|t|)"])
  
  mf <- model.frame(m)
  x  <- mf[[urban_var]]
  
  # ======== 你问“在哪改标尺？”就在这里 
  # effect_scale = "1" -> per 1 unit
  # effect_scale = "IQR" -> per IQR
  dd <- pick_delta(urban_var, x, effect_scale = "IQR")
  if (is.null(dd)) return(NULL)
  
  eff <- effect_with_delta(beta, se, dd$delta, z = z975)
  
  tibble(
    exposure = urban_var,
    beta = beta, se = se, p = pval,
    delta = dd$delta, effect_scale = dd$scale_used,
    RR = eff$rr, RR_l = eff$rr_l, RR_u = eff$rr_u,
    PC = eff$pc, PC_l = eff$pc_l, PC_u = eff$pc_u,
    x_IQR = dd$x_IQR, x_SD = dd$x_SD
  )
}

fit_urban_direct_median2group_one <- function(data, urban_var, ns_value = 4, min_n = 30, min_each_group = 10) {
  if (!urban_var %in% names(data)) return(NULL)
  required <- c(
    "n","dead_year","total.population","GDP_avg_quintile",
    "urbanization_rate",
    "popolation_age_65_Ratio","Male_Ratio",
    "t2m_mean_summer","t2m_mean_winter","t2m_sd_summer","t2m_sd_winter",
    "climate_zone_name_code", urban_var
  )
  miss <- setdiff(required, names(data))
  if (length(miss) > 0) return(NULL)
  
  dat <- data[, required, drop = FALSE] %>%
    prep_common() %>%
    filter(complete.cases(.))
  
  if (nrow(dat) < min_n) return(NULL)
  if (nlevels(dat$climate_zone_name_code) < 2) return(NULL)
  if (!is.numeric(dat[[urban_var]])) return(NULL)
  x <- dat[[urban_var]]
  if (sd(x, na.rm = TRUE) == 0) return(NULL)
  
  med <- as.numeric(median(x, na.rm = TRUE))
  dat$x_bin <- ifelse(x > med, 1L, 0L)
  n0 <- sum(dat$x_bin == 0, na.rm = TRUE)
  n1 <- sum(dat$x_bin == 1, na.rm = TRUE)
  if (n0 < min_each_group || n1 < min_each_group) return(NULL)
  
  dat<- dat %>%
    mutate(
      urbanization_rate_quintile = as.factor(dplyr::ntile(urbanization_rate, 5)))
  dat$dead_year <- as.numeric(dat$dead_year)
  dat$GDP_avg_quintile <- as.factor(dat$GDP_avg_quintile)
  dat$climate_zone_name_code <- as.factor(dat$climate_zone_name_code)
  dat$lpopulation <- log(dat$total.population)
  
  fml <- as.formula(paste0(
    "n ~ x_bin",
    " + splines::ns(dead_year, ", ns_value, ")",
    " + as.factor(GDP_avg_quintile)",
    " + offset(lpopulation)",
    "+as.factor(urbanization_rate_quintile)",
    "+popolation_age_65_Ratio","+Male_Ratio",
    " + t2m_mean_summer + t2m_mean_winter + t2m_sd_summer + t2m_sd_winter",
    " + s(climate_zone_name_code, bs='re')"
  ))
  
  m <- tryCatch(mgcv::gam(fml, family = quasipoisson, data = dat, method = "REML"),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)
  
  pt <- as.data.frame(summary(m)$p.table)
  if (!"x_bin" %in% rownames(pt)) return(NULL)
  
  beta <- unname(pt["x_bin", "Estimate"])
  se   <- unname(pt["x_bin", "Std. Error"])
  pval <- unname(pt["x_bin", "Pr(>|t|)"])
  
  RR   <- exp(beta)
  RR_l <- exp(beta - z975 * se)
  RR_u <- exp(beta + z975 * se)
  
  tibble(
    exposure = urban_var,
    cut_type = "median",
    cut_value = med,
    n_used = nrow(dat),
    n_low = n0,
    n_high = n1,
    beta = beta, se = se, p = pval,
    RR = RR, RR_l = RR_l, RR_u = RR_u,
    PC = (RR - 1) * 100,
    PC_l = (RR_l - 1) * 100,
    PC_u = (RR_u - 1) * 100
  )
}

mode_fun <- function(v) {
  v <- v[!is.na(v)]
  if (!length(v)) return(NA)
  ux <- unique(v)
  ux[which.max(tabulate(match(v, ux)))]
}

prep_dat_cc_urban <- function(data, urban_var) {
  required <- c(
    "n","dead_year","total.population","GDP_avg_quintile",
    "urbanization_rate",
    "popolation_age_65_Ratio","Male_Ratio",
    "t2m_mean_summer","t2m_mean_winter","t2m_sd_summer","t2m_sd_winter",
    "climate_zone_name_code", urban_var
  )
  miss <- setdiff(required, names(data))
  if (length(miss) > 0) return(NULL)
  
  dat_cc <- data %>%
    mutate(
      dead_year = as.numeric(dead_year),
      GDP_avg_quintile = as.factor(GDP_avg_quintile),
      climate_zone_name_code = droplevels(as.factor(climate_zone_name_code)),
      lpopulation = safe_log(total.population),
      urbanization_rate_quintile = as.factor(dplyr::ntile(urbanization_rate, 5))
    ) %>%
    dplyr::select(n, dead_year, total.population, lpopulation,
                  GDP_avg_quintile,
                  urbanization_rate_quintile,
                  popolation_age_65_Ratio,Male_Ratio,
                  t2m_mean_summer, t2m_mean_winter, t2m_sd_summer, t2m_sd_winter,
                  climate_zone_name_code,
                  all_of(urban_var)) %>%
    na.omit()
  
  if (nrow(dat_cc) < 30) return(NULL)
  if (nlevels(dat_cc$climate_zone_name_code) < 2) return(NULL)
  if (!is.numeric(dat_cc[[urban_var]])) return(NULL)
  if (sd(dat_cc[[urban_var]], na.rm = TRUE) == 0) return(NULL)
  dat_cc
}

fit_urban_direct_nonlinear_curve_one <- function(data, urban_var,
                                                 ns_value = 4, k_smooth = 3,
                                                 grid_n = 100, q_lo = 0.01, q_hi = 0.99) {
  dat_cc <- prep_dat_cc_urban(data, urban_var)
  if (is.null(dat_cc)) return(NULL)
  
  fml <- as.formula(paste0(
    "n ~ s(", urban_var, ", k=", k_smooth, ")",
    " + ns(dead_year, ", ns_value, ")",
    " + as.factor(GDP_avg_quintile)",
    " + offset(lpopulation)",
    "+as.factor(urbanization_rate_quintile)",
    "+popolation_age_65_Ratio","+Male_Ratio",
    " + t2m_mean_summer + t2m_mean_winter + t2m_sd_summer + t2m_sd_winter",
    " + s(climate_zone_name_code, bs='re')"
  ))
  
  m <- tryCatch(mgcv::gam(fml, family = quasipoisson, data = dat_cc, method = "REML"),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)
  
  x_vec <- dat_cc[[urban_var]]
  x_grid <- seq(quantile(x_vec, q_lo, na.rm = TRUE),
                quantile(x_vec, q_hi, na.rm = TRUE),
                length.out = grid_n)
  x_ref <- as.numeric(median(x_vec, na.rm = TRUE))
  
  newdata <- dat_cc[rep(1, grid_n), , drop = FALSE]
  newdata[[urban_var]] <- x_grid
  newdata$dead_year <- median(dat_cc$dead_year, na.rm = TRUE)
  newdata$total.population <- median(dat_cc$total.population, na.rm = TRUE)
  newdata$lpopulation <- safe_log(newdata$total.population)
  
  newdata$popolation_age_65_Ratio <- mean(dat_cc$popolation_age_65_Ratio, na.rm = TRUE)
  newdata$Male_Ratio <- mean(dat_cc$Male_Ratio, na.rm = TRUE)
  
  newdata$t2m_mean_summer <- mean(dat_cc$t2m_mean_summer, na.rm = TRUE)
  newdata$t2m_mean_winter <- mean(dat_cc$t2m_mean_winter, na.rm = TRUE)
  newdata$t2m_sd_summer   <- mean(dat_cc$t2m_sd_summer, na.rm = TRUE)
  newdata$t2m_sd_winter   <- mean(dat_cc$t2m_sd_winter, na.rm = TRUE)
  newdata$GDP_avg_quintile <- factor(mode_fun(dat_cc$GDP_avg_quintile), levels = levels(dat_cc$GDP_avg_quintile))
  #newdata$urbanization_rate_quintile <- mode_fun(dat_cc$urbanization_rate_quintile)
  
  newdata$climate_zone_name_code <- factor(mode_fun(dat_cc$climate_zone_name_code), levels = levels(dat_cc$climate_zone_name_code))
  
  pr <- predict(m, newdata = newdata, type = "terms", se.fit = TRUE)
  term_name <- grep(paste0("s\\(", urban_var, "\\)"), colnames(pr$fit), value = TRUE)
  if (!length(term_name)) return(NULL)
  term_name <- term_name[1]
  
  fit_term <- as.numeric(pr$fit[, term_name])
  se_term  <- as.numeric(pr$se.fit[, term_name])
  
  new_ref <- newdata[1, , drop = FALSE]
  new_ref[[urban_var]] <- x_ref
  pr_ref <- predict(m, newdata = new_ref, type = "terms", se.fit = TRUE)
  fit_ref <- as.numeric(pr_ref$fit[, term_name])
  
  logRR   <- fit_term - fit_ref
  logRR_l <- logRR - z975 * se_term
  logRR_u <- logRR + z975 * se_term
  
  RR   <- exp(logRR); RR_l <- exp(logRR_l); RR_u <- exp(logRR_u)
  
  tibble(
    exposure = urban_var,
    x = x_grid, x_ref = x_ref,
    logRR = logRR, RR = RR, RR_l = RR_l, RR_u = RR_u,
    PC = (RR - 1)*100, PC_l = (RR_l - 1)*100, PC_u = (RR_u - 1)*100,
    ns = ns_value, k = k_smooth, n_used = nrow(dat_cc)
  )
}

run_urban_direct_all <- function(datasets_all, meta, urban_vars, ns_value = 4) {
  out_lin <- list(); out_bin <- list(); out_nl <- list()
  
  for (nm in names(datasets_all)) {
    dat <- datasets_all[[nm]]
    mi <- meta %>% filter(dataset_name == nm) %>% slice(1)
    
    vars_use <- intersect(urban_vars, names(dat))
    if (!length(vars_use)) next
    
    lin <- bind_rows(lapply(vars_use, function(v) fit_urban_direct_linear_one(dat, v, ns_value)))
    if (nrow(lin)) {
      lin <- lin %>% mutate(dataset_name = nm, disease_base = mi$disease_base,
                            strat_var = mi$strat_var, strat_level = mi$strat_level,
                            q_fdr = p.adjust(p, method = "BH"))
      out_lin[[nm]] <- lin
    }
    
    bin <- bind_rows(lapply(vars_use, function(v) fit_urban_direct_median2group_one(dat, v, ns_value)))
    if (nrow(bin)) {
      bin <- bin %>% mutate(dataset_name = nm, disease_base = mi$disease_base,
                            strat_var = mi$strat_var, strat_level = mi$strat_level,
                            q_fdr = p.adjust(p, method = "BH"))
      out_bin[[nm]] <- bin
    }
    
    nl <- bind_rows(lapply(vars_use, function(v) fit_urban_direct_nonlinear_curve_one(dat, v, ns_value)))
    if (nrow(nl)) {
      nl <- nl %>% mutate(dataset_name = nm, disease_base = mi$disease_base,
                          strat_var = mi$strat_var, strat_level = mi$strat_level)
      out_nl[[nm]] <- nl
    }
  }
  
  list(
    linear = bind_rows(out_lin),
    median2group = bind_rows(out_bin),
    nonlinear = bind_rows(out_nl)
  )
}

urban_direct_out <- run_urban_direct_all(datasets_all, meta, urban_environment_vars, ns_value)

fwrite(urban_direct_out$linear,
       file.path(out_dir, "result_urban_direct_linear_strata.csv"),
       row.names = FALSE)
fwrite(urban_direct_out$median2group,
       file.path(out_dir, "result_urban_direct_median2group_strata.csv"),
       row.names = FALSE)
fwrite(urban_direct_out$nonlinear,
       file.path(out_dir, "result_urban_direct_nonlinear_curve_strata.csv"),
       row.names = FALSE)


#### 8) Effect modification (linear / categorical / nonlinear)####
# ---------- (A) linear modification: n ~ exposure * urban ----------
effect_scale_linear <- "1"  # 你要求“标尺用1”
#
exposure_vars_all <-c("z_mean_PM2.5_TAP","z_mean_SO4","z_mean_NO3",
                      "z_mean_NH4","z_mean_OM","z_mean_BC","score_PM25mix",
                      "z_sum_heat_intensity0__Tmean_92.5_7.5",
                      #"z_sum_heat_intensity_mean2__Tmean_92.5_7.5",
                      #"z_sum_heat_intensity_sum2__Tmean_92.5_7.5",
                      "z_sum_PM2.5_TAP_event_intensity0",
                      "z_sum_heat_event0__Tmean_92.5_7.5",
                      "z_sum_PM2.5_TAP_event0",
                      "composite_event_sum_heat_event0__Tmean_92.5_7.5_sum_PM2.5_TAP_event0",
                      "composite_intensity_sum_heat_intensity0__Tmean_92.5_7.5_sum_PM2.5_TAP_event_intensity0"
                      #"composite_intensity_sum_heat_intensity_mean2__Tmean_92.5_7.5_sum_PM2.5_TAP_event_intensity0",
                      #"composite_intensity_sum_heat_intensity_sum2__Tmean_92.5_7.5_sum_PM2.5_TAP_event_intensity0"
)
urban_environment_vars<-c(
  "score_LandCover_pca","score_Functional_Service_Density_pca",
  "score_Transport_Facility_Density_pca","score_Road_Network_Density_pca",
  "score_total_pca_allvars","score_total_pca_2stage",
  "forest_shrub_grassland","water_wetland"
)
#
result_long_continuous_modification <- tibble()

for (nm in names(datasets_all)) {
  
  data <- datasets_all[[nm]] %>%
    mutate(
      dead_year = as.numeric(dead_year),
      GDP_avg_quintile = as.factor(GDP_avg_quintile),
      urbanization_rate_quintile = as.factor(dplyr::ntile(urbanization_rate, 5)),
      climate_zone_name_code = droplevels(as.factor(climate_zone_name_code)),
      lpopulation = safe_log(total.population)
    )
  
  mi <- meta %>% filter(dataset_name == nm) %>% slice(1)
  
  for (exposure_ in exposure_vars_all) {
    if (!exposure_ %in% names(data)) next
    
    for (urban_ in urban_environment_vars) {
      if (!urban_ %in% names(data)) next
      
      fml <- as.formula(paste0(
        "n ~ ", exposure_, " * ", urban_,
        " + splines::ns(dead_year, ", ns_value, ")",
        " + as.factor(GDP_avg_quintile)",
        " + offset(lpopulation)",
        " + t2m_mean_summer + t2m_mean_winter + t2m_sd_summer + t2m_sd_winter",
        " + s(climate_zone_name_code, bs='re')"
      ))
      
      m <- tryCatch(
        mgcv::gam(fml, family = quasipoisson, data = data, na.action = na.omit, method = "REML"),
        error = function(e) NULL
      )
      if (is.null(m)) next
      
      mf <- model.frame(m)
      if (!exposure_ %in% names(mf) || !urban_ %in% names(mf)) next
      if (!is.numeric(mf[[exposure_]]) || !is.numeric(mf[[urban_]])) next
      
      pol_x <- mf[[exposure_]]
      urb_x <- mf[[urban_]]
      
      dd <- pick_delta(exposure_, pol_x, effect_scale = effect_scale_linear)
      if (is.null(dd)) next
      delta_use <- dd$delta
      
      grid <- make_urban_grid(urb_x)
      if (is.null(grid)) next
      
      cf <- coef(m); V <- vcov(m)
      main_term <- exposure_
      int1 <- paste0(exposure_, ":", urban_)
      int2 <- paste0(urban_, ":", exposure_)
      int_term <- if (int1 %in% names(cf)) int1 else if (int2 %in% names(cf)) int2 else NA_character_
      if (is.na(int_term)) next
      if (!main_term %in% names(cf)) next
      
      b1 <- unname(cf[main_term])
      b2 <- unname(cf[int_term])
      
      v11 <- V[main_term, main_term]
      v22 <- V[int_term, int_term]
      v12 <- V[main_term, int_term]
      
      beta_per1 <- b1 + b2 * grid
      var_per1  <- v11 + (grid^2) * v22 + 2 * grid * v12
      se_per1   <- sqrt(pmax(var_per1, 0))
      
      eff <- effect_with_delta_vec(beta_per1, se_per1, delta = delta_use, z = z975)
      
      out <- tibble(
        dataset_name = nm,
        disease_base = mi$disease_base,
        strat_var = mi$strat_var,
        strat_level = mi$strat_level,
        exposure = exposure_,
        urban = urban_,
        x = grid,
        beta_per1 = beta_per1,
        se_per1 = se_per1,
        effect_scale = dd$scale_used,
        delta = delta_use
      ) %>% bind_cols(eff)
      
      result_long_continuous_modification <- bind_rows(result_long_continuous_modification, out)
    }
  }
}

fwrite(result_long_continuous_modification,
       file.path(out_dir, "result_exposure_continuous_modification_strata.csv"),
       row.names = FALSE)


# ---------- (B) categorical modification: stratify by urban median; within-strata run exposure ----------
result_long_cat_modification <- tibble()

for (nm in names(datasets_all)) {
  
  data0 <- datasets_all[[nm]] %>%
    mutate(
      dead_year = as.numeric(dead_year),
      GDP_avg_quintile = as.factor(GDP_avg_quintile),
      climate_zone_name_code = droplevels(as.factor(climate_zone_name_code)),
      urbanization_rate_quintile = as.factor(dplyr::ntile(urbanization_rate, 5)),
      lpopulation = safe_log(total.population)
    )
  
  mi <- meta %>% filter(dataset_name == nm) %>% slice(1)
  
  for (strate_factor in urban_environment_vars) {
    if (!strate_factor %in% names(data0)) next
    
    med <- median(data0[[strate_factor]], na.rm = TRUE)
    if (!is.finite(med)) next
    
    data <- data0 %>% mutate(.strata = ifelse(.data[[strate_factor]] > med, 1L, 0L))
    
    for (fenzu in 0:1) {
      data1 <- data %>% filter(.strata == fenzu)
      if (nrow(data1) < 50) next
      if (nlevels(droplevels(data1$climate_zone_name_code)) < 2) next
      
      for (exposure_ in exposure_vars_all) {
        if (!exposure_ %in% names(data1)) next
        
        needed <- c(
          "n","dead_year","GDP_avg_quintile","lpopulation",
          "urbanization_rate_quintile",
          "popolation_age_65_Ratio","Male_Ratio",
          "t2m_mean_summer","t2m_mean_winter","t2m_sd_summer","t2m_sd_winter",
          "climate_zone_name_code", exposure_
        )
        needed <- intersect(needed, names(data1))
        dat_mod <- data1[, needed, drop = FALSE] %>% na.omit()
        
        if (nrow(dat_mod) < 30) next
        if (!is.numeric(dat_mod[[exposure_]])) next
        if (sd(dat_mod[[exposure_]], na.rm = TRUE) == 0) next
        if (nlevels(droplevels(dat_mod$climate_zone_name_code)) < 2) next
        
        fml <- as.formula(paste0(
          "n ~ ", exposure_,
          " + splines::ns(dead_year, ", ns_value, ")",
          " + as.factor(GDP_avg_quintile)",
          " + offset(lpopulation)",
          "+as.factor(urbanization_rate_quintile)",
          "+popolation_age_65_Ratio","+Male_Ratio",
          " + t2m_mean_summer + t2m_mean_winter + t2m_sd_summer + t2m_sd_winter",
          " + s(climate_zone_name_code, bs='re')"
        ))
        
        model <- tryCatch(
          mgcv::gam(fml, family = quasipoisson, data = dat_mod, method = "REML"),
          error = function(e) NULL
        )
        if (is.null(model)) next
        
        pt <- as.data.frame(summary(model)$p.table)
        if (!exposure_ %in% rownames(pt)) next
        
        beta <- unname(pt[exposure_, "Estimate"])
        se   <- unname(pt[exposure_, "Std. Error"])
        pval <- unname(pt[exposure_, "Pr(>|t|)"])
        
        mf <- model.frame(model)
        x  <- mf[[exposure_]]
        dd <- pick_delta(exposure_, x, effect_scale = "1")
        if (is.null(dd)) next
        
        eff <- effect_with_delta_vec(beta, se, delta = dd$delta, z = z975)
        
        out <- tibble(
          dataset_name = nm,
          disease_base = mi$disease_base,
          strat_var = mi$strat_var,
          strat_level = mi$strat_level,
          strate_factor = strate_factor,
          cut_type = "median",
          cut_value = med,
          fenzu = fenzu,
          group_label = ifelse(fenzu == 1, "High(>median)", "Low(<=median)"),
          exposure = exposure_,
          beta_per1 = beta,
          se_per1   = se,
          p = pval,
          effect_scale = dd$scale_used,
          delta = dd$delta
        ) %>% bind_cols(eff)
        
        result_long_cat_modification <- bind_rows(result_long_cat_modification, out)
      }
    }
  }
}

fwrite(result_long_cat_modification,
       file.path(out_dir, "result_exposure_cat_modification_strata.csv"),
       row.names = FALSE)


# ---------- (C) nonlinear modification: lpmatrix partial derivative + multiple exposure reference probs ----------
derivative_by_lpmatrix <- function(model, newdata, var, eps) {
  nd1 <- newdata; nd2 <- newdata
  nd1[[var]] <- nd1[[var]] + eps
  nd2[[var]] <- nd2[[var]] - eps
  X1 <- predict(model, newdata = nd1, type = "lpmatrix")
  X2 <- predict(model, newdata = nd2, type = "lpmatrix")
  dX <- (X1 - X2) / (2 * eps)
  beta <- coef(model)
  V <- vcov(model)
  deriv <- as.numeric(dX %*% beta)
  se    <- sqrt(pmax(0, rowSums((dX %*% V) * dX)))
  list(derivative = deriv, se = se)
}

run_nonlinear_effectmod <- function(data, dataset_meta,
                                    exposure_, urban_,
                                    ns_value = 4,
                                    k_smooth = 3,
                                    k_ti = 3,
                                    n_grid = 40,
                                    min_n = 100,
                                    q_lo = 0.01,
                                    q_hi = 0.99,
                                    ref_probs = c(0.10, 0.25, 0.50, 0.75, 0.90),
                                    z = qnorm(0.975)) {
  
  need <- c("n","dead_year","total.population","GDP_avg_quintile",
            "popolation_age_65_Ratio","Male_Ratio",
            "urbanization_rate",
            "t2m_mean_summer","t2m_mean_winter","t2m_sd_summer","t2m_sd_winter",
            "climate_zone_name_code", exposure_, urban_)
  if (any(!need %in% names(data))) return(NULL)
  
  dat_cc <- data %>%
    mutate(
      dead_year = as.numeric(dead_year),
      GDP_avg_quintile = as.factor(GDP_avg_quintile),
      urbanization_rate_quintile = as.factor(dplyr::ntile(urbanization_rate, 5)),
      climate_zone_name_code = droplevels(as.factor(climate_zone_name_code)),
      lpopulation = safe_log(total.population)
    ) %>%
    dplyr::select(all_of(need), lpopulation) %>%
    na.omit()
  
  if (nrow(dat_cc) < min_n) return(NULL)
  if (!is.numeric(dat_cc[[exposure_]]) || !is.numeric(dat_cc[[urban_]])) return(NULL)
  if (sd(dat_cc[[exposure_]], na.rm = TRUE) == 0) return(NULL)
  if (sd(dat_cc[[urban_]], na.rm = TRUE) == 0) return(NULL)
  if (nlevels(dat_cc$climate_zone_name_code) < 2) return(NULL)
  
  fml <- as.formula(paste0(
    "n ~ s(", exposure_, ", k=", k_smooth, ")",
    " + s(", urban_, ", k=", k_smooth, ")",
    " + ti(", exposure_, ", ", urban_, ", k=c(", k_ti, ",", k_ti, "))",
    " + splines::ns(dead_year, ", ns_value, ")",
    " + as.factor(GDP_avg_quintile)",
    #"+as.factor(urbanization_rate_quintile)",
    #"+popolation_age_65_Ratio","+Male_Ratio",
    " + offset(lpopulation)",
    " + t2m_mean_summer + t2m_mean_winter + t2m_sd_summer + t2m_sd_winter",
    " + s(climate_zone_name_code, bs='re')"
  ))
  
  m <- tryCatch(
    mgcv::gam(fml, family = quasipoisson, data = dat_cc, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(m)) return(NULL)
  
  # 你要求“标尺用1”
  delta_use <- 1
  effect_scale_used <- ifelse(grepl("^z_", exposure_), "SD(z_: +1=+1SD)", "1(unit)")
  
  ex <- dat_cc[[exposure_]]
  delta_iqr <- as.numeric(IQR(ex, na.rm = TRUE))
  delta_sd  <- as.numeric(sd(ex, na.rm = TRUE))
  if (!is.finite(delta_iqr) || delta_iqr <= 0) delta_iqr <- NA_real_
  if (!is.finite(delta_sd)  || delta_sd  <= 0) delta_sd  <- NA_real_
  
  urb_vec <- dat_cc[[urban_]]
  lo <- as.numeric(quantile(urb_vec, q_lo, na.rm = TRUE))
  hi <- as.numeric(quantile(urb_vec, q_hi, na.rm = TRUE))
  if (!is.finite(lo) || !is.finite(hi) || lo == hi) return(NULL)
  urb_grid <- seq(lo, hi, length.out = n_grid)
  
  newdata0 <- dat_cc[rep(1, length(urb_grid)), , drop = FALSE]
  newdata0[[urban_]] <- urb_grid
  newdata0$dead_year <- median(dat_cc$dead_year, na.rm = TRUE)
  newdata0$total.population <- median(dat_cc$total.population, na.rm = TRUE)
  newdata0$lpopulation <- safe_log(newdata0$total.population)
  newdata0$t2m_mean_summer <- mean(dat_cc$t2m_mean_summer, na.rm = TRUE)
  newdata0$t2m_mean_winter <- mean(dat_cc$t2m_mean_winter, na.rm = TRUE)
  newdata0$t2m_sd_summer   <- mean(dat_cc$t2m_sd_summer, na.rm = TRUE)
  newdata0$t2m_sd_winter   <- mean(dat_cc$t2m_sd_winter, na.rm = TRUE)
  
  newdata$popolation_age_65_Ratio <- mean(dat_cc$popolation_age_65_Ratio, na.rm = TRUE) 
  newdata$Male_Ratio <- mean(dat_cc$Male_Ratio, na.rm = TRUE)
  
  newdata0$GDP_avg_quintile <- factor(mode_fun(dat_cc$GDP_avg_quintile), levels = levels(dat_cc$GDP_avg_quintile))
  
   newdata$urbanization_rate_quintile <-factor(mode_fun(dat_cc$urbanization_rate_quintile), 
                                levels = levels(dat_cc$urbanization_rate_quintile))
  newdata0$climate_zone_name_code <- factor(mode_fun(dat_cc$climate_zone_name_code), levels = levels(dat_cc$climate_zone_name_code))
  
  sdx <- sd(ex, na.rm = TRUE)
  eps <- if (is.finite(sdx) && sdx > 0) 1e-4 * sdx else 1e-4
  
  out_list <- list()
  
  for (pp in ref_probs) {
    pol_ref <- as.numeric(quantile(ex, pp, na.rm = TRUE))
    if (!is.finite(pol_ref)) next
    
    newdata <- newdata0
    newdata[[exposure_]] <- pol_ref
    
    der <- derivative_by_lpmatrix(m, newdata, var = exposure_, eps = eps)
    d_beta <- der$derivative
    d_se   <- der$se
    
    eff <- effect_with_delta_vec(d_beta, d_se, delta = delta_use, z = z)
    
    out_list[[as.character(pp)]] <- tibble(
      dataset_name = dataset_meta$dataset_name,
      disease_base = dataset_meta$disease_base,
      strat_var = dataset_meta$strat_var,
      strat_level = dataset_meta$strat_level,
      exposure = exposure_,
      urban = urban_,
      x = urb_grid,
      pol_ref = pol_ref,
      pol_ref_p = pp,
      beta_per1 = d_beta,
      se_per1 = d_se,
      effect_scale = effect_scale_used,
      delta = delta_use,
      ex_IQR = delta_iqr,
      ex_SD  = delta_sd
    ) %>% bind_cols(eff)
  }
  
  bind_rows(out_list)
}

result_long_modification_nonlinear <- tibble()

for (nm in names(datasets_all)) {
  dat0 <- datasets_all[[nm]]
  mi <- meta %>% filter(dataset_name == nm) %>% slice(1)
  
  for (exposure_ in exposure_vars_all) {
    if (!exposure_ %in% names(dat0)) next
    
    for (urban_ in urban_environment_vars) {
      if (!urban_ %in% names(dat0)) next
      
      tmp <- run_nonlinear_effectmod(
        data = dat0,
        dataset_meta = mi,
        exposure_ = exposure_,
        urban_ = urban_,
        ns_value = ns_value,
        k_smooth = 3,
        k_ti = 3,
        n_grid = 40,
        min_n = 100,
        ref_probs = c(0.10, 0.25, 0.50, 0.75, 0.90),
        z = z975
      )
      
      if (!is.null(tmp) && nrow(tmp) > 0) {
        result_long_modification_nonlinear <- bind_rows(result_long_modification_nonlinear, tmp)
      }
    }
  }
}

fwrite(result_long_modification_nonlinear,
       file.path(out_dir, "result_exposure_continuous_modification_nonlinear_strata.csv"),
       row.names = FALSE)



#### 9) Mediation (urban -> exposure -> outcome), cluster bootstrap####
make_re_term <- function(df, col = "climate_zone_name_code") {
  if (!col %in% names(df)) return("")
  v <- droplevels(as.factor(df[[col]]))
  if (nlevels(v) >= 2) paste0(" + s(", col, ", bs='re')") else ""
}

bootstrap_sample <- function(data, cluster_col = NULL) {
  n <- nrow(data)
  if (is.null(cluster_col) || !(cluster_col %in% names(data))) {
    idx <- sample.int(n, n, replace = TRUE)
    return(data[idx, , drop = FALSE])
  }
  clus <- unique(data[[cluster_col]])
  clus_boot <- sample(clus, length(clus), replace = TRUE)
  dplyr::bind_rows(lapply(clus_boot, function(cc) data[data[[cluster_col]] == cc, , drop = FALSE]))
}

prep_model_data_med <- function(df, urban_, exposure_) {
  need <- c("n","dead_year","total.population","GDP_avg_quintile",
            "urbanization_rate",
            "popolation_age_65_Ratio","Male_Ratio",
            "t2m_mean_summer","t2m_mean_winter","t2m_sd_summer","t2m_sd_winter",
            "climate_zone_name_code", urban_, exposure_)
  if (any(!need %in% names(df))) return(NULL)
  
  out <- df %>%
    mutate(
      dead_year = as.numeric(dead_year),
      GDP_avg_quintile = as.factor(GDP_avg_quintile),
      urbanization_rate_quintile = as.factor(dplyr::ntile(urbanization_rate, 5)),
      climate_zone_name_code = droplevels(as.factor(climate_zone_name_code)),
      lpopulation = safe_log(total.population)
    ) %>%
    dplyr::select(all_of(need), lpopulation) %>%
    na.omit()
  out
}

bootstrap_indirect_effect <- function(data, urban_, exposure_,
                                      ns_value = 4,
                                      cluster_col = "city_code") {
  
  boot_raw  <- bootstrap_sample(data, cluster_col = cluster_col)
  boot_data <- prep_model_data_med(boot_raw, urban_, exposure_)
  
  if (is.null(boot_data) || nrow(boot_data) < 80) {
    return(c(beta_a=NA,beta_b=NA,beta_c2=NA,
             indirect=NA,direct=NA,total=NA))
  }
  
  if (!is.numeric(boot_data[[urban_]]) || 
      !is.numeric(boot_data[[exposure_]])) {
    return(c(beta_a=NA,beta_b=NA,beta_c2=NA,
             indirect=NA,direct=NA,total=NA))
  }
  
  if (sd(boot_data[[urban_]], na.rm=TRUE)==0 ||
      sd(boot_data[[exposure_]], na.rm=TRUE)==0) {
    return(c(beta_a=NA,beta_b=NA,beta_c2=NA,
             indirect=NA,direct=NA,total=NA))
  }
  
  re_term <- make_re_term(boot_data, "climate_zone_name_code")
  
  # mediator model
  f_med <- as.formula(paste0(
    exposure_, " ~ ", urban_,
    " + splines::ns(dead_year, ", ns_value, ")",
    " + as.factor(GDP_avg_quintile)",
    " + t2m_mean_summer + t2m_mean_winter + t2m_sd_summer + t2m_sd_winter",
    re_term
  ))
  
  med_fit <- tryCatch(
    mgcv::gam(f_med, family=gaussian(),
              data=boot_data, method="REML"),
    error=function(e) NULL
  )
  if (is.null(med_fit)) {
    return(c(beta_a=NA,beta_b=NA,beta_c2=NA,
             indirect=NA,direct=NA,total=NA))
  }
  
  cf_med <- coef(med_fit)
  if (!urban_ %in% names(cf_med)) {
    return(c(beta_a=NA,beta_b=NA,beta_c2=NA,
             indirect=NA,direct=NA,total=NA))
  }
  beta_a <- unname(cf_med[urban_])
  
  # outcome model
  f_out <- as.formula(paste0(
    "n ~ ", urban_, " + ", exposure_,
    " + splines::ns(dead_year, ", ns_value, ")",
    " + as.factor(GDP_avg_quintile)",
    " + offset(lpopulation)",
    " + t2m_mean_summer + t2m_mean_winter + t2m_sd_summer + t2m_sd_winter",
    re_term
  ))
  
  out_fit <- tryCatch(
    mgcv::gam(f_out, family=quasipoisson(),
              data=boot_data, method="REML"),
    error=function(e) NULL
  )
  if (is.null(out_fit)) {
    return(c(beta_a=beta_a,beta_b=NA,beta_c2=NA,
             indirect=NA,direct=NA,total=NA))
  }
  
  cf_out <- coef(out_fit)
  if (!exposure_ %in% names(cf_out) ||
      !urban_ %in% names(cf_out)) {
    return(c(beta_a=beta_a,beta_b=NA,beta_c2=NA,
             indirect=NA,direct=NA,total=NA))
  }
  
  beta_b  <- unname(cf_out[exposure_])
  beta_c2 <- unname(cf_out[urban_])
  
  indirect <- beta_a * beta_b
  direct   <- beta_c2
  total    <- indirect + direct
  
  c(beta_a=beta_a,
    beta_b=beta_b,
    beta_c2=beta_c2,
    indirect=indirect,
    direct=direct,
    total=total)
}


all_pairs <- tidyr::expand_grid(
  dataset_name = meta$dataset_name,
  urban = urban_environment_vars,
  exposure = exposure_vars_all
)
all_pairs<-all_pairs %>% 
  filter(!str_detect(dataset_name,'resp2'))
all_pairs<-all_pairs %>% 
  filter(str_detect(urban,'pca')|str_detect(urban,'forest_shrub_grassland')|str_detect(urban,'water_wetland'))
all_pairs<-all_pairs %>% 
  filter(str_detect(exposure,'92.5_7.5')|str_detect(exposure,'score_PM25mix')|str_detect(exposure,'PM2.5'))
all_pairs<-all_pairs %>% 
  filter(str_detect(exposure,'92.5_7.5')|str_detect(exposure,'score_PM25mix')|str_detect(exposure,'z_sum_PM2.5_TAP_event0')|str_detect(exposure,'z_sum_PM2.5_TAP_event_intensity0')|str_detect(exposure,'z_mean_PM2.5_TAP'))
n_boot <- 400
cluster_col <- "city_code"
urban_scale <- "per1"

res_list <- furrr::future_map(
  seq_len(nrow(all_pairs)),
  function(i) {
    
    dn        <-  all_pairs$dataset_name[i]
    urban_    <- all_pairs$urban[i]
    exposure_ <- all_pairs$exposure[i]
    
    message(sprintf("[Pair %d/%d] %s | U=%s | P=%s",
                    i, nrow(all_pairs), dn, urban_, exposure_))
    
    dat <- datasets_all[[dn]]
    if (is.null(dat)) return(NULL)
    if (!urban_ %in% names(dat) || !exposure_ %in% names(dat)) return(NULL)
    
    # -------- 计算 delta（固定常数列） 
    dat_scale <- prep_model_data_med(dat, urban_, exposure_)
    if (is.null(dat_scale) || nrow(dat_scale) < 80) return(NULL)
    
    delta_per1   <- 1
    delta_per10  <- 10
    delta_perIQR <- as.numeric(IQR(dat_scale[[urban_]], na.rm = TRUE))
    delta_perSD  <- as.numeric(sd(dat_scale[[urban_]],  na.rm = TRUE))
    
    if (!is.finite(delta_perIQR) || delta_perIQR <= 0) delta_perIQR <- NA_real_
    if (!is.finite(delta_perSD)  || delta_perSD  <= 0) delta_perSD  <- NA_real_
    
    # -------- bootstrap 
    boots <- replicate(
      n_boot,
      bootstrap_indirect_effect(
        dat, urban_, exposure_,
        ns_value = ns_value,
        cluster_col = cluster_col
      ),
      simplify = "matrix"
    )
    
    boots <- t(boots)
    
    boots_df <- tibble::as_tibble(boots) %>%
      mutate(
        disease  = dn,
        urban    = urban_,
        exposure = exposure_,
        boot_id  = dplyr::row_number(),
        delta_per1   = delta_per1,
        delta_per10  = delta_per10,
        delta_perIQR = delta_perIQR,
        delta_perSD  = delta_perSD
      )
    
    # -------- 只保留 6 个 beta 列 
    cols_keep <- c("beta_a","beta_b","beta_c2",
                   "indirect","direct","total")
    
    boots_df2 <- boots_df %>%
      dplyr::select(
        disease, urban, exposure, boot_id,
        delta_per1, delta_per10,
        delta_perIQR, delta_perSD,
        all_of(cols_keep)
      ) %>%
      dplyr::filter(if_all(all_of(cols_keep), is.finite))
    
    n_used <- nrow(boots_df2)
    if (n_used < 20) {
      message(sprintf("  -> skipped (n_boot_used=%d)", n_used))
      return(NULL)
    }
    
    # -------- summary 
    means <- boots_df2 %>%
      summarise(across(all_of(cols_keep),
                       ~ mean(.x, na.rm = TRUE)))
    
    cis <- boots_df2 %>%
      summarise(across(
        all_of(cols_keep),
        list(
          l = ~ as.numeric(quantile(.x, 0.025, na.rm = TRUE)),
          u = ~ as.numeric(quantile(.x, 0.975, na.rm = TRUE))
        )
      ))
    
    summary_df <- tibble(
      disease = dn,
      urban = urban_,
      exposure = exposure_,
      n_boot_used = n_used,
      delta_per1   = delta_per1,
      delta_per10  = delta_per10,
      delta_perIQR = delta_perIQR,
      delta_perSD  = delta_perSD
    ) %>%
      bind_cols(
        setNames(as.list(means[1, ]), paste0(names(means), "_mean")),
        as.list(cis[1, ])
      )
    
    message(sprintf("  -> done (n_boot_used=%d)", n_used))
    
    list(draws = boots_df2,
         summary = summary_df)
  },
  .options = furrr::furrr_options(
    seed = TRUE,
    packages = c("dplyr","tibble","mgcv","splines","purrr","tidyr")
  )
)

res_list <- purrr::compact(res_list)

result_mediation_draws <- bind_rows(map(res_list, "draws"))
result_mediation_summary <- bind_rows(map(res_list, "summary"))

fwrite(result_mediation_draws,
       file.path(out_dir, "result_mediation_draws_strata.csv"),
       row.names = FALSE)
fwrite(result_mediation_summary,
       file.path(out_dir, "result_mediation_summary_strata.csv"),
       row.names = FALSE)

cat("\nAll done. Outputs written to:\n", out_dir, "\n")





#场景预测####
diseases2_total2<-readRDS("G:/TMP/data_updated/diseases2_total2.rds")
diseases2_cvd2<-readRDS("E:/PKU_fanyi/PM2.5/data_city/diseases2_cvd2.rds")
diseases2_accident2<-readRDS("E:/PKU_fanyi/PM2.5/data_city/diseases2_accident2.rds")
diseases2_accident2_age1<-readRDS("E:/PKU_fanyi/PM2.5/data_city/diseases2_accident2_age1.rds")
####先生成预测值#####
setwd('G:/TMP')
make_total_future_projection_table <- function(
    total_data,
    urban_exposure,
    city_var = "city",
    province_var = "province",
    year_var = "dead_year",
    future_years = c(2030, 2050, 2070),
    projection_k = 3,
    score_lower = 0,
    score_upper = 100,
    output_csv = NULL
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(mgcv)
    library(purrr)
  })
  
  prep_projection_df <- function(dat, urban_var) {
    dat %>%
      dplyr::select(
        .city = all_of(city_var),
        .province = all_of(province_var),
        dead_year = all_of(year_var),
        value = all_of(urban_var)
      ) %>%
      dplyr::mutate(
        dead_year = as.numeric(dead_year),
        value = as.numeric(value)
      ) %>%
      dplyr::filter(!is.na(.city), !is.na(.province), !is.na(dead_year), !is.na(value))
  }
  
  project_by_city_one_var <- function(df_var,
                                      target_years,
                                      projection_k,
                                      urban_var,
                                      score_lower = 0,
                                      score_upper = 100) {
    agg <- df_var %>%
      dplyr::group_by(.city, .province, dead_year) %>%
      dplyr::summarise(
        value_mean = mean(value, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(.city, dead_year)
    
    agg$.city_key <- paste(agg$.province, agg$.city, sep = "||")
    by_city <- split(agg, agg$.city_key)
    
    pred_list <- list()
    fit_list <- list()
    
    for (kk in names(by_city)) {
      dd <- by_city[[kk]]
      cc <- dd$.city[1]
      prov_cc <- dd$.province[1]
      
      yrs <- as.numeric(dd$dead_year)
      vals <- as.numeric(dd$value_mean)
      u <- length(unique(yrs))
      
      use_gam <- isTRUE(u >= 5)
      
      if (use_gam) {
        k_use <- max(3L, min(projection_k, u - 1L))
        fit <- try(
          mgcv::gam(vals ~ s(yrs, k = k_use),
                    optimizer = c("outer", "bfgs")),
          silent = TRUE
        )
        if (inherits(fit, "try-error")) use_gam <- FALSE
      }
      
      if (u < 2) {
        method_use <- "constant"
        fitted_vals <- rep(mean(vals, na.rm = TRUE), length(vals))
        pred_raw <- rep(mean(vals, na.rm = TRUE), length(target_years))
      } else if (!use_gam) {
        fit <- lm(vals ~ yrs)
        method_use <- "linear"
        fitted_vals <- fitted(fit)
        pred_raw <- predict(fit, newdata = data.frame(yrs = as.numeric(target_years)))
      } else {
        method_use <- "nonlinear"
        fitted_vals <- predict(fit, newdata = data.frame(yrs = yrs), type = "response")
        pred_raw <- predict(fit, newdata = data.frame(yrs = as.numeric(target_years)), type = "response")
      }
      
      pred_clip <- pmin(pmax(as.numeric(pred_raw), score_lower), score_upper)
      
      r2 <- tryCatch(
        cor(vals, fitted_vals, use = "complete.obs")^2,
        error = function(e) NA_real_
      )
      
      pred_list[[kk]] <- data.frame(
        disease_projection_source = "total",
        urban = urban_var,
        .city = cc,
        city = cc,
        .province = prov_cc,
        province = prov_cc,
        year = as.numeric(target_years),
        value_pred_raw = as.numeric(pred_raw),
        value_pred = pred_clip,
        value_clipped = as.logical(abs(pred_clip - as.numeric(pred_raw)) > 1e-10),
        score_lower = score_lower,
        score_upper = score_upper,
        observed_min = min(vals, na.rm = TRUE),
        observed_max = max(vals, na.rm = TRUE),
        method = method_use,
        stringsAsFactors = FALSE
      )
      
      fit_list[[kk]] <- data.frame(
        disease_projection_source = "total",
        urban = urban_var,
        .city = cc,
        city = cc,
        .province = prov_cc,
        province = prov_cc,
        method = method_use,
        n_obs = nrow(dd),
        observed_min = min(vals, na.rm = TRUE),
        observed_max = max(vals, na.rm = TRUE),
        r2 = r2,
        stringsAsFactors = FALSE
      )
    }
    
    list(
      pred = dplyr::bind_rows(pred_list),
      fit = dplyr::bind_rows(fit_list)
    )
  }
  
  vars_use <- intersect(urban_exposure, names(total_data))
  if (length(vars_use) == 0) stop("urban_exposure 在 total_data 中一个也没找到")
  
  pred_all <- list()
  fit_all <- list()
  
  for (uv in vars_use) {
    df_var <- prep_projection_df(total_data, uv)
    
    res_uv <- project_by_city_one_var(
      df_var = df_var,
      target_years = future_years,
      projection_k = projection_k,
      urban_var = uv,
      score_lower = score_lower,
      score_upper = score_upper
    )
    
    pred_all[[uv]] <- res_uv$pred
    fit_all[[uv]] <- res_uv$fit
  }
  
  pred_df <- bind_rows(pred_all)
  fit_df <- bind_rows(fit_all)
  
  out <- list(
    future_projection_values = pred_df,
    future_projection_fit = fit_df
  )
  
  if (!is.null(output_csv)) {
    write.csv(pred_df, output_csv, row.names = FALSE, fileEncoding = "UTF-8")
    fit_csv <- sub("\\.csv$", "_fit.csv", output_csv)
    write.csv(fit_df, fit_csv, row.names = FALSE, fileEncoding = "UTF-8")
  }
  
  invisible(out)
}
#
exposure_vars_all <- c(
  "z_mean_PM2.5_TAP",
  #"z_mean_SO4", "z_mean_NO3",
  #"z_mean_NH4", "z_mean_OM", "z_mean_BC", "score_PM25mix",
  "z_sum_heat_intensity0__Tmean_92.5_7.5"
  #"z_sum_PM2.5_TAP_event_intensity0",
  #"z_sum_heat_event0__Tmean_92.5_7.5",
 # "z_sum_PM2.5_TAP_event0",
 # "composite_event_sum_heat_event0__Tmean_92.5_7.5_sum_PM2.5_TAP_event0",
 # "composite_intensity_sum_heat_intensity0__Tmean_92.5_7.5_sum_PM2.5_TAP_event_intensity0"
)

urban_environment_vars <- c(
  #"forest_shrub_grassland", "water_wetland", "impervious",
  "score_LandCover_pca", "score_Functional_Service_Density_pca", "score_Transport_Facility_Density_pca",
  "score_Road_Network_Density_pca", "score_total_pca_allvars"
  #"Residential_per_totalArea", "Workplace_per_totalArea", "Industrial_production_per_totalArea",
  #"Education_culture_per_totalArea", "Healthcare_per_totalArea", "Government_per_totalArea",
  #"Commerce_services_per_totalArea", 
  #"Leisure_tourism_per_totalArea", "Logistics_delivery_per_totalArea",
  #"Urban_public_transit_per_totalArea",
 # "Intercity_passenger_hub_per_totalArea", "Motorized_transport_support_per_totalArea",
 # "total_len_km_per_totalArea", "major_len_km_per_totalArea", "minor_len_km_per_totalArea"
)
future_projection_table_total <- make_total_future_projection_table(
  total_data = diseases2_total2,
  urban_exposure = urban_environment_vars,
  city_var = "city",
  province_var = "province",
  year_var = "dead_year",
  future_years = c(2030, 2050, 2070),
  projection_k = 3
)
future_projection_table<-future_projection_table_total$future_projection_values
future_projection_fit<-future_projection_table_total$future_projection_fit
#
table(future_projection_table$value_clipped)

future_projection_table %>%
  dplyr::filter(value_clipped) %>%
  dplyr::count(urban, year, method)
#
# 保存生成的 projection 表
write.csv(future_projection_table, "future_projection_table.csv", row.names = FALSE,fileEncoding = 'UTF-8')
write.csv(future_projection_table, 
          "E:/Deskbook/Result_generate/1_Urban_environment/result_description/future_projection_table.csv", 
          row.names = FALSE,fileEncoding = 'UTF-8')

write.csv(future_projection_fit, "future_projection_fit.csv", row.names = FALSE,fileEncoding = 'UTF-8')
####针对文件夹生成z_标化变量和复合暴露变量####
# 0) PATHS & SETTINGS
in_dir  <- "E:/PKU_fanyi/PM2.5/data_city"
out_dir <- "E:/PKU_fanyi/PM2.5/data_city2"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n_workers <- 2
future::plan(future::multisession, workers = n_workers)

# 你已经跑好的 exposure direct 结果（关键！不再重跑）
result_exposure_all_file <- "E:/PKU_fanyi/PM2.5/result/result_long_modification/result_exposure_direct_all.csv"



# 1) Helper: parse disease + strata info from filename
parse_strata_info <- function(fname) {
  nm <- tools::file_path_sans_ext(basename(fname))
  rest <- sub("^diseases2_", "", nm)
  rest2 <- sub("^(total|cvd|resp|accident)\\d*", "", rest)
  rest2 <- sub("^_", "", rest2)
  
  if (rest2 == "") {
    return(tibble(strat_var = NA_character_, strat_level = NA_character_))
  }
  
  # sex1/age2
  if (grepl("^[A-Za-z]+\\d+$", rest2)) {
    return(tibble(
      strat_var = sub("\\d+$", "", rest2),
      strat_level = sub("^[A-Za-z]+", "", rest2)
    ))
  }
  # cluster_1
  if (grepl("^[A-Za-z]+_\\d+$", rest2)) {
    return(tibble(
      strat_var = sub("_\\d+$", "", rest2),
      strat_level = sub("^.*_", "", rest2)
    ))
  }
  # Climate_Alpine Plateau Climate / 其它：按第一个 "_" 切
  tibble(
    strat_var = sub("_.*$", "", rest2),
    strat_level = sub("^[^_]+_", "", rest2)
  )
}

parse_disease_info <- function(fname) {
  nm <- tools::file_path_sans_ext(basename(fname))
  rest <- sub("^diseases2_", "", nm)
  m <- stringr::str_match(rest, "^(total|cvd|resp|accident)(\\d*)")
  if (is.na(m[1,2])) {
    return(tibble(disease_base = NA_character_, disease_version = 0L))
  }
  dv <- m[1,3]
  dv <- ifelse(is.na(dv) || dv == "", "0", dv)
  tibble(disease_base = m[1,2], disease_version = as.integer(dv))
}

# 2) Read all rds files => datasets_all + meta
rds_files <- list.files(in_dir, pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(rds_files) > 0)

meta <- tibble(file = rds_files) %>%
  mutate(dataset_name = tools::file_path_sans_ext(basename(file))) %>%
  bind_cols(
    purrr::map_dfr(rds_files, parse_disease_info),
    purrr::map_dfr(rds_files, parse_strata_info)
  ) 
meta<-meta %>% filter(disease_base!="resp")
datasets_all <- setNames(lapply(meta$file, function(f) {
  x <- readRDS(f)
  as.data.frame(x)
}), meta$dataset_name)

# 4) Read precomputed exposure-direct results (YOUR SCHEMA)
#    and extract overall weights for composite score

read_result_exposure_all <- function(path) {
  stopifnot(file.exists(path))
  
  x <- tryCatch(
    readr::read_csv(path, show_col_types = FALSE,
                    locale = readr::locale(encoding = "GB18030")),
    error = function(e) readr::read_csv(path, show_col_types = FALSE)
  ) %>% tibble::as_tibble()
  
  # ---- 标准化列名：disease -> disease_base
  if ("disease" %in% names(x) && !"disease_base" %in% names(x)) {
    x <- dplyr::rename(x, disease_base = disease)
  }
  
  # ---- 必要列检查
  need <- c("disease_base", "exposure", "beta", "se", "p")
  miss <- setdiff(need, names(x))
  if (length(miss) > 0) {
    stop("result_exposure_all 文件缺少必要列：", paste(miss, collapse = ", "))
  }
  
  # ---- 类型修正（防止读成字符）
  x <- x %>%
    mutate(
      disease_base = as.character(disease_base),
      exposure     = as.character(exposure),
      beta         = as.numeric(beta),
      se           = as.numeric(se),
      p            = as.numeric(p)
    )
  
  # ---- 若没有 q_fdr 就补一个（你现在已有 q_fdr，这里只是兜底）
  if (!"q_fdr" %in% names(x)) {
    x <- x %>%
      group_by(disease_base) %>%
      mutate(q_fdr = p.adjust(p, method = "BH")) %>%
      ungroup()
  } else {
    x <- x %>% mutate(q_fdr = as.numeric(q_fdr))
  }
  
  # ---- 给后续统一接口用：补齐缺失列（全 NA）
  if (!"dataset_name" %in% names(x)) x$dataset_name <- NA_character_
  if (!"strat_var"    %in% names(x)) x$strat_var    <- NA_character_
  if (!"strat_level"  %in% names(x)) x$strat_level  <- NA_character_
  
  x
}

extract_beta_base_overall <- function(result_exposure_all_df, base_exposure_vars) {
  
  # 只取你要用于加权构建 score 的 exposure（通常就是 base_exposure_vars）
  out <- result_exposure_all_df %>%
    filter(exposure %in% base_exposure_vars) %>%
    # 如果同一 disease+exposure 不小心有多条，优先保留 p 最小的一条
    arrange(disease_base, exposure, p) %>%
    group_by(disease_base, exposure) %>%
    slice(1) %>%
    ungroup()
  
  # 给你留足信息（后续调试/追溯方便）
  keep_cols <- intersect(
    c("disease_base","exposure","exposure_type","feature_type",
      "beta","se","p","q_fdr","ns","r2","dev_expl",
      "RR_1","RR_1_l","RR_1_u","PC_1","PC_1_l","PC_1_u"),
    names(out)
  )
  
  out %>% dplyr::select(all_of(keep_cols))
}

# ---- 读取你已有的 overall 结果 
result_exposure_all_precomputed <- read_result_exposure_all(result_exposure_all_file)

# ---- 抽取 overall 权重 
beta_base_overall <- extract_beta_base_overall(
  result_exposure_all_df = result_exposure_all_precomputed,
  base_exposure_vars = base_exposure_vars
)
#统一进行标化
#
ranges <- list(
  c("sum_heat_event0__Tmean_90_10", "mean_BC")
)
#
template_dt <- datasets_all[["diseases2_accident2_age1"]]
vars_to_std <- unique(unlist(lapply(ranges, \(rg) col_range(template_dt, rg[1], rg[2]))))
# ---- 标化函数：z-score（生成新列 z_*** 或覆盖原列）
standardize_one_df <- function(df, vars, prefix = "z_", overwrite = FALSE,
                               skip_binary = TRUE, min_n = 30){
  
  df <- data.table::as.data.table(df)
  
  # 只保留本 df 中存在的列
  v_ok <- intersect(vars, names(df))
  if (length(v_ok) == 0) return(df)
  
  # 只保留 numeric（兼容老 data.table：用 with=FALSE）
  is_num <- vapply(df[, v_ok, with = FALSE], is.numeric, logical(1))
  v_ok <- v_ok[is_num]
  if (length(v_ok) == 0) return(df)
  
  # 可选：跳过 0/1 变量
  if (skip_binary) {
    is_bin <- vapply(v_ok, function(v){
      x <- df[[v]]
      x <- x[is.finite(x)]
      length(x) > 0 && all(x %in% c(0,1))
    }, logical(1))
    v_ok <- v_ok[!is_bin]
    if (length(v_ok) == 0) return(df)
  }
  
  for (v in v_ok) {
    x  <- df[[v]]
    x2 <- x[is.finite(x)]
    if (length(x2) < min_n) next
    
    m <- mean(x2)
    s <- stats::sd(x2)
    if (!is.finite(s) || s == 0) next
    
    z <- (x - m) / s
    
    if (overwrite) {
      df[[v]] <- z
    } else {
      df[[paste0(prefix, v)]] <- z
    }
  }
  
  df
}

# ---- 批量处理 datasets_all
datasets_all <- purrr::imap(datasets_all, ~ standardize_one_df(.x, vars_to_std,
                                                               prefix = "z_", overwrite = FALSE,
                                                               skip_binary = TRUE))

# 5) Build composite scores (weights from OVERALL disease_base)
add_beta_weighted_scores_strata <- function(datasets_all, meta, beta_table_overall,
                                            score_specs,
                                            make_z = TRUE, z_prefix = "z_",
                                            min_n_z = 30,
                                            normalize_by_absw = FALSE,
                                            na_to_zero = FALSE,
                                            weight_filter = NULL) {
  
  score_specs <- score_specs[vapply(score_specs, function(v) length(v) > 0, logical(1))]
  score_names <- names(score_specs)
  
  weights_used <- vector("list", length(datasets_all))
  names(weights_used) <- names(datasets_all)
  
  for (nm in names(datasets_all)) {
    dat <- as.data.frame(datasets_all[[nm]])
    mi <- meta %>% filter(dataset_name == nm) %>% slice(1)
    db <- mi$disease_base
    
    weights_used[[nm]] <- list()
    
    for (sc in score_names) {
      vars0 <- unique(score_specs[[sc]])
      vars  <- intersect(vars0, names(dat))
      
      if (length(vars) == 0) {
        dat[[sc]] <- NA_real_
        if (make_z) dat[[paste0(z_prefix, sc)]] <- NA_real_
        weights_used[[nm]][[sc]] <- tibble(exposure = character(), beta = numeric())
        next
      }
      
      w <- beta_table_overall %>%
        filter(disease_base == db, exposure %in% vars) %>%
        dplyr::select(exposure, beta, p, q_fdr)
      
      if (!is.null(weight_filter)) {
        w <- w %>% filter(!!weight_filter)
      }
      
      if (nrow(w) == 0) {
        dat[[sc]] <- NA_real_
        if (make_z) dat[[paste0(z_prefix, sc)]] <- NA_real_
        weights_used[[nm]][[sc]] <- w %>% dplyr::select(exposure, beta)
        next
      }
      
      X <- as.matrix(dat[, w$exposure, drop = FALSE])
      if (na_to_zero) X[is.na(X)] <- 0
      
      b <- w$beta
      if (normalize_by_absw) {
        denom <- sum(abs(b), na.rm = TRUE)
        if (is.finite(denom) && denom > 0) b <- b / denom
      }
      
      s <- as.numeric(X %*% b)
      dat[[sc]] <- s
      
      if (make_z) {
        s2 <- s[is.finite(s)]
        if (length(s2) >= min_n_z && sd(s2) > 0) {
          dat[[paste0(z_prefix, sc)]] <- (s - mean(s2)) / sd(s2)
        } else {
          dat[[paste0(z_prefix, sc)]] <- NA_real_
        }
      }
      
      weights_used[[nm]][[sc]] <- w %>% dplyr::select(exposure, beta)
    }
    
    datasets_all[[nm]] <- dat
  }
  
  list(datasets_all = datasets_all, weights_used = weights_used)
}

#

#
scores_out <- add_beta_weighted_scores_strata(
  datasets_all = datasets_all,
  meta = meta,
  beta_table_overall = beta_base_overall,
  score_specs = score_specs,
  make_z = TRUE,
  z_prefix = "z_",
  normalize_by_absw = FALSE,
  na_to_zero = FALSE,
  weight_filter = NULL
)

datasets_all <- scores_out$datasets_all
weights_used <- scores_out$weights_used
#saveRDS(weights_used, file.path(out_dir, "weights_used_by_dataset.rds"))

for (i in seq_len(nrow(meta))){
  nm<-meta$dataset_name[i]
  out_file<-file.path(out_dir,paste0(nm,".rds"))
  saveRDS(datasets_all[[nm]],out_file)
}




#分层解析代码
parse_strata_info <- function(fname) {
  nm <- tools::file_path_sans_ext(basename(fname))
  rest <- sub("^diseases2_", "", nm)
  rest2 <- sub("^(total|cvd|resp|accident)\\d*", "", rest)
  rest2 <- sub("^_", "", rest2)
  
  if (rest2 == "") {
    return(tibble::tibble(strat_var = NA_character_, strat_level = NA_character_))
  }
  
  # sex1 / age2
  if (grepl("^[A-Za-z]+\\d+$", rest2)) {
    return(tibble::tibble(
      strat_var = sub("\\d+$", "", rest2),
      strat_level = sub("^[A-Za-z]+", "", rest2)
    ))
  }
  
  # cluster_1
  if (grepl("^[A-Za-z]+_\\d+$", rest2)) {
    return(tibble::tibble(
      strat_var = sub("_\\d+$", "", rest2),
      strat_level = sub("^.*_", "", rest2)
    ))
  }
  
  # Climate_Alpine Plateau Climate / 其它
  tibble::tibble(
    strat_var = sub("_.*$", "", rest2),
    strat_level = sub("^[^_]+_", "", rest2)
  )
}

parse_disease_info <- function(fname) {
  nm <- tools::file_path_sans_ext(basename(fname))
  rest <- sub("^diseases2_", "", nm)
  m <- stringr::str_match(rest, "^(total|cvd|resp|accident)(\\d*)")
  if (is.na(m[1,2])) {
    return(tibble::tibble(disease_base = NA_character_, disease_version = 0L))
  }
  dv <- m[1,3]
  dv <- ifelse(is.na(dv) || dv == "", "0", dv)
  tibble::tibble(
    disease_base = m[1,2],
    disease_version = as.integer(dv)
  )
}

#现在跑预测v2####
# ============================================================
# External running script for urban resilience scenarios v4
# Parallel + stratification parsing
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(data.table)
  library(mgcv)
  library(randomForest)
  library(future)
  library(future.apply)
  library(stringr)
  library(tibble)
})

source("G:/TMP/code_document/run_urban_resilience_scenarios_final_v4.R")

base_dir <- "G:/TMP/code_document/"

long_dir <- file.path(base_dir, "data_updated")
short_dir <- file.path(base_dir, "data_short")

output_long <- file.path(base_dir, "scenario_results_long_v4")
output_short <- file.path(base_dir, "scenario_results_short_v4")

N_WORKERS <- 4

exposure_vars_long <- c(
  "z_mean_PM2.5_TAP",
  "z_sum_heat_intensity0__Tmean_92.5_7.5"
)

exposure_vars_short <- c(
  "z_lag01_PM2.5_TAP",
  "z_heat_intensity_sum2__Tmean_92.5_7.5"
)

urban_exposure <- c(
  "score_LandCover_pca",
  "score_Functional_Service_Density_pca",
  "score_Transport_Facility_Density_pca",
  "score_Road_Network_Density_pca",
  "score_total_pca_allvars"
)

nl_long <- data.table::fread(
  file.path(base_dir, "result_nonlinear_continuous_modification.csv"),
  encoding = "UTF-8"
) %>%
  as.data.frame() %>%
  mutate(
    disease = trimws(as.character(disease)),
    exposure = if ("exposure" %in% names(.)) trimws(as.character(exposure)) else exposure,
    urban = trimws(as.character(urban)),
    disease = dplyr::case_when(
      disease == "total" ~ "total2",
      disease == "accident" ~ "accident2",
      TRUE ~ disease
    )
  )

nl_short <- data.table::fread(
  file.path(base_dir, "result_short_nonlinear_continuous_modification.csv"),
  encoding = "UTF-8"
) %>%
  as.data.frame() %>%
  mutate(
    disease = if ("disease" %in% names(.)) trimws(as.character(disease)) else NA_character_,
    urban = trimws(as.character(urban))
  )

future_projection_table <- read.csv(
  file.path(base_dir, "future_projection_table.csv"),
  fileEncoding = "UTF-8"
)

pca_projection_params <- readr::read_rds(
  file.path(base_dir, "pca_projection_params.rds")
)

# ============================================================
# 先检查长期文件识别
# ============================================================

file_info_long_check <- make_file_info_v4(
  rds_dir = long_dir,
  time_scale = "long",
  include_unstratified = TRUE,
  include_stratified = TRUE
)

print(file_info_long_check %>%
        select(source_file, disease, disease_base, disease_suffix, strat_var, strat_level, is_stratified))

# write.csv(
#   file_info_long_check,
#   file.path(output_long, "CHECK_file_info_long.csv"),
#   row.names = FALSE,
#   fileEncoding = "UTF-8"
# )

# ============================================================
# Long-term scenarios
# ============================================================

res_long <- run_af_from_rds_folder_v4(
  rds_dir = long_dir,
  time_scale = "long",
  exposure_vars = exposure_vars_long,
  urban_exposure = urban_exposure,
  result_nonlinear_continuous_modification = nl_long,
  future_projection_table = future_projection_table,
  pca_projection_params = pca_projection_params,
  output_dir = output_long,
  
  include_unstratified = TRUE,
  include_stratified = TRUE,
  
  run_future = TRUE,
  run_future_short = FALSE,
  run_threshold = TRUE,
  run_policy_integrated_be = TRUE,
  
  policy_enhancement_delta = 0.30,
  
  run_compact_city = TRUE,
  run_complete_community = TRUE,
  run_active_mobility = TRUE,
  
  n_boot = 300,
  ntree = 300,
  seed = 2026,
  n_workers = N_WORKERS,
  
  death_reference_type = "avg_annual_observed",
  base_year = 2020
)
# ============================================================
# 如果是list
# ============================================================

# ============================================================
# 先检查短期文件识别
# ============================================================

file_info_short_check <- make_file_info_v4(
  rds_dir = short_dir,
  time_scale = "short",
  include_unstratified = TRUE,
  include_stratified = TRUE
)

print(file_info_short_check %>%
        select(source_file, disease, disease_base, disease_suffix, strat_var, strat_level, is_stratified))

# write.csv(
#   file_info_short_check,
#   file.path(output_short, "CHECK_file_info_short.csv"),
#   row.names = FALSE,
#   fileEncoding = "UTF-8"
# )

# ============================================================
# Short-term scenarios
# ============================================================

res_short <- run_af_from_rds_folder_v4(
  rds_dir = short_dir,
  time_scale = "short",
  file_pattern = "^diseases2_total_short\\.rds$",
  exposure_vars = exposure_vars_short,
  urban_exposure = urban_exposure,
  result_nonlinear_continuous_modification = nl_short,
  future_projection_table = NULL,
  pca_projection_params = pca_projection_params,
  output_dir = output_short,
  
  include_unstratified = TRUE,
  include_stratified = FALSE,
  
  run_future = FALSE,
  run_future_short = FALSE,
  run_threshold = TRUE,
  run_policy_integrated_be = TRUE,
  
  policy_enhancement_delta = 0.30,
  
  run_compact_city = TRUE,
  run_complete_community = TRUE,
  run_active_mobility = TRUE,
  
  n_boot = 1000,
  seed = 2026,
  n_workers = N_WORKERS,
  
  death_reference_type = "avg_annual_observed"
)

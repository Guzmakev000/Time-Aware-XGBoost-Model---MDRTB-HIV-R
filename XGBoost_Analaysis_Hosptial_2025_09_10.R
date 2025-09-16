################################################################################
############ Latent Class Growth Analysis of BDQ and ART Data ##################
################################################################################

## ---------------------------------------------------------------------------
## Utilities & Libraries
## ---------------------------------------------------------------------------
table_na <- function(..., useNA = "ifany") base::table(..., useNA = useNA)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(xgboost); library(pROC); library(caret); library(Metrics)
  library(PRROC);  library(boot);  library(Matrix); library(betareg)
  library(writexl)
})

set.seed(123)

## ---------------------------------------------------------------------------
## Load & Merge (assumes Files/*.RData paths exist and contain mentioned objects)
## ---------------------------------------------------------------------------
load("Files/BDQ_Week.RData")
load("Files/PRX3.RData")
load("Files/PRX_Aim_Hosp.RData")
load("Files/PRX_Admission.RData")
load("Files/PRX_BDQStart.RData")

BDQ_Week <- BDQ_Week %>%
  left_join(PRX3 %>% select(
    ID, Age, Gender, BMI, Education, Housing, Informal_Settlement,
    Marital_Status, Children_num, Employ_yn, Income_Rand, Grants,
    Imprisonment_Hx, KarnofskyScore, VL_Undetectable, VL_Copies,
    Alcohol, Smoking, SubstanceUse
  ), by = "ID") %>%
  left_join(PRX_Aim_Hosp %>% select(ID, AIM, Hospitalized), by = "ID") %>%
  left_join(PRX_BDQStart %>% select(ID, BDQStart_Difference_days), by = "ID") %>%
  filter(EOT_Outcome != 6)

# Remove duplicate cols before merging PRX_Admission; keep ID
common_vars <- intersect(names(PRX_Admission), names(BDQ_Week))
BDQ_Week <- BDQ_Week %>%
  select(-all_of(setdiff(common_vars, "ID"))) %>%
  left_join(PRX_Admission, by = "ID")

rm(PRX_Admission, PRX3, PRX_BDQStart)

## Outcome: 0 = Successful, 1 = Unsuccessful
BDQ_Week <- BDQ_Week %>%
  mutate(Outcome = case_when(
    EOT_Outcome %in% c(1, 2) ~ 0,
    EOT_Outcome %in% c(3, 4, 5) ~ 1,
    TRUE ~ NA_real_
  ))

## Optional: strata (unused, retained for reference)
BDQ_Week$Strata <- NA_integer_
BDQ_Week$Strata[BDQ_Week$ID >= 150001 & BDQ_Week$ID <= 150301] <- 1L
BDQ_Week$Strata[BDQ_Week$ID >= 150302 & BDQ_Week$ID <= 160000] <- 2L

## Censoring
BDQ_Week <- BDQ_Week %>%
  mutate(
    Censored = ifelse(EOT_Outcome %in% c(4, 5), 1, 0),
    Week_Censored = Week_Outcome
  )

## Mean adherence across weeks (for QC only; not used downstream)
MAX_WEEKS <- 24
all_week_cols <- paste0("Week_", 1:31)
have_week_cols <- intersect(all_week_cols, names(BDQ_Week))
BDQ_Week <- BDQ_Week %>%
  rowwise() %>%
  mutate(total_BDQ_adherence = mean(c_across(all_of(have_week_cols)), na.rm = TRUE)) %>%
  ungroup()

## VL_Undetectable cleanup to 0/1 numeric
BDQ_Week <- BDQ_Week %>%
  mutate(
    VL_Undetectable = ifelse(VL_Undetectable == "Yes" | VL_Copies < 200, "Yes", "No"),
    VL_Undetectable = factor(VL_Undetectable, levels = c("Yes","No"))
  )
BDQ_Week$VL_Undetectable <- as.integer(BDQ_Week$VL_Undetectable == "Yes")

## Marital dummies
BDQ_Week <- BDQ_Week %>%
  mutate(
    Marital_Status = case_when(
      Marital_Status == 1 ~ "Married",
      Marital_Status == 2 ~ "Partnered",
      Marital_Status == 3 ~ "Single",
      TRUE ~ NA_character_
    ),
    Marital_Status = factor(Marital_Status, levels = c("Single","Married","Partnered")),
    Married_dummy   = as.integer(Marital_Status == "Married"),
    Partnered_dummy = as.integer(Marital_Status == "Partnered"),
    Single_dummy    = as.integer(Marital_Status == "Single")
  ) %>%
  select(-Marital_Status)

## ---------------------------------------------------------------------------
## Hospitalization spans (days + weekly markers)
## ---------------------------------------------------------------------------
BDQ_Week <- BDQ_Week %>% mutate(First_Date = as.Date(First_Date))

BDQ_Week <- BDQ_Week %>%
  mutate(
    hosp1_days = ifelse(!is.na(First_Date) & !is.na(Discharge_Date),
                        as.numeric(Discharge_Date - First_Date), 0),
    hosp2_days = ifelse(!is.na(Readmission_Date) & !is.na(Readmission_Discharge_Date),
                        as.numeric(Readmission_Discharge_Date - Readmission_Date), 0),
    hosp3_days = ifelse(!is.na(Readmission_2_Date) & !is.na(Readmission_Discharge_2_Date),
                        as.numeric(Readmission_Discharge_2_Date - Readmission_2_Date), 0),
    total_days_hospitalized = ifelse(Hospitalized %in% c(0, 2), 0,
                                     hosp1_days + hosp2_days + hosp3_days),
    total_days_hospitalized = pmax(total_days_hospitalized, 0)
  )

for (wk in 1:MAX_WEEKS) {
  BDQ_Week[[paste0("Hosp_Week_", wk)]] <- 0L
  in1 <- !is.na(BDQ_Week$First_Date) & !is.na(BDQ_Week$Discharge_Date) &
    wk >= 1 &
    wk <= (floor(as.numeric(BDQ_Week$Discharge_Date - BDQ_Week$First_Date) / 7) + 1)
  in2 <- !is.na(BDQ_Week$Readmission_Date) & !is.na(BDQ_Week$Readmission_Discharge_Date) &
    wk >= (floor(as.numeric(BDQ_Week$Readmission_Date - BDQ_Week$First_Date) / 7) + 1) &
    wk <= (floor(as.numeric(BDQ_Week$Readmission_Discharge_Date - BDQ_Week$First_Date) / 7) + 1)
  in3 <- !is.na(BDQ_Week$Readmission_2_Date) & !is.na(BDQ_Week$Readmission_Discharge_2_Date) &
    wk >= (floor(as.numeric(BDQ_Week$Readmission_2_Date - BDQ_Week$First_Date) / 7) + 1) &
    wk <= (floor(as.numeric(BDQ_Week$Readmission_Discharge_2_Date - BDQ_Week$First_Date) / 7) + 1)
  BDQ_Week[[paste0("Hosp_Week_", wk)]] <- as.integer(in1 | in2 | in3)
}

## ---------------------------------------------------------------------------
## Time-series engineered features (lean: last-3 window only)
## ---------------------------------------------------------------------------
.row_slope <- function(x){
  x <- as.numeric(x); idx <- seq_along(x); keep <- is.finite(x)
  if (sum(keep) < 2) return(NA_real_)
  xi <- idx[keep]; yi <- x[keep]
  cx <- xi - mean(xi); cy <- yi - mean(yi)
  sum(cx * cy) / sum(cx * cx)
}
.row_var <- function(x){
  x <- as.numeric(x); x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  stats::var(x)
}
.tail_vec <- function(x,m){ x <- as.numeric(x); n <- length(x); if(n<1) return(numeric(0)); x[max(1,n-m+1):n] }

compute_ts_features_long <- function(df, max_weeks = MAX_WEEKS){
  ids <- df$ID
  out <- vector("list", length(ids) * max_weeks)
  ctr <- 0L
  for (i in seq_along(ids)) {
    all_weeks <- as.numeric(df[i, paste0("Week_", 1:max_weeks), drop = TRUE])
    for (k in 1:max_weeks) {
      w <- all_weeks[1:k]
      keep <- is.finite(w)
      w3 <- .tail_vec(w, min(3, k))
      mean_last3  <- if (is.atomic(w3)) mean(w3, na.rm = TRUE) else NA_real_
      slope_last3 <- if (is.atomic(w3)) .row_slope(w3) else NA_real_
      var_last3   <- if (is.atomic(w3)) .row_var(w3) else NA_real_
      ctr <- ctr + 1L
      out[[ctr]] <- data.frame(
        ID = ids[i], WeeksUsed = k,
        mean_last3, slope_last3, var_last3,
        stringsAsFactors = FALSE
      )
    }
  }
  dplyr::bind_rows(out)
}

BDQ_TS <- compute_ts_features_long(BDQ_Week, max_weeks = MAX_WEEKS)

## ---------------------------------------------------------------------------
## IPCW (simple pragmatic model) + truncation
## ---------------------------------------------------------------------------
#Finding variables that predict censoring using univariate logistic regression 
#Plan to applied at early censoring (Week 2 for our analysis)
candidate_vars_ipwm <- c(
  "Age",           # baseline
  "Gender",        # baseline
  "Week_1",        # early adherence
  "Hosp_Week_1",   # early hospitalization
  "KarnofskyScore",     # clinical baseline
  "Education",
  "Housing",
  "Single_dummy",
  "Married_dummy",
  "Partnered_dummy",
  "Income_Rand",
  "Grants",
  "Imprisonment_Hx",
  "VL_Undetectable",
  "Alcohol",
  "Smoking",
  "SubstanceUse",
  "BDQStart_Difference_days",
  "AIM", 
  "Hospitalized"
)

# Keep only complete cases for selected variables
ipcw_data <- BDQ_Week %>% 
  select(Censored, all_of(candidate_vars_ipwm)) %>% 
  drop_na()

# Fit logistic model
ipcw_model <- glm(Censored ~ ., data = ipcw_data, family = binomial)

# View summary
summary(ipcw_model)

table_na(BDQ_Week$Censored, useNA = "ifany")
table_na(BDQ_Week$Week_Censored)
table_na(BDQ_Week$ID[is.na(BDQ_Week$Week_Censored)])
table_na(BDQ_Week$Outcome[is.na(BDQ_Week$Week_Censored)])
table_na(BDQ_Week$EOT_Date[is.na(BDQ_Week$Week_Censored)])

#Input missing values for Week_Censored
table_na(BDQ_Week$Week_Censored) #1 missing 
BDQ_Week <- BDQ_Week %>%
  mutate(Week_Censored_in = ifelse(is.na(Week_Censored), 24, Week_Censored))
table_na(BDQ_Week$Week_Censored, useNA = "ifany")
#Input missing alues for income_rand
BDQ_Week <- BDQ_Week %>%
  mutate(Income_Rand_in = ifelse(is.na(Income_Rand), 0, Income_Rand))
#Need to impute VL undetectable 
BDQ_Week <- BDQ_Week %>%
  mutate(VL_Undetectable_in = ifelse(is.na(VL_Undetectable), "No", VL_Undetectable)) 


# Fit IPCW logistic model using early predictors
ipcw_model <- glm(Censored ~ Week_1 + Hosp_Week_1 + Income_Rand_in + VL_Undetectable_in + Hospitalized,
                  data = BDQ_Week, family = binomial)

BDQ_Week <- BDQ_Week %>%
  mutate(
    Censor_Prob = predict(ipcw_model, type = "response"),
    IPCW = 1 / pmax(1 - Censor_Prob, 1e-6)
  )
BDQ_Week$IPCW <- pmin(BDQ_Week$IPCW, quantile(BDQ_Week$IPCW, 0.90, na.rm = TRUE))

## ---------------------------------------------------------------------------
## Temporal train/valid split (70/30 within AIM, by First_Date)
## ---------------------------------------------------------------------------
set.seed(123)
df_split <- BDQ_Week %>%
  select(ID, AIM, Outcome, First_Date) %>%
  distinct() %>%
  mutate(
    First_Date = as.Date(First_Date),
    First_Date = if_else(is.na(First_Date),
                         min(First_Date, na.rm = TRUE) - 1,
                         First_Date)
  ) %>%
  arrange(AIM, First_Date)

cut_table <- df_split %>%
  group_by(AIM) %>%
  summarize(
    cutoff = as.Date(quantile(as.numeric(First_Date), probs = 0.70, na.rm = TRUE),
                     origin = "1970-01-01"),
    .groups = "drop"
  )

df_with_cut <- df_split %>%
  left_join(cut_table, by = "AIM") %>%
  mutate(split_group = if_else(First_Date <= cutoff, "train", "valid"))

train_ids <- df_with_cut %>% filter(split_group == "train") %>% pull(ID)
valid_ids <- df_with_cut %>% filter(split_group == "valid") %>% pull(ID)

BDQ_train <- BDQ_Week %>% filter(ID %in% train_ids)
BDQ_valid <- BDQ_Week %>% filter(ID %in% valid_ids)
stopifnot(length(intersect(BDQ_train$ID, BDQ_valid$ID)) == 0)

cat("\nOutcome distribution (train vs valid):\n")
print(prop.table(table(BDQ_train$Outcome)))
print(prop.table(table(BDQ_valid$Outcome)))

## ---------------------------------------------------------------------------
## MODEL CONFIG
## ---------------------------------------------------------------------------
covariates <- c(
  "Age","Gender","BMI","Single_dummy","Married_dummy",
  "KarnofskyScore","VL_Undetectable","Alcohol",
  "Hospitalized","BDQStart_Difference_days"
)

MAX_WEEKS   <- 24
MODEL_WEEKS <- 0:24
TUNE_WEEKS  <- 0:24

N_FOLDS          <- 10
N_RANDOM_TRIALS  <- 100
NROUNDS_LIMIT    <- 1000
EARLY_STOPPING   <- 100
SEED_BASE        <- 2025

EARLY_FOCUS_WEEKS <- 1:4
EARLY_WEIGHT_MULT <- 4.0

## ---------------------------------------------------------------------------
## Build horizon-specific design frame (lean; no interactions)
## ---------------------------------------------------------------------------
build_lagged_train <- function(BDQ_df, n_weeks, covariates) {
  stopifnot(n_weeks >= 0)
  if (n_weeks == 0) {
    model_vars <- c("ID","Outcome","IPCW", covariates)
    out <- BDQ_df %>%
      dplyr::select(dplyr::all_of(model_vars)) %>%
      tidyr::drop_na(dplyr::all_of(c("Outcome","IPCW",covariates)))
    return(out)
  }
  
  # Keep only current and two most recent adherence weeks 
  recent_idx  <- seq.int(max(1, n_weeks - 2), n_weeks)
  week_vars   <- paste0("Week_", recent_idx)
  have_weeks  <- intersect(week_vars, names(BDQ_df))
  
  model_vars  <- c("ID","Outcome","IPCW", covariates, have_weeks)
  base <- BDQ_df %>% dplyr::select(dplyr::all_of(model_vars))
  
  # Engineered features 
  feats_k <- BDQ_TS %>%
    dplyr::filter(WeeksUsed == n_weeks) %>%
    dplyr::select(-WeeksUsed)
  
  keep_eng <- grep("(last3)$|(last3_)", names(feats_k), value = TRUE)
  feats_k  <- dplyr::select(feats_k, dplyr::any_of(c("ID", keep_eng)))
  if (n_weeks < 3) {
    feats_k <- dplyr::select(feats_k, -tidyselect::contains("last3"))
  }
  
  out <- base %>% dplyr::left_join(feats_k, by = "ID")
  
  # Hospitalization flags
  hosp_this_col <- paste0("Hosp_Week_", n_weeks)
  hosp_prev_col <- if (n_weeks > 1) paste0("Hosp_Week_", n_weeks - 1) else NULL
  
  h_this <- if (hosp_this_col %in% names(BDQ_df)) BDQ_df[[hosp_this_col]] else rep(NA_integer_, nrow(BDQ_df))
  h_prev <- if (!is.null(hosp_prev_col) && hosp_prev_col %in% names(BDQ_df)) BDQ_df[[hosp_prev_col]] else rep(NA_integer_, nrow(BDQ_df))
  
  if (any(!is.na(h_this)) || any(!is.na(h_prev))) {
    hosp_df <- dplyr::tibble(
      ID = BDQ_df$ID,
      Hosp_ThisWeek = dplyr::case_when(
        h_this == 1 ~ 1L,
        h_this == 0 ~ 0L,
        TRUE        ~ NA_integer_
      ),
      Hosp_LastWeek = dplyr::case_when(
        h_prev == 1 ~ 1L,
        h_prev == 0 ~ 0L,
        TRUE        ~ NA_integer_
      )
    ) %>%
      dplyr::mutate(
        Hosp_ThisOrLast = as.integer(
          (dplyr::coalesce(Hosp_ThisWeek, 0L) == 1L) |
            (dplyr::coalesce(Hosp_LastWeek, 0L) == 1L)
        )
      )
    out <- out %>% dplyr::left_join(hosp_df, by = "ID")
  }
  
  # Required columns
  required_cols <- c("Outcome","IPCW", covariates, have_weeks)
  if ("Hosp_ThisOrLast" %in% names(out)) {
    required_cols <- c(required_cols, "Hosp_ThisOrLast")
  }
  
  out <- tidyr::drop_na(out, dplyr::all_of(required_cols))
  out
}

## ---------------------------------------------------------------------------
## Combined tuning matrix with early-horizon upweighting
## ---------------------------------------------------------------------------
build_combined_tuning <- function(BDQ_train, covariates,
                                  tune_weeks = TUNE_WEEKS,
                                  early_weeks = EARLY_FOCUS_WEEKS,
                                  early_weight_mult = EARLY_WEIGHT_MULT) {
  stopifnot(early_weight_mult >= 1)
  
  built <- lapply(tune_weeks, function(k) {
    dfk <- build_lagged_train(BDQ_train, k, covariates)
    if (!nrow(dfk)) return(NULL)
    
    dfk$.rowid <- seq_len(nrow(dfk))
    dfk <- dfk[is.finite(dfk$IPCW), , drop = FALSE]
    if (!nrow(dfk)) return(NULL)
    
    mf <- model.frame(Outcome ~ . - ID - IPCW - .rowid, data = dfk, na.action = na.omit)
    if (!nrow(mf)) return(NULL)
    
    X <- model.matrix(~ . - 1, data = mf)
    y <- as.numeric(model.response(mf))
    
    kept_rowids <- model.frame(~ .rowid, data = dfk, na.action = na.omit)$`.rowid`[as.integer(rownames(mf))]
    w <- dfk$IPCW[kept_rowids]
    
    cap <- stats::quantile(w, 0.90, na.rm = TRUE)
    w[w > cap] <- cap
    if (k %in% early_weeks) w <- w * early_weight_mult
    w <- w / mean(w, na.rm = TRUE)
    
    list(k = k, X = X, y = y, w = as.numeric(w),
         id = dfk$ID[kept_rowids])
  })
  
  built <- Filter(Negate(is.null), built)
  stopifnot(length(built) > 0)
  
  all_cols <- Reduce(union, lapply(built, function(b) colnames(b$X)))
  X_list <- lapply(built, function(b){
    miss <- setdiff(all_cols, colnames(b$X))
    if (length(miss)) {
      add <- matrix(0, nrow = nrow(b$X), ncol = length(miss))
      colnames(add) <- miss
      b$X <- cbind(b$X, add)
    }
    b$X[, all_cols, drop = FALSE]
  })
  
  X <- do.call(rbind, X_list)
  y <- unlist(lapply(built, `[[`, "y"))
  w <- unlist(lapply(built, `[[`, "w"))
  id_per_row <- unlist(lapply(built, `[[`, "id"))
  
  list(
    X = Matrix(X, sparse = TRUE),
    y = as.numeric(y),
    w = as.numeric(w),
    cols = all_cols,
    id = id_per_row
  )
}

drop_all_na_cols <- function(df, keep = c("Outcome","IPCW","ID")) {
  nz <- vapply(df, function(x) !all(is.na(x)), logical(1))
  # Always keep key cols if present
  nz[names(nz) %in% keep] <- TRUE
  df[, nz, drop = FALSE]
}

drop_all_constant_cols <- function(df, keep = c("Outcome","IPCW","ID")) {
  is_const <- vapply(df, function(x) {
    # treat factors/characters as constant if they have 0-1 unique non-NA values
    if (is.factor(x) || is.character(x)) {
      ux <- unique(x[!is.na(x)])
      return(length(ux) <= 1L)
    }
    y <- x[is.finite(x)]  # ignore NA/Inf for const check
    if (length(y) == 0L) return(TRUE)
    length(unique(y)) <= 1L
  }, logical(1))
  # never drop key columns
  is_const[names(is_const) %in% keep] <- FALSE
  df[, !is_const, drop = FALSE]
}


## ---------------------------------------------------------------------------
## Global tuning with ID-grouped, stratified folds
## ---------------------------------------------------------------------------
tune_data    <- build_combined_tuning(BDQ_train, covariates)
dtrain_tune  <- xgb.DMatrix(data = tune_data$X, label = tune_data$y, weight = tune_data$w)

# ID-level label: any positive across that person’s rows?
id_label <- tapply(tune_data$y, tune_data$id, function(v) as.integer(any(v == 1)))
id_vec   <- names(id_label)

# Stratified folds on IDs (balanced by id_label)
set.seed(SEED_BASE)
k_safe <- max(2, min(N_FOLDS, sum(id_label == 0), sum(id_label == 1)))
id_folds <- caret::createFolds(factor(id_label), k = k_safe, list = TRUE, returnTrain = FALSE)

# Map ID-based folds back to row indices in the stacked matrix
folds <- lapply(id_folds, function(id_idx){
  test_ids <- id_vec[id_idx]
  which(tune_data$id %in% test_ids)
})

# Sanity: each test fold must have both classes (by rows)
y_all <- tune_data$y
ok <- sapply(folds, function(idx) length(unique(y_all[idx])) == 2)
if (!all(ok)) {
  bad <- which(!ok)
  stop(sprintf("Degenerate test folds (single class) after grouping: %s", paste(bad, collapse=", ")))
}

# Random search space
rand_grid <- function(n = N_RANDOM_TRIALS, seed = SEED_BASE) {
  set.seed(seed)
  data.frame(
    eta               = 10^runif(n, log10(0.02), log10(0.12)),
    max_depth         = sample(2:8, n, replace = TRUE),
    min_child_weight  = sample(c(1:8), n, replace = TRUE),
    subsample         = runif(n, 0.6, 1),
    colsample_bytree  = runif(n, 0.6, 1),
    lambda            = 10^runif(n, log10(0.10), log10(10)),
    alpha             = 10^runif(n, log10(0.01), log10(1)),
    gamma             = sample(c(0, 0.25, 0.5, 1), n, replace = TRUE),
    max_delta_step    = sample(c(0, 1, 2), n, replace = TRUE)
  )
}
grid <- rand_grid()

## ---- Run CV with robust metric extraction (AUPRC primary) ----
cv_summaries <- lapply(seq_len(nrow(grid)), function(i){
  p <- grid[i, ]
  params <- list(
    objective        = "binary:logistic",
    eval_metric      = c("aucpr","auc"),
    eta              = p$eta,
    max_depth        = p$max_depth,
    min_child_weight = p$min_child_weight,
    subsample        = p$subsample,
    colsample_bytree = p$colsample_bytree,
    lambda           = p$lambda,  # or reg_lambda
    alpha            = p$alpha,   # or reg_alpha
    gamma            = p$gamma,
    scale_pos_weight = 1,         # IPCW already provided as DMatrix weight
    tree_method      = "hist",
    max_bin          = 256,
    verbosity        = 0
  )
  
  set.seed(SEED_BASE + i)
  cv <- xgb.cv(
    params                = params,
    data                  = dtrain_tune,
    nrounds               = NROUNDS_LIMIT,
    folds                 = folds,
    stratified            = FALSE,   # we already grouped/stratified by ID
    verbose               = FALSE,
    early_stopping_rounds = EARLY_STOPPING,
    showsd                = TRUE
  )
  
  ev <- cv$evaluation_log
  
  get_metric <- function(ev, metric = c("aucpr","auc")) {
    metric <- match.arg(metric)
    pats <- c(
      paste0("^test-",  metric, "-mean$"),
      paste0("^test_",  metric, "_mean$"),
      paste0("^test\\.",metric,"\\.mean$")
    )
    for (pat in pats) {
      cols <- grep(pat, names(ev), value = TRUE)
      if (length(cols) >= 1L) return(max(ev[[cols[1]]], na.rm = TRUE))
    }
    NA_real_
  }
  
  data.frame(
    trial          = i,
    best_iteration = cv$best_iteration,
    aucpr_mean     = get_metric(ev, "aucpr"),
    auc_mean       = get_metric(ev, "auc")
  )
})

cv_summaries <- bind_rows(cv_summaries)

best_idx     <- with(cv_summaries, order(-aucpr_mean, -auc_mean, -best_iteration))[1]
best_trial   <- cv_summaries$trial[best_idx]
best_iter <- cv_summaries$best_iteration[best_idx]

nrounds_use <- min(
  max(50, ceiling(3 * best_iter)),  
  NROUNDS_LIMIT                     
)

GLOBAL_NROUNDS <- nrounds_use

best_params  <- as.list(grid[best_trial, ])

print(best_params)

GLOBAL_PARAMS <- list(
  objective         = "binary:logistic",
  eval_metric       = c("aucpr","auc"),
  eta               = best_params$eta,
  max_depth         = best_params$max_depth,
  min_child_weight  = best_params$min_child_weight,
  subsample         = best_params$subsample,
  colsample_bytree  = best_params$colsample_bytree,
  lambda            = best_params$lambda,  # or reg_lambda
  alpha             = best_params$alpha,   # or reg_alpha
  scale_pos_weight  = 1,                   # DO NOT double-weight with IPCW
  tree_method       = "hist",
  max_bin           = 256,
  verbosity         = 0
)

message("Global tuned (trial ", best_trial, "):")
print(cv_summaries[best_idx, ])
message("Global best nrounds (guarded): ", GLOBAL_NROUNDS)

## Train global model
global_model <- xgb.train(params = GLOBAL_PARAMS, data = dtrain_tune, nrounds = GLOBAL_NROUNDS, verbose = 0)
preds_prob_global <- predict(global_model, tune_data$X)

pick_thr_f1 <- function(prob, y, grid = seq(0.02, 0.98, by = 0.001),
                        posrate_bounds = c(0.05, 0.95),
                        min_spec = 0.50, min_sens = 0.50) {
  y <- as.integer(y)
  calc_row <- function(thr) {
    pred <- as.integer(prob >= thr)
    tp <- sum(pred==1 & y==1); fp <- sum(pred==1 & y==0)
    fn <- sum(pred==0 & y==1); tn <- sum(pred==0 & y==0)
    prec <- if ((tp+fp)>0) tp/(tp+fp) else 0
    rec  <- if ((tp+fn)>0) tp/(tp+fn) else 0
    f1   <- if ((prec+rec)>0) 2*prec*rec/(prec+rec) else 0
    spec <- if ((tn+fp)>0) tn/(tn+fp) else 0
    posrate <- mean(pred==1)
    c(thr=thr, f1=f1, sens=rec, spec=spec, posrate=posrate)
  }
  M <- t(vapply(grid, calc_row, numeric(5)))
  ok <- which(M[, "posrate"] >= posrate_bounds[1] &
                M[, "posrate"] <= posrate_bounds[2] &
                M[, "spec"]    >= min_spec &
                M[, "sens"]    >= min_sens)
  best <- if (length(ok)) {
    f1max <- max(M[ok, "f1"])
    cand  <- ok[M[ok, "f1"] >= f1max - 1e-12]
    cand[which.max(M[cand, "sens"])]
  } else which.max(M[, "f1"])
  as.numeric(M[best, "thr"])
}

best_thresh <- pick_thr_f1(preds_prob_global, tune_data$y)
print(best_thresh)

BLOCKLIST_REGEX <- paste(
  "^ID$",                 # identifier
  "^IPCW$",               # sample weight
  "^\\.rowid$",           # helper index
  "^Censored$",           # censoring outcome
  "^Censor_Prob$",        # model-derived prob of censoring
  "^Week_Censored",       # any Week_Censored* columns
  "^Week_Censored_in$",
  "^Income_Rand_in$",
  "^VL_Undetectable_in$",
  "^Hosp_ThisWeek$",
  "^Hosp_LastWeek$",
  "^Hosp_ThisOrLast$",
  sep = "|"
)

###############################################################################
### Validation Data ####
event_counts_valid <- lapply(0:24, function(k) {
  # Build horizon-specific validation frame
  lagged_valid <- build_lagged_train(BDQ_valid, k, covariates) |>
    drop_all_na_cols() |>
    drop_all_constant_cols()
  
  # If nothing left, return an empty summary row for this k
  if (!nrow(lagged_valid)) {
    return(data.frame(
      WeeksUsed  = k, N_valid = 0, Events = NA_integer_, NonEvents = NA_integer_,
      EventRate  = NA_real_, TwoClasses = 0L, stringsAsFactors = FALSE
    ))
  }
  
  # Align to modeling columns/rows; drop ID/IPCW from RHS
  mf_va <- try(
    stats::model.frame(Outcome ~ . - ID - IPCW,
                       data = lagged_valid, na.action = stats::na.omit),
    silent = TRUE
  )
  
  # If model.frame failed or dropped everything, return empty summary
  if (inherits(mf_va, "try-error") || !nrow(mf_va)) {
    return(data.frame(
      WeeksUsed  = k, N_valid = 0, Events = NA_integer_, NonEvents = NA_integer_,
      EventRate  = NA_real_, TwoClasses = 0L, stringsAsFactors = FALSE
    ))
  }
  
  y_valid <- model.response(mf_va)
  n_all   <- length(y_valid)
  n_pos   <- sum(y_valid == 1, na.rm = TRUE)
  n_neg   <- sum(y_valid == 0, na.rm = TRUE)
  posrate <- if (n_all > 0) n_pos / n_all else NA_real_
  two_cls <- as.integer(length(unique(y_valid)) == 2)
  
  data.frame(
    WeeksUsed  = k,
    N_valid    = n_all,
    Events     = n_pos,
    NonEvents  = n_neg,
    EventRate  = posrate,
    TwoClasses = two_cls,
    stringsAsFactors = FALSE
  )
})

event_counts_valid <- dplyr::bind_rows(event_counts_valid)

print(event_counts_valid)

#### Training data ####
event_counts_train <- lapply(0:24, function(k) {
  # Build horizon-specific training frame and drop columns with all NA
  lagged_train <- build_lagged_train(BDQ_train, k, covariates) |> drop_all_na_cols()
  # Align to the same rows used in modeling (via model.frame)
  mf_tr   <- stats::model.frame(Outcome ~ . - ID - IPCW, data = lagged_train, na.action = stats::na.omit)
  y_train <- model.response(mf_tr)
  
  n_all   <- length(y_train)
  n_pos   <- sum(y_train == 1, na.rm = TRUE)
  n_neg   <- sum(y_train == 0, na.rm = TRUE)
  posrate <- if (n_all > 0) n_pos / n_all else NA_real_
  
  data.frame(
    WeeksUsed   = k,
    N_train     = n_all,
    Events      = n_pos,
    NonEvents   = n_neg,
    EventRate   = posrate,
    TwoClasses  = as.integer(length(unique(y_train)) == 2),
    stringsAsFactors = FALSE
  )
})

events_by_week_train <- dplyr::bind_rows(event_counts_train)

# Look at the training table
print(events_by_week_train)

#### Helper ######################
cv_platt_calibrate <- function(prob, y, K = 5, seed = 1) {
  stopifnot(length(prob) == length(y))
  prob <- pmin(pmax(prob, 1e-6), 1 - 1e-6)
  set.seed(seed)
  f <- caret::createFolds(y, k = K, list = TRUE)
  p_cal <- rep(NA_real_, length(y))
  for (k in seq_along(f)) {
    idx_te <- f[[k]]; idx_tr <- setdiff(seq_along(y), idx_te)
    logit_p <- qlogis(prob[idx_tr])
    fit <- glm(y[idx_tr] ~ logit_p, family = binomial())
    p_cal[idx_te] <- plogis(predict(fit, newdata = data.frame(logit_p = qlogis(prob[idx_te]))))
  }
  p_cal
}

## ---------------------------------------------------------------------------
## Main week-by-week training & evaluation 
## ---------------------------------------------------------------------------
results <- data.frame(
  WeeksUsed = MODEL_WEEKS,
  AUC = NA, AUC_LCI = NA, AUC_UCI = NA,
  Sensitivity = NA, Sensitivity_LCI = NA, Sensitivity_UCI = NA,
  Specificity = NA, Specificity_LCI = NA, Specificity_UCI = NA,
  F1_Score = NA, F1_LCI = NA, F1_UCI = NA,
  PPV = NA, PPV_LCI = NA, PPV_UCI = NA,
  NPV = NA, NPV_LCI = NA, NPV_UCI = NA,
  Brier_Cal = NA, Brier_LCI = NA, Brier_UCI = NA,
  Accuracy = NA, Accuracy_LCI = NA, Accuracy_UCI = NA,
  Brier_Null = NA, Brier_RelativeReduction = NA,
  MAE_Cal = NA, RMSE_Cal = NA, AUPRC = NA,
  SingleClassValid = NA,
  stringsAsFactors = FALSE
)

# Save artifacts for specific weeks
SAVE_WEEKS <- c(0, 4, 10)
wk_artifacts <- list()
preds_calibrated_wk4 <- NULL
y_valid_wk4          <- NULL
week4_valid          <- NULL           
preds_prob_wk4_ens      <- NULL        
y_valid_wk4_aligned     <- NULL
preds_prob_wk10_ens     <- NULL       
y_valid_wk10_aligned    <- NULL
best_thresh <- NA_real_ 

for (n_weeks in MODEL_WEEKS) {
  # Build TRAIN / VALID frames
  lagged_train <- build_lagged_train(BDQ_train, n_weeks, covariates) %>% drop_all_na_cols()
  lagged_valid <- build_lagged_train(BDQ_valid, n_weeks, covariates) %>% drop_all_na_cols()
  
  if (!"IPCW" %in% names(lagged_train)) lagged_train$IPCW <- 1
  if (!"IPCW" %in% names(lagged_valid)) lagged_valid$IPCW <- 1
  if (n_weeks == 0) {
    lagged_train$IPCW <- 1
    lagged_valid$IPCW <- 1
  }
  
  mf_tr   <- model.frame(Outcome ~ . - ID - IPCW, data = lagged_train, na.action = na.omit)
  X_train <- model.matrix(~ . - Outcome - 1, data = mf_tr)
  y_train <- model.response(mf_tr); kept_tr <- as.integer(rownames(mf_tr))
  
  mf_va   <- model.frame(Outcome ~ . - ID - IPCW, data = lagged_valid, na.action = na.omit)
  X_valid <- model.matrix(~ . - Outcome - 1, data = mf_va)
  y_valid <- model.response(mf_va)
  
  kept_va <- as.integer(rownames(mf_va)) 
  week_valid_aligned <- lagged_valid[kept_va, , drop = FALSE]
  
  ## IPCW weights + mild sex/label balancing
  w_ipcw <- lagged_train$IPCW[kept_tr]
  w_ipcw[!is.finite(w_ipcw)] <- 1
  
  sex_kept <- lagged_train$Gender[kept_tr]
  MALE_TARGET <- 0.5; FEMALE_TARGET <- 0.5
  tab_sex <- table(sex_kept); p_m <- as.numeric(tab_sex["0"])/sum(tab_sex); p_f <- as.numeric(tab_sex["1"])/sum(tab_sex)
  w_sex <- rep(1, length(sex_kept))
  if (!is.na(p_m) && p_m > 0) w_sex[sex_kept == 0] <- (MALE_TARGET   / p_m) 
  if (!is.na(p_f) && p_f > 0) w_sex[sex_kept == 1] <- (FEMALE_TARGET / p_f) 
  
  y_kept   <- y_train
  key <- paste0("g", sex_kept, "_y", y_kept)
  counts <- table(key); target <- mean(counts)
  w_labsex <- target / as.numeric(counts[key])
  
  keep_ipcw <- !is.na(w_ipcw)
  w_train <- w_ipcw[keep_ipcw] * w_sex[keep_ipcw] * w_labsex[keep_ipcw]
  w_train <- pmin(w_train, stats::quantile(w_train, 0.90, na.rm = TRUE))
  w_train <- w_train / mean(w_train, na.rm = TRUE)
  
  # Align columns
  all_cols <- union(colnames(X_train), colnames(X_valid))
  miss <- setdiff(all_cols, colnames(X_train))
  if (length(miss)) {
    add <- matrix(0, nrow = nrow(X_train), ncol = length(miss)); colnames(add) <- miss
    X_train <- cbind(X_train, add)
  }
  miss <- setdiff(all_cols, colnames(X_valid))
  if (length(miss)) {
    add <- matrix(0, nrow = nrow(X_valid), ncol = length(miss)); colnames(add) <- miss
    X_valid <- cbind(X_valid, add)
  }
  X_train <- X_train[, all_cols, drop = FALSE]
  X_valid <- X_valid[, all_cols, drop = FALSE]
  
  # Blocklist (no identifiers/weights as features)
  bad <- grep(BLOCKLIST_REGEX, colnames(X_train), ignore.case = TRUE, value = TRUE)
  if (length(bad)) {
    X_train <- X_train[, setdiff(colnames(X_train), bad), drop = FALSE]
    X_valid <- X_valid[, setdiff(colnames(X_valid), bad), drop = FALSE]
  }
  common_cols <- intersect(colnames(X_train), colnames(X_valid))
  X_train <- X_train[, common_cols, drop = FALSE]
  X_valid <- X_valid[, common_cols, drop = FALSE]
  
  # XGBoost
  dtrain <- xgb.DMatrix(X_train, label = y_train, weight = w_train)
  dvalid <- xgb.DMatrix(X_valid, label = y_valid, weight = rep(1, length(y_valid)))
  
  neg_wk <- sum(y_train == 0); pos_wk <- sum(y_train == 1)
  spw_wk <- if (pos_wk > 0) neg_wk / pos_wk else 1
  
  local_params <- GLOBAL_PARAMS
  local_params$scale_pos_weight <- spw_wk
  local_params$max_delta_step <- 1
  local_params$gamma <- 0.5
  
  set.seed(123 + n_weeks)
  watchlist <- list(valid = dvalid)
  model <- xgb.train(
    params = local_params, data = dtrain,
    nrounds = min(GLOBAL_NROUNDS * 2, NROUNDS_LIMIT),
    watchlist = watchlist, early_stopping_rounds = 1000, verbose = 0
  )
  best_iter <- if (!is.null(model$best_iteration)) model$best_iteration else GLOBAL_NROUNDS
  
  # small ensemble for stability
  if (n_weeks %in% SAVE_WEEKS) wk_models <- list()
  
  SEEDS <- 2025 + 0:4
  pred_mat <- sapply(SEEDS, function(s, .save = n_weeks %in% SAVE_WEEKS){
    set.seed(s + n_weeks)
    m <- xgb.train(
      params = local_params, data = dtrain,
      nrounds = min(GLOBAL_NROUNDS * 2, NROUNDS_LIMIT),
      watchlist = watchlist, early_stopping_rounds = 100, verbose = 0
    )
    if (.save) wk_models[[as.character(s)]] <<- m  
    bi <- if (!is.null(m$best_iteration)) m$best_iteration else GLOBAL_NROUNDS
    predict(m, newdata = dvalid, ntreelimit = bi)
  })
  
  # collapse ensemble to per-row mean probabilities
  if (is.null(dim(pred_mat))) {
    pred_mat <- matrix(pred_mat, ncol = 1)
  }
  if (nrow(pred_mat) != length(y_valid)) {
    stop(sprintf("Length mismatch: nrow(pred_mat)=%d, length(y_valid)=%d",
                 nrow(pred_mat), length(y_valid)))
  }
  preds_prob <- rowMeans(pred_mat)  
  
  # Save artifacts for Week 4 / Week 10 
  if (n_weeks %in% SAVE_WEEKS) {
    wk_artifacts[[as.character(n_weeks)]] <- list(
      y_valid    = as.integer(y_valid),
      preds_prob = preds_prob,
      X_valid    = X_valid,
      model_cols = colnames(X_train),
      models     = wk_models
    )
  }
  
  if (n_weeks == 4) {
    preds_calibrated_wk4  <- preds_calibrated
    y_valid_wk4           <- truth
    week4_valid           <- week_valid_aligned  
    preds_prob_wk4_ens    <- preds_prob
    y_valid_wk4_aligned   <- truth
    best_thresh           <- threshold            
  }
  
  if (n_weeks == 10) {
    preds_prob_wk10_ens   <- preds_prob
    y_valid_wk10_aligned  <- truth
  }
  
  ## ---- CALIBRATE FIRST ---------------------------------------------------
  truth <- as.integer(y_valid)                  
  gender_valid <- model.frame(~ Gender, data = mf_va)$Gender
  
  # Sex-specific Platt (logistic) calibration on the VALID fold
  preds_f <- preds_prob[gender_valid == 1]; truth_f <- truth[gender_valid == 1]
  preds_m <- preds_prob[gender_valid == 0]; truth_m <- truth[gender_valid == 0]
  
  preds_calibrated <- numeric(length(preds_prob))
  if (length(truth_f) >= 5 && length(unique(truth_f)) >= 2) {
    cal_f <- glm(truth_f ~ preds_f, family = binomial())
    preds_calibrated[gender_valid == 1] <- predict(cal_f, type = "response")
  } else {
    preds_calibrated[gender_valid == 1] <- preds_f
  }
  if (length(truth_m) >= 5 && length(unique(truth_m)) >= 2) {
    cal_m <- glm(truth_m ~ preds_m, family = binomial())
    preds_calibrated[gender_valid == 0] <- predict(cal_m, type = "response")
  } else {
    preds_calibrated[gender_valid == 0] <- preds_m
  }

  prob_for_metrics <- preds_calibrated  
  
  ## Optionally save week-specific artifacts
  if (n_weeks %in% SAVE_WEEKS) {
    wk_artifacts[[as.character(n_weeks)]] <- list(
      y_valid    = as.integer(truth),
      preds_prob = prob_for_metrics,   
      X_valid    = X_valid,
      model_cols = colnames(X_train),
      models     = wk_models
    )
  }
  
  ## ---- THRESHOLD (on calibrated probs) -----------------------------------
  has_two <- length(unique(truth)) == 2
  res_row <- which(results$WeeksUsed == n_weeks)
  
  # --- Constrained Youden 
  threshold <- if (has_two) {
     roc_obj <- pROC::roc(truth, prob_for_metrics, quiet = TRUE, direction = "<")
     grid <- pROC::coords(roc_obj, x = roc_obj$thresholds, input = "threshold",
                          ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
     df <- as.data.frame(grid)
     df$sens <- as.numeric(df$sensitivity); df$spec <- as.numeric(df$specificity)
     df$J <- df$sens + df$spec - 1
     cand <- subset(df, sens >= 0.7 & spec >= 0.7)
     if (nrow(cand)) cand$threshold[which.max(cand$J)] else best_thresh
   } else best_thresh
  
  if (!is.finite(threshold)) threshold <- 0.5
  
  ## ---- CLASSIFY + REPORT ---------------------
  preds_class <- as.integer(prob_for_metrics >= threshold)
  
  tp <- sum(preds_class==1 & truth==1); fp <- sum(preds_class==1 & truth==0)
  fn <- sum(preds_class==0 & truth==1); tn <- sum(preds_class==0 & truth==0)
  
  sens <- if ((tp+fn)>0) tp/(tp+fn) else NA_real_
  spec <- if ((tn+fp)>0) tn/(tn+fp) else NA_real_
  ppv  <- if ((tp+fp)>0) tp/(tp+fp) else NA_real_
  npv  <- if ((tn+fn)>0) tn/(tn+fn) else NA_real_
  f1   <- if (is.finite(ppv) && is.finite(sens) && (ppv+sens)>0) 2*ppv*sens/(ppv+sens) else NA_real_
  acc  <- mean(preds_class == truth)
  
  # AUC / AUPRC also on calibrated probs
  auc_val <- NA_real_; pr_auc <- NA_real_
  if (has_two) {
    roc_obj <- pROC::roc(truth, prob_for_metrics, quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))
    pr <- tryCatch(
      PRROC::pr.curve(scores.class0 = prob_for_metrics[truth==1],
                      scores.class1 = prob_for_metrics[truth==0], curve = FALSE),
      error = function(e) NULL
    )
    if (!is.null(pr)) pr_auc <- pr$auc.integral
  }
  
  # Calibration metrics
  brier_cal  <- mean((prob_for_metrics - truth)^2)
  mae_cal    <- Metrics::mae(truth, prob_for_metrics)
  rmse_cal   <- Metrics::rmse(truth, prob_for_metrics)
  brier_null <- mean((mean(truth) - truth)^2)
  brier_rr   <- (brier_null - brier_cal) / brier_null
  
  ## ---- 4) BOOTSTRAP CIs------------------
  BOOT_R <- 1000
  if (has_two) {
    df_boot <- data.frame(truth = truth, prob = prob_for_metrics)
    boot_metrics <- function(dat, idx) {
      t_truth <- dat$truth[idx]; t_prob <- dat$prob[idx]
      if (length(unique(t_truth)) < 2L) return(rep(NA_real_, 9))
      t_pred <- as.integer(t_prob >= threshold)  # SAME threshold on calibrated probs
      tp <- sum(t_pred==1 & t_truth==1); fp <- sum(t_pred==1 & t_truth==0)
      fn <- sum(t_pred==0 & t_truth==1); tn <- sum(t_pred==0 & t_truth==0)
      sens_b <- if ((tp+fn)>0) tp/(tp+fn) else NA_real_
      spec_b <- if ((tn+fp)>0) tn/(tn+fp) else NA_real_
      ppv_b  <- if ((tp+fp)>0) tp/(tp+fp) else NA_real_
      npv_b  <- if ((tn+fn)>0) tn/(tn+fn) else NA_real_
      f1_b   <- if (is.finite(ppv_b) && is.finite(sens_b) && (ppv_b + sens_b)>0) 2*ppv_b*sens_b/(ppv_b + sens_b) else NA_real_
      acc_b  <- mean(t_pred == t_truth)
      auc_b  <- suppressWarnings(as.numeric(pROC::auc(pROC::roc(t_truth, t_prob, quiet = TRUE))))
      brier_b <- mean((t_prob - t_truth)^2)
      pr_b <- tryCatch(
        PRROC::pr.curve(scores.class0 = t_prob[t_truth==1],
                        scores.class1 = t_prob[t_truth==0], curve = FALSE)$auc.integral,
        error = function(e) NA_real_
      )
      c(auc_b, sens_b, spec_b, f1_b, ppv_b, npv_b, acc_b, brier_b, pr_b)
    }
    set.seed(123 + n_weeks)
    boot_out <- boot::boot(df_boot, boot_metrics, R = BOOT_R, strata = df_boot$truth)
    q <- function(col) stats::quantile(boot_out$t[,col], c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    ci_auc   <- q(1); ci_sens <- q(2); ci_spec <- q(3); ci_f1 <- q(4)
    ci_ppv   <- q(5); ci_npv  <- q(6); ci_acc  <- q(7); ci_brier <- q(8); ci_auprc <- q(9)
  } else {
    ci_auc <- ci_sens <- ci_spec <- ci_f1 <- ci_ppv <- ci_npv <- ci_acc <- ci_brier <- ci_auprc <- c(NA_real_, NA_real_)
  }
  
  ## ---- WRITE OUT---------------------------------------------
  results$SingleClassValid[res_row] <- as.integer(!has_two)
  results[res_row, c("AUC","Sensitivity","Specificity","Accuracy","F1_Score",
                     "Brier_Cal","Brier_Null","Brier_RelativeReduction",
                     "MAE_Cal","RMSE_Cal","AUPRC","PPV","NPV")] <-
    c(auc_val, sens, spec, acc, f1, brier_cal, brier_null, brier_rr, mae_cal, rmse_cal, pr_auc, ppv, npv)
  
  results[res_row, c("AUC_LCI","AUC_UCI")]                 <- ci_auc
  results[res_row, c("Sensitivity_LCI","Sensitivity_UCI")] <- ci_sens
  results[res_row, c("Specificity_LCI","Specificity_UCI")] <- ci_spec
  results[res_row, c("F1_LCI","F1_UCI")]                   <- ci_f1
  results[res_row, c("PPV_LCI","PPV_UCI")]                 <- ci_ppv
  results[res_row, c("NPV_LCI","NPV_UCI")]                 <- ci_npv
  results[res_row, c("Accuracy_LCI","Accuracy_UCI")]       <- ci_acc
  results[res_row, c("Brier_LCI","Brier_UCI")]             <- ci_brier
  results[res_row, c("AUPRC_LCI","AUPRC_UCI")]             <- ci_auprc
}

print(results)

##############################################
## Save results
##############################################
dir.create("Results", showWarnings = FALSE, recursive = TRUE)
writexl::write_xlsx(results, path = "Results/model_performance_by_week.xlsx")

## ---------------------------------------------------------------------------
## Importance ): compute on Week 0 / Week 4 if stored
## ---------------------------------------------------------------------------
avg_importance <- function(models, feat_names) {
  imps <- lapply(models, function(m){
    imp <- xgb.importance(feature_names = feat_names, model = m)
    imp <- imp %>% dplyr::select(Feature, Gain, Cover, Frequency)
    dplyr::full_join(tibble::tibble(Feature = feat_names), imp, by = "Feature") %>%
      dplyr::mutate(dplyr::across(c(Gain, Cover, Frequency), ~tidyr::replace_na(.x, 0)))
  }) %>% dplyr::bind_rows(.id = "model_id")
  
  imps %>%
    dplyr::group_by(Feature) %>%
    dplyr::summarise(
      Gain = mean(Gain), Cover = mean(Cover), Frequency = mean(Frequency),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Gain_Pct = 100 * Gain / ifelse(sum(Gain) > 0, sum(Gain), 1),
      Rank = dplyr::dense_rank(dplyr::desc(Gain))
    ) %>%
    dplyr::arrange(dplyr::desc(Gain)) %>%
    dplyr::select(Rank, Feature, Gain, Gain_Pct, Cover, Frequency)
}

# Print snapshots at weeks of interest
for (wk in c(0, 4)) {
  art <- wk_artifacts[[as.character(wk)]]
  if (!is.null(art) && !is.null(art$models)) {
    cat("\nVariable importance snapshot at Week", wk, ":\n")
    vi <- avg_importance(art$models, art$model_cols)
    print(utils::head(vi, 20), row.names = FALSE)
    # optional: save to CSV
    utils::write.csv(vi, file = paste0("VI_week_", wk, ".csv"), row.names = FALSE)
  }
}

##### Support ##################################################################
make_calib_df <- function(prob, y, bins = 5, min_bin_n = 1) {
  stopifnot(length(prob) == length(y))
  df <- tibble::tibble(prob = prob, y = y) |>
    dplyr::filter(!is.na(prob), !is.na(y)) |>
    dplyr::mutate(
      prob = pmin(pmax(prob, 1e-7), 1 - 1e-7),
      bin  = cut(prob, breaks = seq(0, 1, length.out = bins + 1),
                 include.lowest = TRUE, right = FALSE)
    ) |>
    dplyr::group_by(bin) |>
    dplyr::summarise(
      n         = dplyr::n(),
      mean_pred = mean(prob),
      observed  = mean(y == 1),
      .groups   = "drop"
    ) |>
    dplyr::arrange(mean_pred)
  
  df <- dplyr::filter(df, n >= min_bin_n)
  df
}

get_subgroup_metrics <- function(preds, truth, subgroup, thr = 0.5) {
  sub_data <- data.frame(pred = preds, truth = truth, group = subgroup)
  sub_data |>
    dplyr::group_by(group) |>
    dplyr::summarise(
      n = dplyr::n(),
      AUC = {
        y <- truth; p <- pred
        if (length(unique(y)) == 2) as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE))) else NA_real_
      },
      Brier = mean((pred - truth)^2),
      Sensitivity = {
        y <- truth; c <- as.integer(pred >= thr)
        tp <- sum(c==1 & y==1); fn <- sum(c==0 & y==1)
        if ((tp+fn)>0) tp/(tp+fn) else NA_real_
      },
      Specificity = {
        y <- truth; c <- as.integer(pred >= thr)
        tn <- sum(c==0 & y==0); fp <- sum(c==1 & y==0)
        if ((tn+fp)>0) tn/(tn+fp) else NA_real_
      },
      .groups = "drop"
    )
}

# Robust threshold chooser with guardrails
choose_threshold_robust <- function(
    probs, truth,
    sens_target = 0.70,
    min_spec    = 0.60,
    posrate_bounds = c(0.05, 0.95),
    grid = seq(0.01, 0.99, by = 0.01)
){
  stopifnot(length(probs) == length(truth))
  truth <- as.integer(truth)
  
  m <- vapply(grid, function(t){
    pred <- as.integer(probs >= t)
    tp <- sum(pred==1 & truth==1)
    fp <- sum(pred==1 & truth==0)
    fn <- sum(pred==0 & truth==1)
    tn <- sum(pred==0 & truth==0)
    
    sens <- if ((tp+fn)>0) tp/(tp+fn) else NA_real_
    spec <- if ((tn+fp)>0) tn/(tn+fp) else NA_real_
    ppv  <- if ((tp+fp)>0) tp/(tp+fp) else NA_real_
    f1   <- if (is.finite(ppv) && is.finite(sens) && (ppv+sens)>0) 2*ppv*sens/(ppv+sens) else NA_real_
    posr <- mean(pred==1)
    
    c(sens, spec, ppv, f1, posr)
  }, numeric(5))
  
  df <- data.frame(
    threshold = grid,
    Sensitivity = m[1,], Specificity = m[2,],
    PPV = m[3,], F1 = m[4,], PosRate = m[5,]
  )
  
  ok <- which(
    !is.na(df$F1) &
      df$Sensitivity >= sens_target &
      df$Specificity >= min_spec &
      df$PosRate >= posrate_bounds[1] &
      df$PosRate <= posrate_bounds[2]
  )
  
  pick <- if (length(ok)) {
    ok[which.max(df$F1[ok])]
  } else {
    ok2 <- which(!is.na(df$F1) &
                   df$Specificity >= min_spec &
                   df$PosRate >= posrate_bounds[1] &
                   df$PosRate <= posrate_bounds[2])
    if (length(ok2)) {
      ok2[which.max(df$F1[ok2])]
    } else {
      J <- with(df, Sensitivity + Specificity - 1)
      J[is.na(J)] <- -Inf
      cand <- which(J == max(J))
      cand_in_bounds <- cand[df$PosRate[cand] >= posrate_bounds[1] & df$PosRate[cand] <= posrate_bounds[2]]
      if (length(cand_in_bounds)) {
        cand <- cand_in_bounds[which.min(abs(df$Sensitivity[cand_in_bounds] - sens_target))]
      } else {
        cand <- cand[which.min(abs(df$Sensitivity[cand] - sens_target))]
      }
      cand[1]
    }
  }
  
  list(
    threshold = df$threshold[pick],
    metrics   = df[pick, ],
    curve     = df
  )
}

  ##############################################
  ## Calibration data 
  ##############################################
  if (!is.null(preds_calibrated_wk4)) {
    calib_data_wk4 <- make_calib_df(preds_calibrated_wk4, y_valid_wk4, bins = 5)
  } else warning("Week 4 calibration inputs are NULL.")
  
  theme_nejm <- function() {
    theme_classic(base_family = "Helvetica") +
      theme(
        plot.title = element_blank(),
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 14),
        axis.line  = element_line(size = 0.6, color = "black")
      )
  }
  
  # Plot Week 4 calibration (overall)
  if (exists("calib_data_wk4")) {
    calib_plot_wk4 <- ggplot(calib_data_wk4, aes(x = mean_pred, y = observed)) +
      geom_line(linewidth = 0.7, color = "black") +
      geom_point(size = 2.5, shape = 16, color = "black") +
      geom_abline(slope = 1, intercept = 0,
                  linetype = "dashed", color = "black", linewidth = 0.6) +
      labs(x = "Predicted Probability", y = "Observed Proportion") +
      theme_nejm()
    print(calib_plot_wk4)
  }
  
  tune_threshold_constrained <- function(prob, y, min_sens = 0.6, min_spec = 0.66, fallback = best_thres) {
    stopifnot(length(prob) == length(y))
    ok <- is.finite(prob) & !is.na(y)
    prob <- prob[ok]; y <- as.integer(y[ok])
    if (length(unique(y)) < 2L) return(fallback)
    
    roc_obj <- pROC::roc(response = y, predictor = prob, quiet = TRUE, direction = "<")
    coords  <- pROC::coords(roc_obj, x = roc_obj$thresholds, input = "threshold",
                            ret = c("threshold","sensitivity","specificity"), transpose = FALSE)
    df <- as.data.frame(coords)
    df$sens <- as.numeric(df$sensitivity)
    df$spec <- as.numeric(df$specificity)
    cand <- subset(df, sens >= min_sens & spec >= min_spec)
    if (!nrow(cand)) return(fallback)
    
    j <- cand$sens + cand$spec - 1
    cand$threshold[which.max(j)]
  }
  
  ##############################################
  ## Fairness assessment
  ##############################################
  
  if (!is.null(week4_valid)) {
    week4_valid_a <- week4_valid %>%
      mutate(
        Sex = dplyr::case_when(
          Gender == 0 ~ "Male",
          Gender == 1 ~ "Female",
          TRUE ~ NA_character_
        ),
        AgeGroup = ifelse(Age <= 35, "≤35", ">35")
      )
    
    # --- helper: compute metrics at a threshold ---
    .conf_metrics <- function(prob, y, thr = 0.5) {
      stopifnot(length(prob) == length(y))
      if (length(unique(y)) < 2) {
        return(tibble::tibble(
          Threshold = thr, Sensitivity = NA_real_, Specificity = NA_real_,
          PPV = NA_real_, NPV = NA_real_, F1 = NA_real_, Accuracy = NA_real_
        ))
      }
      pred <- as.integer(prob >= thr)
      tp <- sum(pred == 1 & y == 1)
      fp <- sum(pred == 1 & y == 0)
      tn <- sum(pred == 0 & y == 0)
      fn <- sum(pred == 0 & y == 1)
      sens <- ifelse((tp + fn) > 0, tp/(tp + fn), NA_real_)
      spec <- ifelse((tn + fp) > 0, tn/(tn + fp), NA_real_)
      ppv  <- ifelse((tp + fp) > 0, tp/(tp + fp), NA_real_)
      npv  <- ifelse((tn + fn) > 0, tn/(tn + fn), NA_real_)
      f1   <- ifelse((2*tp + fp + fn) > 0, 2*tp/(2*tp + fp + fn), NA_real_)
      acc  <- (tp + tn)/length(y)
      tibble::tibble(
        Threshold = thr, Sensitivity = sens, Specificity = spec,
        PPV = ppv, NPV = npv, F1 = f1, Accuracy = acc
      )
    }
    
    # --- helper: tune threshold within a subgroup (default: maximize Youden's J) ---
    tune_threshold <- function(prob, y, grid = seq(0.05, 0.95, by = 0.01),
                               target = c("youden", "f1", "balanced_acc")) {
      target <- match.arg(target)
      if (length(unique(y)) < 2) return(0.5)  # fallback if only one class present
      eval <- vapply(grid, function(t) {
        m <- .conf_metrics(prob, y, thr = t)
        switch(target,
               youden = (m$Sensitivity + m$Specificity - 1),
               f1 = m$F1,
               balanced_acc = mean(c(m$Sensitivity, m$Specificity), na.rm = TRUE)
        )
      }, numeric(1))
      grid[which.max(replace(eval, is.na(eval), -Inf))]
    }
    
    # Data aliases
    prob <- preds_calibrated_wk4
    y    <- y_valid_wk4
    sex  <- week4_valid_a$Sex
    
    # Indices per sex
    idx_f <- which(sex == "Female" & !is.na(prob) & !is.na(y))
    idx_m <- which(sex == "Male"   & !is.na(prob) & !is.na(y))
    
    # Learn per-sex thresholds (change target="f1" or "balanced_acc" if preferred)
    thr_female <- if (length(idx_f) >= 5) tune_threshold(prob[idx_f], y[idx_f], target = "youden") else best_thresh
    thr_male   <- if (length(idx_m) >= 5) tune_threshold(prob[idx_m], y[idx_m], target = "youden") else best_thresh
    thr_female <- if (length(idx_f) >= 5) tune_threshold_constrained(prob[idx_f], y[idx_f]) else 0.5
    thr_male   <- if (length(idx_m) >= 5) tune_threshold_constrained(prob[idx_m], y[idx_m]) else 0.5
    
    # Compute subgroup metrics at their own thresholds
    fair_female <- .conf_metrics(prob[idx_f], y[idx_f], thr = thr_female) %>%
      dplyr::mutate(group = "Female", n = length(idx_f))
    fair_male   <- .conf_metrics(prob[idx_m], y[idx_m], thr = thr_male) %>%
      dplyr::mutate(group = "Male", n = length(idx_m))
    
    sex_fairness_split <- dplyr::bind_rows(fair_female, fair_male) %>%
      dplyr::select(group, n, Threshold, Sensitivity, Specificity, PPV, NPV, F1, Accuracy)
    
    # (Optional) keep your original "single-threshold" results for comparison
    thr_global <- if (exists("best_thresh") && is.finite(best_thresh)) best_thresh else 0.5
    sex_fairness_global <- get_subgroup_metrics(prob, y, sex, thr = thr_global)
    
    # Age fairness can stay as-is
    age_fairness <- get_subgroup_metrics(prob, y, week4_valid_a$AgeGroup,
                                         thr = thr_global)
    
    cat("Per-sex learned thresholds:\n")
    print(tibble::tibble(Sex = c("Female","Male"),
                         Threshold = c(thr_female, thr_male)))
    
    cat("\nFairness by Sex (each at its own threshold):\n")
    print(sex_fairness_split)
    
    cat("\nFairness by Sex (single global threshold = ", signif(thr_global, 4), "):\n", sep = "")
    print(sex_fairness_global)
    
    cat("\nFairness by Age Group (global threshold):\n")
    print(age_fairness)
  }
  
  ##############################################
  ## Add 95% bootstrap CIs to fairness metrics
  ##############################################
  # Unified bootstrap for thresholded metrics + AUC
  boot_ci_all <- function(prob, y, thr, R = 1000, seed = 2025,
                          stratified = TRUE, conf = 0.95) {
    stopifnot(length(prob) == length(y))
    ok <- !is.na(prob) & !is.na(y)
    prob <- prob[ok]; y <- y[ok]
    n <- length(y)
    
    nm_th <- c("Sensitivity","Specificity","PPV","NPV","F1","Accuracy")
    
    # Point estimates for thresholded metrics (at fixed thr)
    pt_th <- .conf_metrics(prob, y, thr)
    est_th <- if (n > 0) as.numeric(pt_th[1, nm_th]) else rep(NA_real_, length(nm_th))
    
    # Point estimate for AUC (needs both classes)
    has_two_classes <- length(unique(y)) >= 2L
    auc_est <- if (has_two_classes) {
      as.numeric(pROC::auc(pROC::roc(response = y, predictor = prob, quiet = TRUE, direction = "<")))
    } else NA_real_
    
    # Point estimate for calibrated Brier (Platt; fallback to raw probs on failure)
    brier_cal_est <- tryCatch({
      p_cal <- stats::predict(stats::glm(y ~ prob, family = binomial()), type = "response")
      mean((p_cal - y)^2)
    }, error = function(e) {
      mean((prob - y)^2)
    })
    
    set.seed(seed)
    
    # Precompute class indices for stratified bootstrap
    idx0 <- which(y == 0)
    idx1 <- which(y == 1)
    draw_idx <- function() {
      if (stratified && length(idx0) > 0 && length(idx1) > 0) {
        c(sample(idx0, length(idx0), replace = TRUE),
          sample(idx1, length(idx1), replace = TRUE))
      } else {
        sample(seq_len(n), n, replace = TRUE)
      }
    }
    
    mats <- replicate(R, {
      idb <- draw_idx()
      yb  <- y[idb]; pb <- prob[idb]
      
      # Thresholded metrics (need both classes for sens/spec/etc.)
      if (length(unique(yb)) >= 2L) {
        mb <- .conf_metrics(pb, yb, thr)
        th_vals <- as.numeric(mb[1, nm_th])
        auc_b <- as.numeric(pROC::auc(pROC::roc(yb, pb, quiet = TRUE, direction = "<")))
      } else {
        th_vals <- rep(NA_real_, length(nm_th))
        auc_b   <- NA_real_
      }
      
      c(th_vals, auc_b)   # length = 6 + 1
    })
    
    brier_vec <- vapply(seq_len(R), function(i) {
      idb <- draw_idx()
      yb  <- y[idb]; pb <- prob[idb]
      p_cal_b <- tryCatch({
        if (length(unique(yb)) >= 2L) {
          stats::predict(stats::glm(yb ~ pb, family = binomial()), type = "response")
        } else {
          pb  # fallback: no calibration possible if only one class in resample
        }
      }, error = function(e) pb)
      mean((p_cal_b - yb)^2)
    }, numeric(1))
    
    # Percentile CIs
    alpha <- (1 - conf)/2
    qs <- apply(mats, 1L, stats::quantile, probs = c(alpha, 1 - alpha), na.rm = TRUE)
    qs_brier <- stats::quantile(brier_vec, probs = c(alpha, 1 - alpha), na.rm = TRUE)
    
    out_th <- tibble::tibble(
      Metric    = nm_th,
      Est       = est_th,
      LCI       = qs[1, seq_along(nm_th)],
      UCI       = qs[2, seq_along(nm_th)],
      Threshold = pt_th$Threshold
    )
    
    out_auc <- tibble::tibble(
      Metric    = "AUC",
      Est       = auc_est,
      LCI       = qs[1, length(nm_th) + 1L],
      UCI       = qs[2, length(nm_th) + 1L],
      Threshold = NA_real_
    )
    
    out_brier <- tibble::tibble(
      Metric    = "Brier_Cal",
      Est       = brier_cal_est,
      LCI       = qs_brier[1],
      UCI       = qs_brier[2],
      Threshold = NA_real_
    )
    
    dplyr::bind_rows(out_th, out_auc, out_brier)
  }
  
  if (!is.null(week4_valid)) {
    week4_valid_a <- week4_valid %>%
      dplyr::mutate(
        Sex = dplyr::case_when(
          Gender == 0 ~ "Male",
          Gender == 1 ~ "Female",
          TRUE ~ NA_character_
        ),
        AgeGroup = ifelse(Age <= 35, "≤35", ">35")
      )
    
    # aliases
    prob <- preds_calibrated_wk4
    y    <- y_valid_wk4
    sex  <- week4_valid_a$Sex
    
    # indices
    idx_f <- which(sex == "Female" & !is.na(prob) & !is.na(y))
    idx_m <- which(sex == "Male"   & !is.na(prob) & !is.na(y))
    
    # learn thresholds per sex
    thr_female <- if (length(idx_f) >= 5) tune_threshold_constrained(prob[idx_f], y[idx_f]) else best_thresh
    thr_male   <- if (length(idx_m) >= 5) tune_threshold_constrained(prob[idx_m], y[idx_m]) else best_thresh 
    
    fair_female <- .conf_metrics(prob[idx_f], y[idx_f], thr = thr_female) %>%
      dplyr::mutate(group = "Female", n = length(idx_f))
    fair_male   <- .conf_metrics(prob[idx_m], y[idx_m], thr = thr_male) %>%
      dplyr::mutate(group = "Male", n = length(idx_m))
    sex_fairness_split <- dplyr::bind_rows(fair_female, fair_male) %>%
      dplyr::select(group, n, Threshold, Sensitivity, Specificity, PPV, NPV, F1, Accuracy)
    
    female_ci_all <- boot_ci_all(prob[idx_f], y[idx_f], thr = thr_female, R = 2000) %>%
      dplyr::mutate(group = "Female", n = length(idx_f)) %>%
      dplyr::relocate(group, n, Threshold)
    male_ci_all   <- boot_ci_all(prob[idx_m], y[idx_m], thr = thr_male,   R = 2000) %>%
      dplyr::mutate(group = "Male",   n = length(idx_m)) %>%
      dplyr::relocate(group, n, Threshold)
    
    sex_fairness_split_ci <- dplyr::bind_rows(female_ci_all, male_ci_all)
    
    thr_global <- if (exists("best_thresh") && is.finite(best_thresh)) best_thresh else 0.5
    
    # Sex (global thr)
    sex_global_ci <- dplyr::bind_rows(
      boot_ci_all(prob[idx_f], y[idx_f], thr = thr_global, R = 2000) %>%
        dplyr::mutate(group = "Female", n = length(idx_f)),
      boot_ci_all(prob[idx_m], y[idx_m], thr = thr_global, R = 2000) %>%
        dplyr::mutate(group = "Male",   n = length(idx_m))
    ) %>% dplyr::relocate(group, n, Threshold)
    
    # Age groups (global thr)
    age <- week4_valid_a$AgeGroup
    idx_young <- which(age == "≤35" & !is.na(prob) & !is.na(y))
    idx_old   <- which(age == ">35" & !is.na(prob) & !is.na(y))
    
    age_fairness_ci <- dplyr::bind_rows(
      boot_ci_all(prob[idx_young], y[idx_young], thr = thr_global, R = 2000) %>%
        dplyr::mutate(group = "≤35", n = length(idx_young)),
      boot_ci_all(prob[idx_old],   y[idx_old],   thr = thr_global, R = 2000) %>%
        dplyr::mutate(group = ">35", n = length(idx_old))
    ) %>% dplyr::relocate(group, n, Threshold)
    
    # ---- Printouts ----
    cat("Per-sex learned thresholds:\n")
    print(tibble::tibble(Sex = c("Female","Male"),
                         Threshold = c(thr_female, thr_male)))
    
    cat("\nFairness by Sex (each at its own threshold) — point estimates:\n")
    print(sex_fairness_split)
    
    cat("\nFairness by Sex (each at its own threshold) — 95% CIs (incl. AUC):\n")
    print(sex_fairness_split_ci)
    
    cat("\nFairness by Sex (single global threshold = ", signif(thr_global, 4), ") — 95% CIs (incl. AUC):\n", sep = "")
    print(sex_global_ci)
  }
  
  ##############################################
  ## Fairness by Age Group with Per-Group Thresholds
  ##############################################
  if (!is.null(week4_valid)) {
    # Build age groups (≤35 vs >35)
    week4_valid_a <- week4_valid %>%
      dplyr::mutate(AgeGroup = ifelse(Age <= 35, "≤35", ">35"))
    
    # Aliases
    prob <- preds_calibrated_wk4   # calibrated preds you saved earlier
    y    <- y_valid_wk4
    age  <- week4_valid_a$AgeGroup
    
    # Indices per age group
    idx_young <- which(age == "≤35" & !is.na(prob) & !is.na(y))
    idx_old   <- which(age == ">35" & !is.na(prob) & !is.na(y))
    
    # Learn per-group thresholds using constrained Youden (sens/spec ≥ 0.70)
    thr_young <- if (length(idx_young) >= 5) {
      tune_threshold_constrained(prob[idx_young], y[idx_young], min_sens = 0.70, min_spec = 0.70, fallback = 0.5)
    } else 0.5
    
    thr_old <- if (length(idx_old) >= 5) {
      tune_threshold_constrained(prob[idx_old], y[idx_old], min_sens = 0.70, min_spec = 0.70, fallback = 0.5)
    } else 0.5
    
    # -------- Point estimates at each group's own threshold --------
    fair_young_pt <- .conf_metrics(prob[idx_young], y[idx_young], thr = thr_young) %>%
      dplyr::mutate(group = "≤35", n = length(idx_young), Threshold = thr_young) %>%
      dplyr::relocate(group, n, Threshold)
    
    fair_old_pt <- .conf_metrics(prob[idx_old], y[idx_old], thr = thr_old) %>%
      dplyr::mutate(group = ">35", n = length(idx_old), Threshold = thr_old) %>%
      dplyr::relocate(group, n, Threshold)
    
    age_fairness_split <- dplyr::bind_rows(fair_young_pt, fair_old_pt) %>%
      dplyr::select(group, n, Threshold, Sensitivity, Specificity, PPV, NPV, F1, Accuracy)
    
    cat("\nFairness by Age Group (each at its own constrained-Youden threshold) — point estimates:\n")
    print(age_fairness_split)
    
    # -------- 95% bootstrap CIs for each group's own threshold --------
    set.seed(2025)
    young_ci_all <- boot_ci_all(prob[idx_young], y[idx_young], thr = thr_young, R = 2000) %>%
      dplyr::mutate(group = "≤35", n = length(idx_young)) %>%
      dplyr::relocate(group, n, Threshold)
    
    old_ci_all <- boot_ci_all(prob[idx_old], y[idx_old], thr = thr_old, R = 2000) %>%
      dplyr::mutate(group = ">35", n = length(idx_old)) %>%
      dplyr::relocate(group, n, Threshold)
    
    age_fairness_split_ci <- dplyr::bind_rows(
      young_ci_all %>% dplyr::mutate(Threshold = ifelse(Metric %in% c("AUC","Brier_Cal"), NA_real_, thr_young)),
      old_ci_all   %>% dplyr::mutate(Threshold = ifelse(Metric %in% c("AUC","Brier_Cal"), NA_real_, thr_old))
    )
    
    cat("\nFairness by Age Group (each at its own constrained-Youden threshold) — 95% CIs (incl. AUC & Brier_Cal):\n")
    print(age_fairness_split_ci)
  }
  
  
  ##############################################
  ## ROC Curves for Weeks 4 
  ##############################################
  # Week 4 ROC using aligned labels + ensemble preds
  if (exists("y_valid_wk4_aligned") && exists("preds_prob_wk4_ens")) {
    roc_4 <- pROC::roc(response = y_valid_wk4_aligned,
                       predictor = preds_prob_wk4_ens,
                       quiet = TRUE)
  }
  print(roc_4)
  
  # Plot both
 
  if (exists("roc_4")) {
    # Convert to percent coords
    x4  <- (1 - roc_4$specificities)  * 100
    y4  <- (roc_4$sensitivities)      * 100
   
    # pad & sort so polylines draw correctly and reach the corners
    pad_sort <- function(x, y) {
      x2 <- c(0, x, 100); y2 <- c(0, y, 100)
      o  <- order(x2, y2)
      list(x = x2[o], y = y2[o])
    }
    p4  <- pad_sort(x4, y4)
    
    # Panel
    par(family="Helvetica", mar=c(5,5,2,2),
        lwd=2, cex.lab=1.6, cex.axis=1.4,
        xaxs="i", yaxs="i", lend="butt")  # 'butt' removes rounded line ends
    
    plot(NA, NA, type="n",
         xlim=c(0,100), ylim=c(0,100), asp=1, axes=FALSE,
         xlab="1 - Specificity", ylab="Sensitivity")
    
    # Draw axes exactly at origin + ticks
    segments(0,0, 100,0, lwd=2)  # x-axis line
    segments(0,0, 0,100, lwd=2)  # y-axis line
    axis(1, at=seq(0,100,20), labels=seq(0,100,20),
         pos=0, lwd=0, lwd.ticks=2, cex.axis=1.4)
    axis(2, at=seq(0,100,20), labels=seq(0,100,20),
         pos=0, lwd=0, lwd.ticks=2, las=1, cex.axis=1.4)
    
    # 45° reference
    abline(a=0, b=1, lty=2, col="black", lwd=2)
    
    # Curves: lines + dots
    lines(p4$x,  p4$y,  col="#1f78b4", lwd=3)
    
    # Legend
    legend("bottomright", inset=c(0,0.03),
           legend=c(sprintf("Week 4   AUC = %.3f", as.numeric(pROC::auc(roc_4)))),
           col=c("#1f78b4"), lwd=4, bty="n", cex=1.4)
  }
  
  
  ##############################################
  ## 9) Combined ggplot panels (optional)
  ##############################################
  results_long <- results %>%
    dplyr::select(
      WeeksUsed,
      Sensitivity, Sensitivity_LCI, Sensitivity_UCI,
      Specificity, Specificity_LCI, Specificity_UCI,
      F1_Score, F1_LCI, F1_UCI
    ) %>%
    tidyr::pivot_longer(cols = -WeeksUsed, names_to = "Metric", values_to = "Value") %>%
    dplyr::mutate(
      MetricType = dplyr::case_when(
        grepl("_LCI$", Metric) ~ "LCI",
        grepl("_UCI$", Metric) ~ "UCI",
        TRUE ~ "Main"
      ),
      Metric = gsub("_LCI|_UCI", "", Metric)
    ) %>%
    tidyr::pivot_wider(names_from = MetricType, values_from = Value)
  
  # ------ Find peak F1 in the first 10 weeks (earliest if tie) ------
  sub10 <- results %>%
    dplyr::filter(WeeksUsed >= 1, WeeksUsed <= 10, is.finite(F1_Score))
  
  if (nrow(sub10)) {
    max_f1 <- max(sub10$F1_Score, na.rm = TRUE)
    # earliest week among tied maxima
    peak_week <- min(sub10$WeeksUsed[sub10$F1_Score == max_f1])
    peak_f1   <- round(results$F1_Score[match(peak_week, results$WeeksUsed)], 3)
  } else {
    # fallback: overall peak (earliest if tie)
    max_f1 <- max(results$F1_Score, na.rm = TRUE)
    peak_week <- min(results$WeeksUsed[results$F1_Score == max_f1])
    peak_f1   <- round(results$F1_Score[match(peak_week, results$WeeksUsed)], 3)
  }
  
  # ------ Plot Sensitivity, Specificity, F1; mark peak week ------
  results_long_main <- results_long %>%
    dplyr::filter(Metric %in% c("Sensitivity", "Specificity", "F1_Score"))
  
  x_min <- min(results$WeeksUsed, na.rm = TRUE)
  x_max <- max(results$WeeksUsed, na.rm = TRUE)
  x_breaks <- seq(x_min, x_max, by = 1) 
  
  annotation_label <- paste0("Peak F1\nWeek ", peak_week, " (", peak_f1, ")")
  
  p2 <- ggplot(results_long_main, aes(x = WeeksUsed, y = Main, color = Metric, group = Metric)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_vline(xintercept = peak_week, linetype = "dashed", color = "black") +
    annotate("text",
             x = pmin(peak_week + 0.3, x_max - 0.3),  # nudge right, avoid clipping
             y = 0.2, label = annotation_label,
             hjust = 0, vjust = 1, size = 4.5, family = "Helvetica") +
    scale_color_manual(values = c("Sensitivity" = "darkgreen",
                                  "Specificity" = "firebrick",
                                  "F1_Score"   = "darkorange")) +
    scale_x_continuous(limits = c(x_min, x_max +0.5),
                       breaks = x_breaks,
                       minor_breaks = NULL,
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1),
                       breaks = seq(0, 1, 0.1),
                       expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    labs(x = "Weeks of Adherence Data Used", y = "Model Performance Metric") +
    theme_classic(base_family = "Helvetica") +
    theme(
      axis.title = element_text(size = 16, color = "black", family = "Helvetica"),
      axis.text  = element_text(size = 16, color = "black", family = "Helvetica"),
      axis.line  = element_line(color = "black"),
      legend.title = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 14, color = "black", family = "Helvetica"),
      plot.title = element_blank()
    )
  
  print(p2)
  
  x_min <- min(results$WeeksUsed, na.rm = TRUE)
  x_max <- max(results$WeeksUsed, na.rm = TRUE)
  x_breaks <- seq(x_min, x_max, by = 1)   # use by = 2 for fewer ticks
  
  p1 <- ggplot(results, aes(x = WeeksUsed, y = AUC)) +
    geom_line(color = "steelblue", linewidth = 1) +
    geom_point(color = "steelblue", size = 2, shape = 16) +
    geom_errorbar(aes(ymin = AUC_LCI, ymax = AUC_UCI),
                  width = 0.2, color = "steelblue", alpha = 0.7) +
    geom_hline(yintercept = 0.80, linetype = "dashed", color = "gray40") +
    scale_x_continuous(limits = c(x_min, x_max +0.5),
                       breaks = x_breaks, minor_breaks = NULL,
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1),
                       breaks = seq(0, 1, 0.1),
                       expand = c(0, 0)) +
    labs(x = "Weeks of Adherence Data Used",
         y = "AUROC") +
    theme_classic(base_family = "Helvetica") +
    theme(
      axis.title = element_text(size = 16, color = "black", family = "Helvetica"),
      axis.text  = element_text(size = 16, color = "black", family = "Helvetica"),
      axis.line  = element_line(color = "black"),
      legend.position = "none",
      plot.title = element_blank(),
      plot.margin = margin(t = 15, r = 15, b = 15, l = 15)
    )
  
  print(p1)
  
  #Combined plots 
  library(cowplot)
  pad_right <- 10
  pad_left <- 20
  p1_pad <- p1 + theme(plot.margin = margin(t = 40, r = pad_right, b = 0, l = pad_left))
  p2_pad <- p2 + theme(plot.margin = margin(t = 10, r = pad_right, b = 5, l = pad_left))
  
  combined <- plot_grid(
    p1_pad, p2_pad,
    ncol = 1, align = "v",
    labels = c("A", "B"),
    label_fontfamily = "Helvetica",
    label_fontface = "plain",
    label_size = 16,
    label_x = 0.01,
    label_y = c(1, 1),   # keep at the very top
    hjust = 0, vjust = 2.5
  )
  
  print(combined)

  
  

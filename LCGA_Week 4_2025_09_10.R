# =============================================================================
# LCGA of Bedaquiline Adherence (Weeks 1–4 only)
#   - Fits flexmix mixture of linear trajectories: BDQ ~ Week | ID
#   - Evaluates K = 1..5 with multiple random starts
#   - Outputs: fit metrics, class assignments, class sizes, and diagnostic plots
# Author: Kevin J. Guzman
# Updated: 2025-09-11
# =============================================================================

suppressPackageStartupMessages({
  library(flexmix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(tibble)
  library(stringr)
})

set.seed(20250911)

# -----------------------------------------------------------------------------
# I/O
# -----------------------------------------------------------------------------
# Input files (adjust paths if needed)
in_bdq <- "Files/BDQ_Week.RData"
in_prx <- "Files/PRX3.RData"

# Output directory
out_dir <- "outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
ent_entropy <- function(P) {
  # Normalized entropy (higher = better separation), computed on ID-level posteriors
  # P: data.frame or matrix with columns p1..pK for each ID (rows)
  stopifnot(is.matrix(P) || is.data.frame(P))
  P <- as.matrix(P)
  K <- ncol(P); n <- nrow(P)
  if (K < 1 || n < 1) return(NA_real_)
  # add small epsilon to avoid log(0)
  eps <- 1e-12
  1 + sum(P * log(P + eps)) / (n * log(K))
}

safe_flexmix <- function(formula, data, k, nrep = 50L, verbose = FALSE) {
  # Wrapper to run flexmix robustly
  tryCatch(
    flexmix(formula, data = data, k = k,
            model = FLXMRglm(family = "gaussian"),
            control = list(nrep = nrep, verb = if (verbose) 1 else 0)),
    error = function(e) {
      message(sprintf("   [warn] flexmix failed for K=%s: %s", k, e$message))
      return(NULL)
    }
  )
}

# -----------------------------------------------------------------------------
# Load & prepare data
# -----------------------------------------------------------------------------
load(in_bdq)  # expects object: BDQ_Week
load(in_prx)  # expects object: PRX3

# Merge minimal demographics (optional) and drop transfer/other EOT outcome = 6
BDQ_Week <- BDQ_Week %>%
  left_join(PRX3 %>% select(ID, Age, Gender, BMI), by = "ID") %>%
  filter(EOT_Outcome != 6)

# Optional outcome (not used for LCGA itself; retained for downstream joins)
BDQ_Week <- BDQ_Week %>%
  mutate(
    Outcome = dplyr::case_when(
      EOT_Outcome %in% c(1, 2) ~ 1L,             # Successful
      EOT_Outcome %in% c(3, 4, 5) ~ 0L,          # Unsuccessful
      TRUE ~ NA_integer_
    )
  )

# Keep only Weeks 1–4
week_vars <- paste0("Week_", 1:4)
have_weeks <- intersect(week_vars, names(BDQ_Week))
if (length(have_weeks) < 2L) {
  stop("Need at least 2 week columns among Weeks 1–4 for LCGA.")
}

# Long format for LCGA
BDQ_long <- BDQ_Week %>%
  select(ID, all_of(have_weeks)) %>%
  pivot_longer(cols = starts_with("Week_"), names_to = "Week", values_to = "BDQ") %>%
  mutate(
    Week = as.numeric(str_remove(Week, "Week_")),
    BDQ  = as.numeric(BDQ)
  ) %>%
  filter(Week %in% 1:4) %>%
  filter(is.finite(BDQ)) %>%
  arrange(ID, Week)

# -----------------------------------------------------------------------------
# Fit LCGA on Weeks 1–4 only
# -----------------------------------------------------------------------------
K_grid      <- 1:5         
nrep_starts <- 50L         # Multiple random starts for stability

models   <- vector("list", length(K_grid))
names(models) <- as.character(K_grid)

fit_tbl <- tibble(
  K              = integer(),
  BIC            = numeric(),
  LogLik         = numeric(),
  Entropy        = numeric(),
  VLMR           = numeric(),
  VLMR_p         = numeric(),
  n_components   = integer()
)

assignments_all <- tibble()

message("Fitting LCGA using Weeks 1–4 (K = 1..5) ...")
for (k in K_grid) {
  message(sprintf(" - K = %d", k))
  model_k <- safe_flexmix(BDQ ~ Week | ID, data = BDQ_long, k = k, nrep = nrep_starts)
  if (is.null(model_k)) next
  
  models[[as.character(k)]] <- model_k
  
  # Basic metrics
  bic_k    <- BIC(model_k)
  ll_k     <- as.numeric(logLik(model_k))
  
  # Posterior probabilities per *row* (each ID-week obs)
  post_row <- posterior(model_k)  # n_obs x k matrix
  # Bind ID to rows, then average posteriors across repeated measures per ID
  post_id <- cbind(ID = BDQ_long$ID, as.data.frame(post_row)) %>%
    group_by(ID) %>%
    summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
    arrange(ID)
  
  # Compute normalized entropy at the ID level
  P <- as.matrix(post_id[ , setdiff(names(post_id), "ID"), drop = FALSE])
  ent_k <- ent_entropy(P)
  
  # Most likely class per ID & max posterior
  class_vec <- apply(P, 1, which.max)
  maxprob   <- apply(P, 1, max)
  
  assign_df <- tibble(
    ID           = post_id$ID,
    K            = k,
    Class        = class_vec,
    MaxPosterior = maxprob
  )
  assignments_all <- bind_rows(assignments_all, assign_df)
  
  # Simple VLMR-like statistic (approximate; df = 1 used as heuristic)
  vlmr   <- NA_real_
  vlmr_p <- NA_real_
  if (k > 1 && !is.null(models[[as.character(k - 1)]])) {
    ll_prev <- as.numeric(logLik(models[[as.character(k - 1)]]))
    if (is.finite(ll_prev)) {
      vlmr   <- 2 * (ll_k - ll_prev)
      vlmr_p <- 1 - pchisq(vlmr, df = 1)
    }
  }
  
  fit_tbl <- add_row(
    fit_tbl,
    K = k,
    BIC = bic_k,
    LogLik = ll_k,
    Entropy = ent_k,
    VLMR = vlmr,
    VLMR_p = vlmr_p,
    n_components = k
  )
}

# -----------------------------------------------------------------------------
# Summaries & exports
# -----------------------------------------------------------------------------
# Class sizes & proportions at each K
class_summary <- assignments_all %>%
  count(K, Class, name = "N") %>%
  group_by(K) %>%
  mutate(Proportion = N / sum(N)) %>%
  ungroup() %>%
  arrange(K, Class)

# Wide layout (optional)
class_summary_wide <- class_summary %>%
  pivot_longer(cols = c(N, Proportion),
               names_to = "Metric", values_to = "Value") %>%
  mutate(VarName = paste0("Class_", Class, "_", ifelse(Metric == "N", "N", "P"))) %>%
  select(K, VarName, Value) %>%
  pivot_wider(names_from = VarName, values_from = Value) %>%
  arrange(K)

# Choose “best” K by (i) highest entropy (>= 0.80 preferred), (ii) lowest BIC
valid_k <- class_summary %>%
  group_by(K) %>%
  summarise(min_prop = min(Proportion, na.rm = TRUE), .groups = "drop")

fit_tbl <- fit_tbl %>%
  left_join(valid_k, by = "K") %>%
  mutate(
    Entropy_flag = if_else(Entropy >= 0.80, 1L, 0L),
    Size_flag    = if_else(min_prop > 0.10, 1L, 0L)
  ) %>%
  arrange(desc(Entropy_flag), desc(Size_flag), BIC)

best_row <- fit_tbl %>% slice(1)
best_K   <- best_row$K

message(sprintf(
  "Selected K = %s (Entropy=%.3f, BIC=%.1f, Min class proportion=%.2f)",
  best_K, best_row$Entropy, best_row$BIC, best_row$min_prop
))

# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------
p_bic <- ggplot(fit_tbl, aes(x = K, y = BIC)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_vline(xintercept = best_K, linetype = "dashed") +
  labs(title = "BIC across K (Weeks 1–4)", x = "Number of Classes (K)", y = "BIC") +
  theme_minimal(base_size = 12)

p_entropy <- ggplot(fit_tbl, aes(x = K, y = Entropy)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.80, linetype = "dotted") +
  geom_vline(xintercept = best_K, linetype = "dashed") +
  labs(title = "Normalized Entropy across K (Weeks 1–4)",
       x = "Number of Classes (K)", y = "Entropy (ID-level)") +
  theme_minimal(base_size = 12)

p_classes <- class_summary %>%
  ggplot(aes(x = factor(K), y = Proportion, fill = factor(Class))) +
  geom_col(position = "stack") +
  labs(title = "Class Composition by K (Weeks 1–4)",
       x = "K", y = "Proportion", fill = "Class") +
  theme_minimal(base_size = 12)

# -----------------------------------------------------------------------------
# Best-K assignments merged back to baseline (optional)
# -----------------------------------------------------------------------------
best_assign <- assignments_all %>% filter(K == best_K)
  table(best_assign$Class)


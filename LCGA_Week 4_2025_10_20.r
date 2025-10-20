suppressPackageStartupMessages({
  library(flexmix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(tibble)
  library(stringr)
  library(tidyLPA)
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
## Entropy helper calculation 
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

## Safe flexmix wrapper
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


# Helper to get a posterior-probability matrix from a flexmix model
.safe_post <- function(m) {
   pp <- flexmix::posterior(m)
   if (is.list(pp)) {
      as.matrix(pp$posterior)
   } else {
      as.matrix(pp)
   }
}

# Helper to calculate boostrap log-likelihood ratio test for flexmix models
boot_LRT_flexmix <- function(m_small, m_large, data, form,
                             nboot = 200, nrep = 10, seed = 123,
                             control = list(iter.max = 200)) {
   stopifnot(inherits(m_small, "flexmix"), inherits(m_large, "flexmix"))
   set.seed(seed)
   
   # Observed LR
   LR_obs <- 2 * (as.numeric(logLik(m_large)) - as.numeric(logLik(m_small)))
   
   # Posterior probs from smaller model (for stratified resampling)
   P_small <- flexmix::posterior(m_small)
   if (is.list(P_small)) P_small <- as.matrix(P_small$posterior)
   small_class <- max.col(P_small, ties.method = "first")
   
   # Class proportions under H0 (smaller model)
   class_props <- prop.table(table(small_class))
   classes <- sort(as.integer(names(class_props)))
   n <- nrow(data)
   
   # Preindex rows by class
   rows_by_class <- lapply(classes, function(k) which(small_class == k))
   names(rows_by_class) <- as.character(classes)
   
   # k-values from models
   k_small <- m_small@k
   k_large <- m_large@k
   
   LR_boot <- numeric(nboot)
   
   for (b in seq_len(nboot)) {
      ## ---- Stratified resample ----
      draws <- pmax(0L, round(class_props * n))
      while (sum(draws) < n) draws[which.max(class_props)] <- draws[which.max(class_props)] + 1L
      while (sum(draws) > n) draws[which.max(draws)] <- draws[which.max(draws)] - 1L
      
      idx_sim <- unlist(mapply(function(k, m) {
         pool <- rows_by_class[[as.character(k)]]
         if (length(pool) == 0) sample(seq_len(n), m, replace = TRUE)
         else sample(pool, m, replace = TRUE)
      }, k = classes, m = as.integer(draws), SIMPLIFY = FALSE), use.names = FALSE)
      
      sim_data <- data[idx_sim, , drop = FALSE]
      
      ## ---- Refit both models ----
      m_small_b <- flexmix::flexmix(form, data = sim_data, k = k_small,
                                    control = c(control, list(nrep = nrep)))
      m_large_b <- flexmix::flexmix(form, data = sim_data, k = k_large,
                                    control = c(control, list(nrep = nrep)))
      
      LR_boot[b] <- 2 * (as.numeric(logLik(m_large_b)) - as.numeric(logLik(m_small_b)))
   }
   
   pval <- mean(LR_boot >= LR_obs)
   
   list(
      LR_obs = LR_obs,
      pval = pval,
      LR_boot = LR_boot,
      k_small = k_small,
      k_large = k_large,
      class_props = class_props
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

# Outcomes 
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
### Saving Results of LCGA Model #######
K_grid      <- 1:5         
nrep_starts <- 50L         # Multiple random starts for stability, here we use 50 

models   <- vector("list", length(K_grid))
names(models) <- as.character(K_grid)

fit_tbl <- tibble(
  K              = integer(),
  BIC            = numeric(),
  LogLik         = numeric(),
  Entropy        = numeric(),
  n_components   = integer()
)

assignments_all <- tibble()

### LCGA Model ############################
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
  
  fit_tbl <- add_row(
    fit_tbl,
    K = k,
    BIC = bic_k,
    LogLik = ll_k,
    Entropy = ent_k,
    n_components = k
  )
}

m1 <- flexmix(BDQ ~ Week | ID, data = BDQ_long, k = 1)
m2 <- flexmix(BDQ ~ Week | ID, data = BDQ_long, k = 2)
m3 <- flexmix(BDQ ~ Week | ID, data = BDQ_long, k = 3)
m4 <- flexmix(BDQ ~ Week | ID, data = BDQ_long, k = 4)
m5 <- flexmix(BDQ ~ Week | ID, data = BDQ_long, k = 5)
m6 <- flexmix(BDQ ~ Week | ID, data = BDQ_long, k = 6)

res12 <- boot_LRT_flexmix(
  m_small = m1, m_large = m2,
  data = BDQ_long,
  form = BDQ ~ Week | ID,
  nboot = 1000, nrep = 15, seed = 42
)


res23 <- boot_LRT_flexmix(
  m_small = m2, m_large = m3,
  data = BDQ_long,
  form = BDQ ~ Week | ID,
  nboot = 1000, nrep = 15, seed = 42
)


res34 <- boot_LRT_flexmix(
  m_small = m3, m_large = m4,
  data = BDQ_long,
  form = BDQ ~ Week | ID,
  nboot = 1000, nrep = 15, seed = 42
)

res45 <- boot_LRT_flexmix(
  m_small = m4, m_large = m5,
  data = BDQ_long,
  form = BDQ ~ Week | ID,
  nboot = 1000, nrep = 15, seed = 42
)

res12$LR_obs
res12$pval
res23$LR_obs
res23$pval
res34$LR_obs
res34$pval
res45$LR_obs
res45$pval

# LMR calcuatlion using flexmix for 1 v. 2 
n   <- 282
ll1 <- as.numeric(logLik(m1))
p1  <- m1@df
ll2 <- as.numeric(logLik(m2))
p2  <- m2@df

tidyLPA::calc_lrt(
  n = n,
  null_ll = ll1, null_param = p1, null_classes = 1,
  alt_ll  = ll2, alt_param  = p2, alt_classes  = 2
)

#LMR calcualtion using fleximis for 2 v 3 
n   <- 282
ll2 <- as.numeric(logLik(m2))
p2  <- m2@df
ll3 <- as.numeric(logLik(m3))
p3  <- m3@df
tidyLPA::calc_lrt(
  n = n,
  null_ll = ll2, null_param = p2, null_classes = 2,
  alt_ll  = ll3, alt_param  = p3, alt_classes  = 3
)

#LMR calculation using flexmix for 3 v 4
n   <- 282
ll3 <- as.numeric(logLik(m3))
p3  <- m3@df
ll4 <- as.numeric(logLik(m4))
p4  <- m4@df
tidyLPA::calc_lrt(
  n = n,
  null_ll = ll3, null_param = p3, null_classes = 3,
  alt_ll  = ll4, alt_param  = p4, alt_classes  =4)

#LMR calculation using flexmix for 4 v 5
n   <- 282
ll4 <- as.numeric(logLik(m4))
p4  <- m4@df
ll5 <- as.numeric(logLik(m5))
p5  <- m5@df
tidyLPA::calc_lrt(
  n = n,
  null_ll = ll4, null_param = p4, null_classes = 4,
  alt_ll  = ll5, alt_param  = p5, alt_classes  =5)


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


# =============================================================================
# Wine Quality (Red): LOO Additive + HSIC Causal Discovery (with Misoriented)
# =============================================================================

suppressPackageStartupMessages({
  need <- c("mgcv","ggplot2","dplyr","fastICA","gridExtra")
  miss <- need[!sapply(need, requireNamespace, quietly = TRUE)]
  if (length(miss)) {
    install.packages(miss, repos = "https://cloud.r-project.org")
  }
  library(mgcv); library(ggplot2); library(dplyr); library(fastICA); library(gridExtra)
})

set.seed(1)

# -----------------------------------------------------------------------------
# 0) IO setup
# -----------------------------------------------------------------------------
if (!dir.exists("data"))    dir.create("data", recursive = TRUE)
if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)
if (!dir.exists("results")) dir.create("results", recursive = TRUE)

# -----------------------------------------------------------------------------
# 1) Data (keep header names with spaces exactly as-is)
# -----------------------------------------------------------------------------
if (!file.exists("data/winequality-red.csv")) {
  cat("Downloading winequality-red.csv ...\n")
  download.file(
    "https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-red.csv",
    "data/winequality-red.csv", quiet = TRUE
  )
  cat("Download complete.\n")
}
wine_raw <- read.csv("data/winequality-red.csv", sep = ";", check.names = FALSE)
stopifnot(is.data.frame(wine_raw), nrow(wine_raw) > 0)

cat("Columns found:\n  ", paste(colnames(wine_raw), collapse = ", "), "\n")
wine <- as.data.frame(scale(wine_raw))
cat(sprintf("Wine Quality (red): %d samples x %d variables\n\n", nrow(wine), ncol(wine)))

# -----------------------------------------------------------------------------
# 2) HSIC utilities
# -----------------------------------------------------------------------------
rbf_kernel <- function(x, sigma) {
  d <- as.matrix(dist(as.numeric(x)))
  exp(-(d^2) / (2 * sigma^2))
}

# "legacy" HSIC (not used by default; here for completeness)
hsic_legacy <- function(x, y, sigma = NULL) {
  n <- length(x)
  if (is.null(sigma)) {
    dxy <- as.matrix(dist(cbind(as.numeric(x), as.numeric(y))))
    sigma <- median(dxy[lower.tri(dxy)])
    if (!is.finite(sigma) || sigma <= 0) sigma <- 1
  }
  Kx <- rbf_kernel(x, sigma); Ky <- rbf_kernel(y, sigma)
  H  <- diag(n) - matrix(1/n, n, n)
  sum((H %*% Kx %*% H) * Ky) / n^2
}

# "paper" HSIC: separate medians + (n-1)^2 denominator
hsic_paper <- function(x, y) {
  n  <- length(x)
  dx <- as.matrix(dist(as.numeric(x)))
  dy <- as.matrix(dist(as.numeric(y)))
  sx <- median(dx[lower.tri(dx)]); if (!is.finite(sx) || sx <= 0) sx <- 1
  sy <- median(dy[lower.tri(dy)]); if (!is.finite(sy) || sy <= 0) sy <- 1
  Kx <- rbf_kernel(x, sx); Ky <- rbf_kernel(y, sy)
  H  <- diag(n) - matrix(1/n, n, n)
  sum((H %*% Kx %*% H) * Ky) / (n - 1)^2
}

make_hsic <- function(mode = c("legacy","paper")) {
  mode <- match.arg(mode)
  if (mode == "legacy") {
    function(x, y, sigma = NULL) hsic_legacy(x, y, sigma)
  } else {
    function(x, y, sigma = NULL) hsic_paper(x, y)
  }
}

# -----------------------------------------------------------------------------
# 3) Literature-based reference graph (for eval only; NOT used in tuning)
# -----------------------------------------------------------------------------
get_wine_ground_truth <- function(colnames_vec) {
  vars <- colnames_vec
  A <- matrix(0, length(vars), length(vars), dimnames = list(vars, vars))
  need <- c("fixed acidity","volatile acidity","citric acid","residual sugar",
            "chlorides","free sulfur dioxide","total sulfur dioxide","density",
            "pH","sulphates","alcohol","quality")
  missing <- setdiff(need, vars)
  if (length(missing) > 0) stop("Missing expected variables: ", paste(missing, collapse = ", "))

  A["fixed acidity",        "pH"]                   <- 1
  A["volatile acidity",     "pH"]                   <- 1
  A["citric acid",          "pH"]                   <- 1
  A["citric acid",          "fixed acidity"]        <- 1
  A["residual sugar",       "density"]              <- 1
  A["alcohol",              "density"]              <- 1  # inverse in raw scale
  A["free sulfur dioxide",  "total sulfur dioxide"] <- 1
  A["alcohol",              "quality"]              <- 1
  A["volatile acidity",     "quality"]              <- 1
  A["sulphates",            "quality"]              <- 1
  A
}
truth <- get_wine_ground_truth(colnames(wine))

# -----------------------------------------------------------------------------
# 4) LOO discovery (+ DAG enforcement via minimum-score edge removal)
# -----------------------------------------------------------------------------
detect_cycle <- function(adjacency) {
  n <- nrow(adjacency); color <- rep(0, n)
  cyc <- NULL
  dfs <- function(v) {
    color[v] <<- 1
    for (u in which(adjacency[v, ] == 1)) {
      if (color[u] == 1) { cyc <<- c(u, v); return(TRUE) }
      if (color[u] == 0) { if (dfs(u)) return(TRUE) }
    }
    color[v] <<- 2; FALSE
  }
  for (v in 1:n) if (color[v] == 0 && dfs(v)) break
  cyc
}

enforce_acyclicity <- function(adjacency, weights, show_progress = TRUE) {
  n <- nrow(adjacency); maxit <- n * n; it <- 0
  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = maxit, style = 3)
  }
  repeat {
    cyc <- detect_cycle(adjacency)
    if (is.null(cyc) || it >= maxit) break
    v <- cyc[2]
    inc <- which(adjacency[, v] == 1)
    if (!length(inc)) break
    wts <- weights[v, inc]
    rm  <- inc[which.min(wts)]
    adjacency[rm, v] <- 0
    it <- it + 1
    if (show_progress) setTxtProgressBar(pb, it)
  }
  if (show_progress) close(pb)
  adjacency
}

discover_loo <- function(data,
                         hsic_mode = c("legacy","paper"),
                         sigma = NULL,
                         threshold_percentile = 97.5,
                         enforce_dag = TRUE,
                         show_progress = TRUE) {
  hsic_fun <- make_hsic(hsic_mode)
  p <- ncol(data); n <- nrow(data); vars <- colnames(data)
  score <- matrix(0, p, p, dimnames = list(vars, vars))

  total_pairs <- p * (p - 1)
  if (show_progress) {
    cat(sprintf("Scoring %d ordered pairs with HSIC...\n", total_pairs))
    pb <- txtProgressBar(min = 0, max = total_pairs, style = 3); k_done <- 0L
  }

  for (i in 1:p) for (j in 1:p) {
    if (i == j) { if (show_progress) { k_done <- k_done + 1L; setTxtProgressBar(pb, k_done) }; next }
    preds <- setdiff(1:p, c(i, j))
    if (length(preds) > 0) {
      k_basis <- min(5, floor(n / 40))  # small, stable bases; REML smoothing
      rhs <- paste0("s(`", vars[preds], "`, k=", k_basis, ")", collapse = " + ")
      fml <- as.formula(paste0("`", vars[i], "` ~ ", rhs))
      mdl <- tryCatch(
        suppressWarnings(gam(fml, data = data, method = "REML", gamma = 1.4)),
        error = function(e) lm(as.formula(paste0("`", vars[i], "` ~ ",
                                                 paste0("`", vars[preds], "`", collapse = " + "))),
                               data = data)
      )
      ri <- residuals(mdl)
    } else {
      ri <- data[[i]] - mean(data[[i]])
    }
    score[i, j] <- hsic_fun(ri, data[[j]], sigma = sigma)  # large => j -> i
    if (show_progress) { k_done <- k_done + 1L; setTxtProgressBar(pb, k_done) }
  }
  if (show_progress) close(pb)

  # Fixed high-percentile threshold on pooled HSIC scores (no reference used)
  pos <- score[score > 0]
  thr <- as.numeric(quantile(pos, threshold_percentile / 100, na.rm = TRUE))
  cat(sprintf("Threshold (%.1f%% percentile) = %.6f\n", threshold_percentile, thr))

  adj <- matrix(0, p, p, dimnames = list(vars, vars))
  for (i in 1:p) for (j in 1:p) if (i != j && score[i, j] > thr) adj[j, i] <- 1
  cat(sprintf("Edges selected before DAG enforcement: %d\n", sum(adj)))

  # Enforce DAG by removing the minimum-score edge on each detected cycle
  if (enforce_dag) {
    cat("Enforcing acyclicity (removing minimum-score edges on cycles)...\n")
    adj <- enforce_acyclicity(adj, score, show_progress = TRUE)
    cat(sprintf("Edges after DAG enforcement: %d\n", sum(adj)))
  }

  list(adjacency = adj, scores = score, threshold = thr)
}

# -----------------------------------------------------------------------------
# 5) LiNGAM baseline (fastICA heuristic)
# -----------------------------------------------------------------------------
run_lingam <- function(data) {
  t0 <- Sys.time()
  n <- ncol(data)
  dag <- matrix(0, n, n, dimnames = list(colnames(data), colnames(data)))
  cat("Running LiNGAM baseline...\n")
  try({
    ica <- fastICA(as.matrix(data), n.comp = n)
    W <- ica$W
    B <- solve(W); diag(B) <- 0
    B[abs(B) < 0.05] <- 0
    dag <- (B != 0) * 1
  }, silent = TRUE)
  list(dag = dag, time = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

# -----------------------------------------------------------------------------
# 6) Directed metrics (includes Misoriented, SHD, MSE)
# -----------------------------------------------------------------------------
metrics_directed <- function(est, tru) {
  stopifnot(all(dim(est) == dim(tru)))
  n <- nrow(tru)
  tp <- sum(est == 1 & tru == 1)
  fp <- sum(est == 1 & tru == 0)
  fn <- sum(est == 0 & tru == 1)
  tn <- sum(est == 0 & tru == 0) - n  # exclude diagonals
  precision <- ifelse(tp + fp > 0, tp/(tp+fp), 0)
  recall    <- ifelse(tp + fn > 0, tp/(tp+fn), 0)
  f1        <- ifelse(precision + recall > 0, 2*precision*recall/(precision+recall), 0)
  accuracy  <- (tp + tn) / (n * (n - 1))
  shd       <- sum(abs(est - tru))
  mse       <- mean((est - tru)^2)
  # Directed misorientation: edge j->i in truth but i->j (and not j->i) in estimate
  misoriented <- sum(tru == t(est) & tru != est) / 2
  c(Precision = precision, Recall = recall, F1_Score = f1,
    Graph_Accuracy = accuracy, Misoriented = misoriented, SHD = shd, MSE = mse)
}

# -----------------------------------------------------------------------------
# 7) Smooth components (robust to spaces in names; univariate plots)
# -----------------------------------------------------------------------------
make_syn <- function(df) { setNames(df, make.names(colnames(df), unique = TRUE)) }

fit_uni_gam <- function(df_syn, resp_syn, pred_syn, k = 8) {
  df2 <- df_syn[, c(resp_syn, pred_syn)]
  df2 <- df2[is.finite(df2[[1]]) & is.finite(df2[[2]]), , drop = FALSE]
  gam(reformulate(sprintf("s(%s, k=%d)", pred_syn, k), response = resp_syn),
      data = df2, method = "REML", select = TRUE)
}

plot_uni <- function(df_syn, resp_syn, pred_syn, xlab, ylab) {
  k <- max(6, min(12, floor(nrow(df_syn)/150)))
  m <- fit_uni_gam(df_syn, resp_syn, pred_syn, k = k)
  rng <- range(df_syn[[pred_syn]], na.rm = TRUE)
  newd <- setNames(data.frame(seq(rng[1], rng[2], length.out = 200)), pred_syn)
  pr   <- predict(m, newdata = newd, type = "link", se.fit = TRUE)
  edf  <- sum(m$edf, na.rm = TRUE)
  dd   <- data.frame(x = newd[[1]], y = pr$fit, se = pr$se.fit)
  ggplot(dd, aes(x, y)) +
    geom_ribbon(aes(ymin = y - 1.96*se, ymax = y + 1.96*se), alpha = 0.15) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.6) +
    labs(x = xlab, y = ylab, subtitle = sprintf("EDF = %.1f", edf)) +
    theme_minimal(base_size = 11) +
    theme(plot.subtitle = element_text(size = 9))
}

save_smooths <- function(orig_df, show_progress = TRUE) {
  if (show_progress) {
    cat("Creating smooth plots (3 panels)...\n")
    pb <- txtProgressBar(min = 0, max = 3, style = 3); k <- 0L
  }
  # map original to syntactic
  syn_df  <- make_syn(orig_df)
  map_syn <- setNames(colnames(syn_df), colnames(orig_df))

  get_name <- function(patts) {
    nms <- names(map_syn)
    for (p in patts) {
      h <- grep(p, nms, ignore.case = TRUE, perl = TRUE)
      if (length(h)) return(map_syn[[ nms[h[1]] ]])
    }
    NA_character_
  }

  rq <- get_name(c("^quality$","\\bqual"))
  rd <- get_name(c("^density$","\\bdens"))
  rp <- get_name(c("^pH$","^ph$","p\\.?h\\b"))
  pv <- get_name(c("volatile\\.?\\s*acidity","\\bvolatile\\b"))
  pa <- get_name(c("^alcohol$"))
  pc <- get_name(c("citric\\.?\\s*acid","^citric$"))

  need <- c(rq, rd, rp, pv, pa, pc)
  if (any(!nzchar(need))) stop("Could not resolve needed columns for smooths.")

  p1 <- plot_uni(syn_df, rq, pv, "Volatile acidity (std.)", "Quality (partial)")
  if (show_progress) { k <- k + 1L; setTxtProgressBar(pb, k) }
  p2 <- plot_uni(syn_df, rd, pa, "Alcohol (std.)",          "Density (partial)")
  if (show_progress) { k <- k + 1L; setTxtProgressBar(pb, k) }
  p3 <- plot_uni(syn_df, rp, pc, "Citric acid (std.)",      "pH (partial)")
  if (show_progress) { k <- k + 1L; setTxtProgressBar(pb, k); close(pb) }

  combo <- grid.arrange(p1, p2, p3, ncol = 3)
  ggsave("figures/wine_smooths.pdf", combo, width = 9.5, height = 3.3, dpi = 300)
  ggsave("figures/wine_smooths.png", combo, width = 9.5, height = 3.3, dpi = 300)
  message("Saved: figures/wine_smooths.pdf and figures/wine_smooths.png")
}

# -----------------------------------------------------------------------------
# 8) Run: discover (fixed percentile) → baseline → metrics → save → smooth plots
# -----------------------------------------------------------------------------
HSIC_MODE <- "paper"          # matches manuscript choice
THRESHOLD_PERCENTILE <- 97.5  # fixed high-percentile on pooled HSIC scores
ENFORCE_DAG <- TRUE           # break cycles by removing minimum-score edge

cat("Running LOO discovery (fixed high-percentile; no tuning)...\n")
t0 <- Sys.time()
disc <- discover_loo(
  wine,
  hsic_mode = HSIC_MODE,
  sigma = NULL,                               # median heuristic inside HSIC
  threshold_percentile = THRESHOLD_PERCENTILE,
  enforce_dag = ENFORCE_DAG,
  show_progress = TRUE
)
time_prop <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

cat(sprintf("Fixed threshold (percentile): %.1f%% | Numeric cutoff: %.6f\n",
            THRESHOLD_PERCENTILE, disc$threshold))
cat(sprintf("Edges selected (final): %d\n", sum(disc$adjacency)))

ling <- run_lingam(wine)

m_prop <- metrics_directed(disc$adjacency, truth)
m_ling <- metrics_directed(ling$dag, truth)

# ---- build final table (INCLUDES Misoriented_Edges) --------------------------
results <- data.frame(
  Method = c("Proposed LOO", "LiNGAM"),
  Precision_Directed = c(m_prop["Precision"], m_ling["Precision"]),
  Recall_Directed    = c(m_prop["Recall"],    m_ling["Recall"]),
  F1_Score_Directed  = c(m_prop["F1_Score"],  m_ling["F1_Score"]),
  Graph_Accuracy     = c(m_prop["Graph_Accuracy"], m_ling["Graph_Accuracy"]),
  Misoriented_Edges  = c(m_prop["Misoriented"],    m_ling["Misoriented"]),
  SHD                = c(m_prop["SHD"],            m_ling["SHD"]),
  MSE                = c(m_prop["MSE"],            m_ling["MSE"]),
  Runtime_s          = c(time_prop, ling$time),
  check.names = FALSE
)

cat("\n=== Directed Metrics (vs literature-based reference; not used for tuning) ===\n")
print(format(results, digits = 3), row.names = FALSE)

# Save results + adjacencies + threshold details
write.csv(results, "results/wine_metrics.csv", row.names = FALSE)
write.csv(disc$adjacency, "results/adjacency_proposed.csv")
write.csv(ling$dag,       "results/adjacency_lingam.csv")
saveRDS(list(scores = disc$scores, threshold = disc$threshold, percentile = THRESHOLD_PERCENTILE),
        "results/proposed_scores_threshold.rds")

# Smooth plots (the panel you want)
cat("\nFitting post-hoc univariate GAMs for smooth plots...\n")
save_smooths(wine, show_progress = TRUE)

cat("\nAll done.\nOutputs written to:\n  results/wine_metrics.csv\n  results/adjacency_proposed.csv\n  results/adjacency_lingam.csv\n  results/proposed_scores_threshold.rds\n  figures/wine_smooths.pdf (and .png)\n")

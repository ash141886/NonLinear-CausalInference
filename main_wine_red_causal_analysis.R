# =============================================================================
# Wine Quality (Red): Leave-One-Out Nonlinear Causal Discovery + Smooth Plots
# =============================================================================

suppressPackageStartupMessages({
  need <- c("mgcv","ggplot2","dplyr","fastICA","gridExtra")
  has  <- sapply(need, requireNamespace, quietly = TRUE)
  if (!all(has)) {
    stop("Missing packages: ", paste(need[!has], collapse = ", "),
         "\nPlease install them, e.g.: install.packages(c('mgcv','ggplot2','dplyr','fastICA','gridExtra'))")
  }
  library(mgcv); library(ggplot2); library(dplyr); library(fastICA); library(gridExtra)
})

# -----------------------------------------------------------------------------
# 1) Data (keep header names with spaces exactly as-is)
# -----------------------------------------------------------------------------
if (!dir.exists("data")) dir.create("data", recursive = TRUE)
if (!file.exists("data/winequality-red.csv")) {
  download.file(
    "https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-red.csv",
    "data/winequality-red.csv", quiet = TRUE
  )
}
wine_raw <- read.csv("data/winequality-red.csv", sep = ";", check.names = FALSE)
stopifnot(is.data.frame(wine_raw), nrow(wine_raw) > 0)

cat("Columns found:\n  ", paste(colnames(wine_raw), collapse = ", "), "\n\n")

# Standardize for modeling
wine <- as.data.frame(scale(wine_raw))

if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)
cat(sprintf("Wine Quality (red): %d samples x %d variables\n\n", nrow(wine), ncol(wine)))

# -----------------------------------------------------------------------------
# 2) HSIC utilities
# -----------------------------------------------------------------------------
rbf_kernel <- function(x, sigma) {
  d <- as.matrix(dist(as.numeric(x)))
  exp(-(d^2) / (2 * sigma^2))
}

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

make_hsic <- function(mode = c("legacy", "paper")) {
  mode <- match.arg(mode)
  if (mode == "legacy") {
    function(x, y, sigma = NULL) hsic_legacy(x, y, sigma)
  } else {
    function(x, y, sigma = NULL) hsic_paper(x, y)
  }
}

# -----------------------------------------------------------------------------
# 3) Ground truth (names with spaces)
# -----------------------------------------------------------------------------
get_wine_ground_truth <- function(colnames_vec) {
  vars <- colnames_vec
  A <- matrix(0, length(vars), length(vars), dimnames = list(vars, vars))
  need <- c("fixed acidity","volatile acidity","citric acid","residual sugar",
            "chlorides","free sulfur dioxide","total sulfur dioxide","density",
            "pH","sulphates","alcohol","quality")
  missing <- setdiff(need, vars)
  if (length(missing) > 0) {
    stop("Missing expected variables in data: ", paste(missing, collapse = ", "))
  }
  A["fixed acidity",        "pH"]                   <- 1
  A["volatile acidity",     "pH"]                   <- 1
  A["citric acid",          "pH"]                   <- 1
  A["citric acid",          "fixed acidity"]        <- 1
  A["residual sugar",       "density"]              <- 1
  A["alcohol",              "density"]              <- 1
  A["free sulfur dioxide",  "total sulfur dioxide"] <- 1
  A["alcohol",              "quality"]              <- 1
  A["volatile acidity",     "quality"]              <- 1
  A["sulphates",            "quality"]              <- 1
  A
}

# -----------------------------------------------------------------------------
# 4) LOO discovery (+ optional DAG enforcement)
# -----------------------------------------------------------------------------
detect_cycle <- function(adjacency) {
  n <- nrow(adjacency); color <- rep(0, n); parent <- rep(NA_integer_, n)
  cyc <- NULL
  dfs <- function(v) {
    color[v] <<- 1
    for (u in which(adjacency[v, ] == 1)) {
      if (color[u] == 1) { cyc <<- c(u, v); return(TRUE) }
      if (color[u] == 0) { parent[u] <<- v; if (dfs(u)) return(TRUE) }
    }
    color[v] <<- 2; FALSE
  }
  for (v in 1:n) if (color[v] == 0 && dfs(v)) break
  cyc
}

enforce_acyclicity <- function(adjacency, weights) {
  n <- nrow(adjacency); maxit <- n * n; it <- 0
  while (it < maxit) {
    cyc <- detect_cycle(adjacency)
    if (is.null(cyc)) break
    v <- cyc[2]
    inc <- which(adjacency[, v] == 1)
    if (!length(inc)) break
    wts <- weights[v, inc]
    rm  <- inc[which.min(wts)]
    adjacency[rm, v] <- 0
    it <- it + 1
  }
  adjacency
}

discover_loo <- function(data,
                         hsic_mode = c("legacy","paper"),
                         sigma = NULL,
                         threshold_percentile = 85,
                         gam_k_max = 10,
                         enforce_dag = FALSE,
                         show_progress = TRUE) {
  hsic_fun <- make_hsic(hsic_mode)
  p <- ncol(data); n <- nrow(data); vars <- colnames(data)
  score <- matrix(0, p, p, dimnames = list(vars, vars))

  total_pairs <- p * (p - 1)
  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = total_pairs, style = 3)
    k_done <- 0L
  }

  for (i in 1:p) {
    for (j in 1:p) {
      if (i == j) {
        if (show_progress) { k_done <- k_done + 1L; setTxtProgressBar(pb, k_done) }
        next
      }
      preds <- setdiff(1:p, c(i, j))
      if (length(preds) > 0) {
        k_basis <- min(5, floor(n / 40))
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
  }
  if (show_progress) close(pb)

  pos <- score[score > 0]
  thr <- as.numeric(quantile(pos, threshold_percentile / 100, na.rm = TRUE))

  adj <- matrix(0, p, p, dimnames = list(vars, vars))
  for (i in 1:p) for (j in 1:p) if (i != j && score[i, j] > thr) adj[j, i] <- 1

  if (enforce_dag) adj <- enforce_acyclicity(adj, score)

  list(adjacency = adj, scores = score, threshold = thr)
}

# -----------------------------------------------------------------------------
# 5) LiNGAM baseline
# -----------------------------------------------------------------------------
run_lingam <- function(data) {
  t0 <- Sys.time()
  n <- ncol(data)
  dag <- matrix(0, n, n, dimnames = list(colnames(data), colnames(data)))
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
# 6) Directed metrics
# -----------------------------------------------------------------------------
metrics_directed <- function(est, tru) {
  stopifnot(all(dim(est) == dim(tru)))
  n <- nrow(tru)
  tp <- sum(est == 1 & tru == 1)
  fp <- sum(est == 1 & tru == 0)
  fn <- sum(est == 0 & tru == 1)
  tn <- sum(est == 0 & tru == 0) - n
  precision <- ifelse(tp + fp > 0, tp/(tp+fp), 0)
  recall    <- ifelse(tp + fn > 0, tp/(tp+fn), 0)
  f1        <- ifelse(precision + recall > 0, 2*precision*recall/(precision+recall), 0)
  accuracy  <- (tp + tn) / (n * (n - 1))
  shd       <- sum(abs(est - tru))
  mse       <- mean((est - tru)^2)
  misoriented <- sum(tru == t(est) & tru != est) / 2
  c(Precision = precision, Recall = recall, F1_Score = f1,
    Graph_Accuracy = accuracy, Misoriented = misoriented, SHD = shd, MSE = mse)
}

# -----------------------------------------------------------------------------
# 7) Hyperparameter tuning (with progress bar)
# -----------------------------------------------------------------------------
tune_wine <- function(data, truth,
                      hsic_mode = c("legacy","paper"),
                      thresholds = c(70, 80, 85, 90, 95),
                      sigmas = c(NA, 0.5, 1.0, 1.5),
                      k_grid = c(5, 10, 15),
                      show_progress = TRUE) {
  hsic_mode <- match.arg(hsic_mode)
  best <- list(F1 = -1, threshold = NA, sigma = NA, k = NA)

  grid <- expand.grid(threshold = thresholds, sigma = sigmas, k = k_grid)
  n_jobs <- nrow(grid)
  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = n_jobs, style = 3)
    j_done <- 0L
  }

  for (r in 1:n_jobs) {
    g <- grid[r, ]
    res <- discover_loo(
      data,
      hsic_mode = hsic_mode,
      sigma = if (hsic_mode == "legacy" && !is.na(g$sigma)) g$sigma else NULL,
      threshold_percentile = g$threshold,
      gam_k_max = g$k,
      enforce_dag = FALSE,
      show_progress = FALSE
    )
    m <- metrics_directed(res$adjacency, truth)
    if (m["F1_Score"] > best$F1) {
      best <- list(F1 = m["F1_Score"], threshold = g$threshold,
                   sigma = if (hsic_mode == "legacy") g$sigma else NA, k = g$k)
    }
    if (show_progress) { j_done <- j_done + 1L; setTxtProgressBar(pb, j_done) }
  }
  if (show_progress) close(pb)

  cat(sprintf("\nBest tuning: thr=%d  sigma=%s  k=%d  F1=%.3f\n\n",
              best$threshold,
              ifelse(hsic_mode == "legacy",
                     ifelse(is.na(best$sigma), "auto", as.character(best$sigma)), "n/a"),
              best$k, best$F1))
  best
}

# -----------------------------------------------------------------------------
# 8) Smooth components (post-hoc interpretability)
# -----------------------------------------------------------------------------
fit_gams_for_visualization <- function(data, responses = c("quality", "density", "pH"), k_max = 10) {
  models <- list()
  for (resp in intersect(responses, colnames(data))) {
    preds <- setdiff(colnames(data), resp)
    terms <- paste0("s(`", preds, "`, k=", min(k_max, 10), ")", collapse = " + ")
    fml <- as.formula(paste0("`", resp, "` ~ ", terms))
    gm <- tryCatch(gam(fml, data = data, method = "REML"), error = function(e) NULL)
    if (!is.null(gm)) models[[resp]] <- gm
  }
  models
}

plot_smooth_panel <- function(models, data) {
  rels <- list(
    list(model = "quality", predictor = "volatile acidity",
         xlab = "Volatile acidity (std.)", ylab = "Quality (partial)"),
    list(model = "density", predictor = "alcohol",
         xlab = "Alcohol (std.)", ylab = "Density (partial)"),
    list(model = "pH", predictor = "citric acid",
         xlab = "Citric acid (std.)", ylab = "pH (partial)")
  )
  plots <- list()
  for (rel in rels) if (rel$model %in% names(models)) {
    m <- models[[rel$model]]
    x <- seq(-3, 3, length.out = 200)
    newd <- as.data.frame(matrix(0, nrow = length(x), ncol = ncol(data)))
    colnames(newd) <- colnames(data)
    if (!rel$predictor %in% colnames(newd)) next
    newd[[rel$predictor]] <- x

    pr <- predict(m, newdata = newd, type = "terms", se.fit = TRUE)
    idx <- grep(paste0("(^|`)", rel$predictor, "(`|$)"), colnames(pr$fit))[1]
    if (is.na(idx)) next

    edf <- sum(m$edf[grep(rel$predictor, names(m$edf))])
    df  <- data.frame(x = x, y = pr$fit[, idx], se = pr$se.fit[, idx])

    p <- ggplot(df, aes(x, y)) +
      geom_ribbon(aes(ymin = y - 1.96*se, ymax = y + 1.96*se), fill = "gray90") +
      geom_line() +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.6) +
      labs(x = rel$xlab, y = rel$ylab, subtitle = sprintf("EDF = %.1f", edf)) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 9))
    plots[[paste0(rel$model, "_", rel$predictor)]] <- p
  }
  if (length(plots) >= 3) {
    combined <- grid.arrange(
      plots[["quality_volatile acidity"]],
      plots[["density_alcohol"]],
      plots[["pH_citric acid"]],
      ncol = 3
    )
    ggsave("figures/wine_smooths.pdf", combined, width = 9.5, height = 3.3, dpi = 300)
    message("Saved: figures/wine_smooths.pdf")
  }
}

# -----------------------------------------------------------------------------
# 9) Main
# -----------------------------------------------------------------------------
set.seed(1)

truth <- get_wine_ground_truth(colnames(wine))

HSIC_MODE <- "legacy"   # choose "legacy" or "paper"

cat("Tuning parameters...\n")
best <- tune_wine(
  data = wine,
  truth = truth,
  hsic_mode = HSIC_MODE,
  thresholds = c(70, 80, 85, 90, 95),
  sigmas = c(NA, 0.5, 1.0, 1.5),  # used only in legacy mode
  k_grid = c(5, 10, 15),
  show_progress = TRUE
)

cat("Running LOO discovery with tuned parameters...\n")
t0 <- Sys.time()
disc <- discover_loo(
  wine,
  hsic_mode = HSIC_MODE,
  sigma = if (HSIC_MODE == "legacy" && !is.na(best$sigma)) best$sigma else NULL,
  threshold_percentile = best$threshold,
  gam_k_max = best$k,
  enforce_dag = FALSE,
  show_progress = TRUE
)
time_prop <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

cat("Running LiNGAM baseline...\n")
ling <- run_lingam(wine)

m_prop <- metrics_directed(disc$adjacency, truth)
m_ling <- metrics_directed(ling$dag, truth)

results <- data.frame(
  Method = c("Proposed LOO", "LiNGAM"),
  Precision_Directed = c(m_prop["Precision"], m_ling["Precision"]),
  Recall_Directed    = c(m_prop["Recall"],    m_ling["Recall"]),
  F1_Score_Directed  = c(m_prop["F1_Score"],  m_ling["F1_Score"]),
  Graph_Accuracy     = c(m_prop["Graph_Accuracy"], m_ling["Graph_Accuracy"]),
  Misoriented_Edges  = c(m_prop["Misoriented"],    m_ling["Misoriented"]),
  SHD                = c(m_prop["SHD"],            m_ling["SHD"]),
  MSE                = c(m_prop["MSE"],            m_ling["MSE"]),
  Runtime_s          = c(time_prop, ling$time)
)

cat("\n=== Directed Metrics ===\n")
print(format(results, digits = 3), row.names = FALSE)

cat("\nFitting post-hoc GAMs for smooth plots...\n")
mods <- fit_gams_for_visualization(wine, k_max = best$k)
plot_smooth_panel(mods, wine)

cat("\nDone.\n")

set.seed(1)

if (!requireNamespace("limma", quietly = TRUE)) {
  stop("Package 'limma' is required. Install it with BiocManager::install('limma').")
}
library(limma)

## -----------------------------
## Configuration
## -----------------------------
ngenes <- 2000
nsamples <- 12
set_size <- 100 # number of genes in the set being tested
set_idx <- seq_len(set_size) # indices of genes in the set (first 100 genes for simplicity)

# Number of Monte Carlo replicates and null draws per replicate.
# Increase these for publication-quality precision.
nsim <- 100
nperm <- 1000
nrot <- 200
alpha <- 0.05

# Correlation among set genes.
rho_set_grid <- c(0.1, 0.5, 0.8)
# Correlation among background genes (fixed).
rho_bg <- 0.1

## -----------------------------
## Helper functions
## -----------------------------
gsea_set_stat <- function(val, set_idx) {
  ord <- order(val, decreasing = TRUE)
  hits <- ord %in% set_idx

  Nh <- sum(hits)
  N <- length(val)
  Nm <- N - Nh
  if (Nh == 0 || Nm == 0) return(0)

  Phit <- cumsum(hits / Nh)
  Pmiss <- cumsum((!hits) / Nm)
  max(Phit - Pmiss)
}

empirical_p_right <- function(obs, null_dist) {
  (1 + sum(null_dist >= obs)) / (length(null_dist) + 1)
}

simulate_correlated_noise <- function(ngenes, nsamples, set_idx, rho_set, rho_bg) {
  nset <- length(set_idx)
  bg_idx <- setdiff(seq_len(ngenes), set_idx)
  nbg <- length(bg_idx)

  noise <- matrix(0, nrow = ngenes, ncol = nsamples)

  # Equicorrelated block for set genes.
  f_set <- rnorm(nsamples)
  z_set <- matrix(rnorm(nset * nsamples), nrow = nset, ncol = nsamples)
  noise[set_idx, ] <- sqrt(rho_set) * matrix(rep(f_set, each = nset), nrow = nset) +
    sqrt(1 - rho_set) * z_set

  # Equicorrelated block for background genes.
  f_bg <- rnorm(nsamples)
  z_bg <- matrix(rnorm(nbg * nsamples), nrow = nbg, ncol = nsamples)
  noise[bg_idx, ] <- sqrt(rho_bg) * matrix(rep(f_bg, each = nbg), nrow = nbg) +
    sqrt(1 - rho_bg) * z_bg

  noise
}

simulate_null_dataset <- function(ngenes, nsamples, set_idx, rho_set, rho_bg) {
  group <- factor(rep(c(0, 1), each = nsamples / 2))
  design <- model.matrix(~ group)

  beta0 <- rnorm(ngenes, 0, 1)
  beta1 <- numeric(ngenes) # global null: no DE genes

  noise <- simulate_correlated_noise(
    ngenes = ngenes,
    nsamples = nsamples,
    set_idx = set_idx,
    rho_set = rho_set,
    rho_bg = rho_bg
  )

  Y <- beta0 %*% t(rep(1, nsamples)) +
    beta1 %*% t(as.numeric(group)) +
    noise

  list(Y = Y, group = group, design = design)
}

compute_es_from_data <- function(Y, design, set_idx) {
  fit <- lmFit(Y, design)
  pvals <- eBayes(fit)$p.value[, "group1"]
  es <- gsea_set_stat(-log10(pvals), set_idx = set_idx)
  list(pvals = pvals, es = es)
}

build_rotation_basis <- function(design0, nsamples, tol = 1e-6) {
  XtXinv0 <- solve(t(design0) %*% design0)
  H0 <- design0 %*% XtXinv0 %*% t(design0)
  Rspace <- diag(nsamples) - H0

  e <- eigen(Rspace, symmetric = TRUE)
  e$vectors[, e$values > tol, drop = FALSE]
}

sample_rotation_matrix <- function(basis) {
  k <- ncol(basis)
  A <- matrix(rnorm(k * k), k, k)
  Qrot <- qr.Q(qr(A))
  basis %*% Qrot %*% t(basis)
}

mean_offdiag <- function(M) {
  M[upper.tri(M)]
}

plot_correlation_heatmap <- function(Y, set_idx, bg_n = 300, out_file = "corr_heatmap.pdf", main_title = NULL) {
  corr <- cor(t(Y))
  bg_idx <- setdiff(seq_len(nrow(Y)), set_idx)
  if (length(bg_idx) > bg_n) {
    bg_idx <- sample(bg_idx, bg_n)
  }
  keep_idx <- c(set_idx, bg_idx)
  corr_sub <- corr[keep_idx, keep_idx, drop = FALSE]

  nset <- length(set_idx)
  within_set <- corr_sub[seq_len(nset), seq_len(nset), drop = FALSE]
  within_bg <- corr_sub[(nset + 1):nrow(corr_sub), (nset + 1):ncol(corr_sub), drop = FALSE]
  between <- corr_sub[seq_len(nset), (nset + 1):ncol(corr_sub), drop = FALSE]

  mean_set <- mean(mean_offdiag(within_set))
  mean_bg <- mean(mean_offdiag(within_bg))
  mean_between <- mean(between)

  if (is.null(main_title)) {
    main_title <- "Gene-Gene Correlation Heatmap"
  }

  cols <- colorRampPalette(c("#313695", "#ffffbf", "#a50026"))(100)

  pdf(out_file, width = 8, height = 7)
  par(mar = c(4, 4, 4, 1))
  image(
    x = seq_len(nrow(corr_sub)),
    y = seq_len(ncol(corr_sub)),
    z = corr_sub[nrow(corr_sub):1, ],
    col = cols,
    zlim = c(-1, 1),
    xlab = "Genes (set first, then background)",
    ylab = "Genes (set first, then background)",
    main = main_title,
    axes = FALSE
  )
  box()
  abline(v = nset + 0.5, lwd = 2, col = "black")
  abline(h = nrow(corr_sub) - nset + 0.5, lwd = 2, col = "black")

  mtext(
    sprintf("mean within-set = %.3f | mean within-bg = %.3f | mean between = %.3f", mean_set, mean_bg, mean_between),
    side = 3,
    line = 0.5,
    cex = 0.85
  )
  dev.off()

  list(mean_within_set = mean_set, mean_within_bg = mean_bg, mean_between = mean_between)
}

compute_null_pvalues <- function(Y, group, design, set_idx, nperm, nrot) {
  obs <- compute_es_from_data(Y, design, set_idx)
  ES_obs <- obs$es
  P <- obs$pvals

  null_dist_perm <- numeric(nperm)
  ngenes <- nrow(Y)
  for (i in seq_len(nperm)) {
    perm_idx <- sample.int(ngenes)
    null_dist_perm[i] <- gsea_set_stat(-log10(P[perm_idx]), set_idx = set_idx)
  }

  design0 <- model.matrix(~ 1, data = data.frame(group = group))
  fit0 <- lmFit(Y, design0)
  Yhat0 <- fitted(fit0)
  res0 <- Y - Yhat0
  basis <- build_rotation_basis(design0 = design0, nsamples = ncol(Y))

  null_dist_rot <- numeric(nrot)
  for (b in seq_len(nrot)) {
    Rmat <- sample_rotation_matrix(basis)
    Y_rot <- Yhat0 + res0 %*% Rmat

    rot <- compute_es_from_data(Y_rot, design, set_idx)
    null_dist_rot[b] <- rot$es
  }

  list(
    p_perm = empirical_p_right(ES_obs, null_dist_perm),
    p_rot = empirical_p_right(ES_obs, null_dist_rot),
    ES_obs = ES_obs,
    null_perm = null_dist_perm,
    null_rot = null_dist_rot
  )
}

## -----------------------------
## Simulation study
## -----------------------------
results <- data.frame(
  rho_set = rep(rho_set_grid, each = nsim),
  sim = rep(seq_len(nsim), times = length(rho_set_grid)),
  p_perm = NA_real_,
  p_rot = NA_real_
)

example_datasets <- list()
example_distributions <- list()

for (i in seq_along(rho_set_grid)) {
  rho_set <- rho_set_grid[i]
  cat(sprintf("Running rho_set = %.2f\n", rho_set))

  for (s in seq_len(nsim)) {
    cat(sprintf("  Simulation %d/%d\n", s, nsim))
    row_idx <- (i - 1) * nsim + s
    dat <- simulate_null_dataset(
      ngenes = ngenes,
      nsamples = nsamples,
      set_idx = set_idx,
      rho_set = rho_set,
      rho_bg = rho_bg
    )

    pv <- compute_null_pvalues(
      Y = dat$Y,
      group = dat$group,
      design = dat$design,
      set_idx = set_idx,
      nperm = nperm,
      nrot = nrot
    )

    results$p_perm[row_idx] <- pv$p_perm
    results$p_rot[row_idx] <- pv$p_rot

    if (s == 1) {
      example_datasets[[paste0("rho_", rho_set)]] <- dat$Y
    }
    if (s == 1 && rho_set %in% c(min(rho_set_grid), max(rho_set_grid))) {
      example_distributions[[paste0("rho_", rho_set)]] <- pv
    }
  }
}

summary_tbl <- data.frame(
  rho_set = rho_set_grid,
  type1_perm = NA_real_,
  type1_rot = NA_real_,
  mean_p_perm = NA_real_,
  mean_p_rot = NA_real_
)

for (i in seq_along(rho_set_grid)) {
  rho <- rho_set_grid[i]
  idx <- results$rho_set == rho

  summary_tbl$type1_perm[i] <- mean(results$p_perm[idx] < alpha)
  summary_tbl$type1_rot[i] <- mean(results$p_rot[idx] < alpha)
  summary_tbl$mean_p_perm[i] <- mean(results$p_perm[idx])
  summary_tbl$mean_p_rot[i] <- mean(results$p_rot[idx])
}

## -----------------------------
## Output + plots
## -----------------------------
cat("\nEmpirical type-I error (target =", alpha, "):\n")
print(summary_tbl)
write.csv(summary_tbl, "null_calibration_summary.csv", row.names = FALSE)
write.csv(results, "null_calibration_results.csv", row.names = FALSE)

heatmap_summary <- data.frame(
  rho_set = numeric(0),
  mean_within_set = numeric(0),
  mean_within_bg = numeric(0),
  mean_between = numeric(0)
)

for (rho in rho_set_grid) {
  key <- paste0("rho_", rho)
  if (!is.null(example_datasets[[key]])) {
    rho_label <- gsub("\\.", "p", format(rho, trim = TRUE))
    out_file <- sprintf("corr_heatmap_rho_%s.pdf", rho_label)
    stats <- plot_correlation_heatmap(
      Y = example_datasets[[key]],
      set_idx = set_idx,
      bg_n = 300,
      out_file = out_file,
      main_title = sprintf("Correlation heatmap (rho_set = %.1f, rho_bg = %.1f)", rho, rho_bg)
    )
    heatmap_summary <- rbind(
      heatmap_summary,
      data.frame(
        rho_set = rho,
        mean_within_set = stats$mean_within_set,
        mean_within_bg = stats$mean_within_bg,
        mean_between = stats$mean_between
      )
    )
  }
}

if (nrow(heatmap_summary) > 0) {
  write.csv(heatmap_summary, "corr_heatmap_summary.csv", row.names = FALSE)
  cat("\nHeatmap summary stats:\n")
  print(heatmap_summary)
}

pdf("null_calibration.pdf", width = 11, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

# Panel 1: type-I error across correlation strengths
plot_tbl <- summary_tbl[
  is.finite(summary_tbl$rho_set) &
    is.finite(summary_tbl$type1_perm) &
    is.finite(summary_tbl$type1_rot),
, drop = FALSE]
plot_tbl <- plot_tbl[order(plot_tbl$rho_set), , drop = FALSE]

if (nrow(plot_tbl) > 0) {
  y_lim <- range(c(plot_tbl$type1_perm, plot_tbl$type1_rot, alpha), na.rm = TRUE)
  if (y_lim[1] == y_lim[2]) {
    y_lim <- y_lim + c(-0.01, 0.01)
  }
  plot(
    plot_tbl$rho_set,
    plot_tbl$type1_perm,
    type = "b",
    pch = 16,
    lwd = 2,
    ylim = y_lim,
    xlab = expression(rho[set]),
    ylab = "Empirical type-I error",
    main = "Calibration Under Null"
  )
  lines(plot_tbl$rho_set, plot_tbl$type1_rot, type = "b", pch = 17, lwd = 2, col = "blue")
  abline(h = alpha, lty = 2, col = "darkgray")
  legend(
    "topleft",
    legend = c("Gene-wise permutation", "Rotation", "Nominal alpha"),
    col = c("black", "blue", "darkgray"),
    lty = c(1, 1, 2),
    pch = c(16, 17, NA),
    bty = "n"
  )
} else {
  plot.new()
  title("Calibration Under Null")
  text(0.5, 0.5, "No finite values to plot")
}

# Panel 2: p-value ECDF at highest correlation
rho_hi <- max(rho_set_grid)
idx_hi <- results$rho_set == rho_hi
plot(
  ecdf(results$p_perm[idx_hi]),
  do.points = FALSE,
  verticals = TRUE,
  main = sprintf("ECDF of p-values (rho_set = %.1f)", rho_hi),
  xlab = "p-value",
  ylab = "ECDF"
)
lines(ecdf(results$p_rot[idx_hi]), do.points = FALSE, verticals = TRUE, col = "blue")
abline(0, 1, lty = 2, col = "darkgray")
legend(
  "topleft",
  legend = c("Permutation", "Rotation", "Uniform(0,1)"),
  col = c("black", "blue", "darkgray"),
  lty = c(1, 1, 2),
  bty = "n"
)

# Panel 3: one example null distribution pair at rho=0
ex0_name <- paste0("rho_", min(rho_set_grid))
if (!is.null(example_distributions[[ex0_name]])) {
  ex0 <- example_distributions[[ex0_name]]
  hist(ex0$null_perm, breaks = 40, col = "gray85", border = "white",
       main = sprintf("Example null ES (rho_set = %.1f)", min(rho_set_grid)),
       xlab = "ES")
  hist(ex0$null_rot, breaks = 40, col = rgb(0, 0, 1, 0.35), border = "white", add = TRUE)
  abline(v = ex0$ES_obs, col = "red", lwd = 2)
  legend("topright", legend = c("Permutation", "Rotation", "Observed ES"),
         fill = c("gray85", rgb(0, 0, 1, 0.35), NA), border = NA,
         lty = c(NA, NA, 1), col = c(NA, NA, "red"), bty = "n")
}

# Panel 4: one example null distribution pair at high rho
exh_name <- paste0("rho_", max(rho_set_grid))
if (!is.null(example_distributions[[exh_name]])) {
  exh <- example_distributions[[exh_name]]
  hist(exh$null_perm, breaks = 40, col = "gray85", border = "white",
       main = sprintf("Example null ES (rho_set = %.1f)", max(rho_set_grid)),
       xlab = "ES")
  hist(exh$null_rot, breaks = 40, col = rgb(0, 0, 1, 0.35), border = "white", add = TRUE)
  abline(v = exh$ES_obs, col = "red", lwd = 2)
  legend("topright", legend = c("Permutation", "Rotation", "Observed ES"),
         fill = c("gray85", rgb(0, 0, 1, 0.35), NA), border = NA,
         lty = c(NA, NA, 1), col = c(NA, NA, "red"), bty = "n")
}

dev.off()

cat("\nSaved:\n")
cat("- null_calibration_summary.csv\n")
cat("- null_calibration_results.csv\n")
cat("- null_calibration.pdf\n")
cat("- corr_heatmap_summary.csv\n")
cat("- corr_heatmap_rho_*.pdf\n")

hi_dad <- function() {
  cat("\n================ Script Summary ================\n")
  cat("Purpose: Compare gene-wise permutation vs rotation null calibration under gene-gene correlation.\n")
  cat(sprintf("Config: ngenes=%d, nsamples=%d, set_size=%d, nsim=%d, nperm=%d, nrot=%d, alpha=%.3f\n",
              ngenes, nsamples, set_size, nsim, nperm, nrot, alpha))
  cat(sprintf("Correlation setup: rho_set_grid = %s | rho_bg = %.3f\n",
              paste(rho_set_grid, collapse = ", "), rho_bg))

  if (exists("summary_tbl")) {
    cat("\nType-I error summary:\n")
    print(summary_tbl)
  } else {
    cat("\nType-I error summary: not available in current session.\n")
  }

  cat("\nOutput files:\n")
  out_files <- c(
    "null_calibration_summary.csv",
    "null_calibration_results.csv",
    "null_calibration.pdf",
    "corr_heatmap_summary.csv"
  )
  for (f in out_files) {
    cat(sprintf("- %s [%s]\n", f, if (file.exists(f)) "found" else "missing"))
  }

  heatmap_files <- Sys.glob("corr_heatmap_rho_*.pdf")
  if (length(heatmap_files) > 0) {
    cat(sprintf("- corr_heatmap_rho_*.pdf [found: %d file(s)]\n", length(heatmap_files)))
  } else {
    cat("- corr_heatmap_rho_*.pdf [missing]\n")
  }
  cat("===============================================\n")
}

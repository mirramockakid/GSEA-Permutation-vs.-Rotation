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
nperm <- 500
nrot <- 500
alpha <- 0.05

# Parallelization controls:
#   "none"         = fully sequential
#   "nsim"         = parallelize simulations within each rho_set value
#   "rho_set_grid" = parallelize across rho_set values
parallel_over <- "nsim"
detected_cores <- parallel::detectCores(logical = FALSE)
if (!is.finite(detected_cores)) {
  detected_cores <- 1L
}
ncores <- max(1L, detected_cores - 1L)

# Correlation among set genes.
rho_set_grid <- c(0.0, 0.2, 0.4, 0.6, 0.8)
# Correlation among background genes (fixed).
rho_bg <- 0.1

# Cache simulation outputs so plotting/output steps can be rerun without
# repeating the Monte Carlo simulation.
use_simulation_cache <- TRUE
simulation_cache_file <- "simulation_cache.rds"

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
  fit <- limma::lmFit(Y, design)
  pvals <- limma::eBayes(fit)$p.value[, "group1"]
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

plot_rotation_correlation_preservation <- function(
    Y,
    set_idx,
    bg_n = 300,
    out_file = "rotation_correlation_preservation.pdf",
    main_title = NULL
) {
  if (is.null(main_title)) {
    main_title <- "Rotation preserves gene-gene correlation structure"
  }

  nsamples_local <- ncol(Y)
  design0 <- matrix(1, nrow = nsamples_local, ncol = 1)
  fit0 <- limma::lmFit(Y, design0)
  Yhat0 <- fitted(fit0)
  res0 <- Y - Yhat0
  basis <- build_rotation_basis(design0 = design0, nsamples = nsamples_local)
  Rmat <- sample_rotation_matrix(basis)
  Y_rot <- Yhat0 + res0 %*% Rmat

  bg_idx <- setdiff(seq_len(nrow(Y)), set_idx)
  if (length(bg_idx) > bg_n) {
    bg_idx <- sample(bg_idx, bg_n)
  }
  keep_idx <- c(set_idx, bg_idx)

  corr_orig <- cor(t(Y[keep_idx, , drop = FALSE]))
  corr_rot <- cor(t(Y_rot[keep_idx, , drop = FALSE]))
  upper_idx <- upper.tri(corr_orig)
  corr_orig_vec <- corr_orig[upper_idx]
  corr_rot_vec <- corr_rot[upper_idx]

  pairwise_corr_of_corr <- suppressWarnings(cor(corr_orig_vec, corr_rot_vec, use = "complete.obs"))
  mean_abs_diff <- mean(abs(corr_orig_vec - corr_rot_vec), na.rm = TRUE)

  cols <- colorRampPalette(c("#313695", "#ffffbf", "#a50026"))(100)

  pdf(out_file, width = 12, height = 4.5)
  par(mfrow = c(1, 3), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))

  image(
    x = seq_len(nrow(corr_orig)),
    y = seq_len(ncol(corr_orig)),
    z = corr_orig[nrow(corr_orig):1, ],
    col = cols,
    zlim = c(-1, 1),
    xlab = "Genes",
    ylab = "Genes",
    main = "Original",
    axes = FALSE
  )
  box()

  image(
    x = seq_len(nrow(corr_rot)),
    y = seq_len(ncol(corr_rot)),
    z = corr_rot[nrow(corr_rot):1, ],
    col = cols,
    zlim = c(-1, 1),
    xlab = "Genes",
    ylab = "Genes",
    main = "Rotated",
    axes = FALSE
  )
  box()

  plot(
    corr_orig_vec,
    corr_rot_vec,
    pch = 16,
    cex = 0.35,
    col = rgb(0, 0, 0, 0.2),
    xlab = "Original pairwise correlation",
    ylab = "Rotated pairwise correlation",
    main = "Correlation-by-correlation"
  )
  abline(0, 1, col = "red", lwd = 2)
  mtext(
    sprintf("r = %.3f | mean abs diff = %.3f", pairwise_corr_of_corr, mean_abs_diff),
    side = 3,
    line = 0.2,
    cex = 0.85
  )

  mtext(main_title, side = 3, outer = TRUE, line = 0.3, cex = 0.95)
  dev.off()

  list(
    pairwise_corr_of_corr = pairwise_corr_of_corr,
    mean_abs_diff = mean_abs_diff,
    n_genes_plotted = length(keep_idx),
    n_pairs = length(corr_orig_vec)
  )
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
  fit0 <- limma::lmFit(Y, design0)
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

parallel_lapply <- function(X, FUN, cores) {
  cores <- max(1L, as.integer(cores))
  if (cores <= 1L || length(X) <= 1L) {
    return(lapply(X, FUN))
  }

  if (.Platform$OS.type == "windows") {
    message("Windows detected: falling back to sequential lapply for this script-level parallel helper.")
    lapply(X, FUN)
  } else {
    parallel::mclapply(X, FUN, mc.cores = cores, mc.preschedule = TRUE)
  }
}

run_one_simulation <- function(rho_set, sim_id, seed, capture_dataset = FALSE, capture_distribution = FALSE) {
  set.seed(seed)

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

  list(
    p_perm = pv$p_perm,
    p_rot = pv$p_rot,
    example_Y = if (isTRUE(capture_dataset)) dat$Y else NULL,
    example_pv = if (isTRUE(capture_distribution)) pv else NULL
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
parallel_over <- match.arg(parallel_over, c("none", "nsim", "rho_set_grid"))
ncores <- max(1L, min(as.integer(ncores), parallel::detectCores()))
if (!is.finite(ncores)) {
  ncores <- 1L
}

cache_config <- list(
  ngenes = ngenes,
  nsamples = nsamples,
  set_size = set_size,
  set_idx = set_idx,
  nsim = nsim,
  nperm = nperm,
  nrot = nrot,
  rho_set_grid = rho_set_grid,
  rho_bg = rho_bg,
  alpha = alpha
)

cache_loaded <- FALSE
if (isTRUE(use_simulation_cache) && file.exists(simulation_cache_file)) {
  cache_obj <- readRDS(simulation_cache_file)
  if (
    is.list(cache_obj) &&
      !is.null(cache_obj$config) &&
      identical(cache_obj$config, cache_config) &&
      !is.null(cache_obj$results) &&
      !is.null(cache_obj$example_datasets) &&
      !is.null(cache_obj$example_distributions)
  ) {
    results <- cache_obj$results
    example_datasets <- cache_obj$example_datasets
    example_distributions <- cache_obj$example_distributions
    cache_loaded <- TRUE
    cat(sprintf("Loaded cached simulation outputs from %s\n", simulation_cache_file))
  } else {
    cat(sprintf("Ignoring %s (cache missing fields or config mismatch)\n", simulation_cache_file))
  }
}

if (!cache_loaded) {
  cat(sprintf("Parallel mode: %s | cores: %d\n", parallel_over, ncores))

  # Deterministic per-simulation seeds so results do not depend on scheduling.
  seed_grid <- matrix(
    sample.int(.Machine$integer.max, length(rho_set_grid) * nsim, replace = FALSE),
    nrow = length(rho_set_grid),
    ncol = nsim
  )

  if (parallel_over == "rho_set_grid") {
    rho_jobs <- lapply(seq_along(rho_set_grid), function(i) {
      list(i = i, rho_set = rho_set_grid[i], seeds = seed_grid[i, ])
    })

    rho_results <- parallel_lapply(rho_jobs, function(job) {
      rho_set <- job$rho_set
      sim_out <- vector("list", nsim)

      for (s in seq_len(nsim)) {
        sim_out[[s]] <- run_one_simulation(
          rho_set = rho_set,
          sim_id = s,
          seed = job$seeds[s],
          capture_dataset = (s == 1),
          capture_distribution = (s == 1 && rho_set %in% c(min(rho_set_grid), max(rho_set_grid)))
        )
      }

      list(
        i = job$i,
        rho_set = rho_set,
        sim_out = sim_out
      )
    }, cores = ncores)

    for (rho_res in rho_results) {
      i <- rho_res$i
      rho_set <- rho_res$rho_set
      cat(sprintf("Completed rho_set = %.2f\n", rho_set))

      for (s in seq_len(nsim)) {
        row_idx <- (i - 1) * nsim + s
        sim_res <- rho_res$sim_out[[s]]
        results$p_perm[row_idx] <- sim_res$p_perm
        results$p_rot[row_idx] <- sim_res$p_rot

        if (s == 1) {
          example_datasets[[paste0("rho_", rho_set)]] <- sim_res$example_Y
          if (!is.null(sim_res$example_pv)) {
            example_distributions[[paste0("rho_", rho_set)]] <- sim_res$example_pv
          }
        }
      }
    }
  } else {
    for (i in seq_along(rho_set_grid)) {
      rho_set <- rho_set_grid[i]
      cat(sprintf("Running rho_set = %.2f\n", rho_set))

      sim_ids <- as.list(seq_len(nsim))
      if (parallel_over == "nsim") {
        sim_out <- parallel_lapply(sim_ids, function(s) {
          run_one_simulation(
            rho_set = rho_set,
            sim_id = s,
            seed = seed_grid[i, s],
            capture_dataset = (s == 1),
            capture_distribution = (s == 1 && rho_set %in% c(min(rho_set_grid), max(rho_set_grid)))
          )
        }, cores = ncores)
      } else {
        sim_out <- lapply(sim_ids, function(s) {
          cat(sprintf("  Simulation %d/%d\n", s, nsim))
          run_one_simulation(
            rho_set = rho_set,
            sim_id = s,
            seed = seed_grid[i, s],
            capture_dataset = (s == 1),
            capture_distribution = (s == 1 && rho_set %in% c(min(rho_set_grid), max(rho_set_grid)))
          )
        })
      }

      for (s in seq_len(nsim)) {
        row_idx <- (i - 1) * nsim + s
        sim_res <- sim_out[[s]]
        results$p_perm[row_idx] <- sim_res$p_perm
        results$p_rot[row_idx] <- sim_res$p_rot

        if (s == 1) {
          example_datasets[[paste0("rho_", rho_set)]] <- sim_res$example_Y
          if (!is.null(sim_res$example_pv)) {
            example_distributions[[paste0("rho_", rho_set)]] <- sim_res$example_pv
          }
        }
      }
    }
  }

  if (isTRUE(use_simulation_cache)) {
    saveRDS(
      list(
        config = cache_config,
        results = results,
        example_datasets = example_datasets,
        example_distributions = example_distributions
      ),
      simulation_cache_file
    )
    cat(sprintf("Saved simulation outputs to %s\n", simulation_cache_file))
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
  p_perm_i <- results$p_perm[idx]
  p_rot_i <- results$p_rot[idx]

  n_missing_perm <- sum(!is.finite(p_perm_i))
  n_missing_rot <- sum(!is.finite(p_rot_i))
  if (n_missing_perm > 0 || n_missing_rot > 0) {
    warning(
      sprintf(
        "rho_set=%.2f has non-finite p-values (perm: %d/%d, rot: %d/%d). Summary uses finite values only.",
        rho, n_missing_perm, length(p_perm_i), n_missing_rot, length(p_rot_i)
      )
    )
  }

  summary_tbl$type1_perm[i] <- mean(p_perm_i < alpha, na.rm = TRUE)
  summary_tbl$type1_rot[i] <- mean(p_rot_i < alpha, na.rm = TRUE)
  summary_tbl$mean_p_perm[i] <- mean(p_perm_i, na.rm = TRUE)
  summary_tbl$mean_p_rot[i] <- mean(p_rot_i, na.rm = TRUE)
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

rotation_preservation_summary <- data.frame(
  rho_set = numeric(0),
  pairwise_corr_of_corr = numeric(0),
  mean_abs_diff = numeric(0),
  n_genes_plotted = integer(0),
  n_pairs = integer(0)
)

for (rho in unique(c(min(rho_set_grid), max(rho_set_grid)))) {
  key <- paste0("rho_", rho)
  if (!is.null(example_datasets[[key]])) {
    rho_label <- gsub("\\.", "p", format(rho, trim = TRUE))
    out_file <- sprintf("rotation_corr_preservation_rho_%s.pdf", rho_label)
    stats <- plot_rotation_correlation_preservation(
      Y = example_datasets[[key]],
      set_idx = set_idx,
      bg_n = 300,
      out_file = out_file,
      main_title = sprintf("Rotation correlation preservation (rho_set = %.1f, rho_bg = %.1f)", rho, rho_bg)
    )
    rotation_preservation_summary <- rbind(
      rotation_preservation_summary,
      data.frame(
        rho_set = rho,
        pairwise_corr_of_corr = stats$pairwise_corr_of_corr,
        mean_abs_diff = stats$mean_abs_diff,
        n_genes_plotted = stats$n_genes_plotted,
        n_pairs = stats$n_pairs
      )
    )
  }
}

if (nrow(rotation_preservation_summary) > 0) {
  write.csv(rotation_preservation_summary, "rotation_corr_preservation_summary.csv", row.names = FALSE)
  cat("\nRotation correlation-preservation summary:\n")
  print(rotation_preservation_summary)
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

if (nrow(plot_tbl) == 0) {
  warning("Panel 'Calibration Under Null' has no finite rows in summary_tbl. Check summary_tbl and results for NA/non-finite p-values.")
}

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
cat("- rotation_corr_preservation_summary.csv\n")
cat("- rotation_corr_preservation_rho_*.pdf\n")
if (isTRUE(use_simulation_cache)) {
  cat(sprintf("- %s\n", simulation_cache_file))
}

set.seed(1)

if (!requireNamespace("limma", quietly = TRUE)) {
  stop("Package 'limma' is required. Install it with BiocManager::install('limma').")
}
library(limma)

## -----------------------------
## Configuration
## -----------------------------
ngenes <- 10000
nsamples <- 8
nperm <- 10000
nrot <- 100
set_idx <- 1:100

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

simulate_expression <- function(ngenes, nsamples, set_idx) {
  group <- factor(rep(c(0, 1), each = nsamples / 2))
  design <- model.matrix(~ group)

  beta0 <- rnorm(ngenes, 0, 1)
  beta1 <- numeric(ngenes)
  beta1[set_idx] <- 1.5

  noise <- matrix(rnorm(ngenes * nsamples), nrow = ngenes, ncol = nsamples)

  Y <- beta0 %*% t(rep(1, nsamples)) +
    beta1 %*% t(as.numeric(group)) +
    noise

  list(Y = Y, group = group, design = design)
}

compute_es_from_data <- function(Y, design, set_idx) {
  fit <- lmFit(Y, design)
  pvals <- eBayes(fit)$p.value[, "group1"]
  es <- gsea_set_stat(-log10(pvals), set_idx = set_idx)
  list(fit = fit, pvals = pvals, es = es)
}

build_rotation_basis <- function(design0, nsamples, tol = 1e-6) {
  XtXinv0 <- solve(t(design0) %*% design0)
  H0 <- design0 %*% XtXinv0 %*% t(design0)
  Rspace <- diag(nsamples) - H0

  e <- eigen(Rspace, symmetric = TRUE)
  basis <- e$vectors[, e$values > tol, drop = FALSE]
  basis
}

sample_rotation_matrix <- function(basis) {
  k <- ncol(basis)
  A <- matrix(rnorm(k * k), k, k)
  Qrot <- qr.Q(qr(A))
  basis %*% Qrot %*% t(basis)
}

empirical_p_right <- function(obs, null_dist) {
  (1 + sum(null_dist >= obs)) / (length(null_dist) + 1)
}

## -----------------------------
## Simulate data + observed ES
## -----------------------------
sim <- simulate_expression(
  ngenes = ngenes, nsamples = nsamples, set_idx = set_idx
  )
Y <- sim$Y
group <- sim$group
design <- sim$design

obs <- compute_es_from_data(Y, design, set_idx)
fit <- obs$fit
P <- obs$pvals
ES_obs <- obs$es

## -----------------------------
## Permutation null (gene-label permutations)
## -----------------------------
null_dist_perm <- numeric(nperm)
for (i in seq_len(nperm)) {
  perm_idx <- sample.int(ngenes)
  null_dist_perm[i] <- gsea_set_stat(-log10(P[perm_idx]), set_idx = set_idx)
}

## -----------------------------
## Rotation null
## -----------------------------
# Use null-model residuals so rotation samples the null distribution for group effect.
design0 <- model.matrix(~ 1, data = data.frame(group = group))
fit0 <- lmFit(Y, design0)
Yhat0 <- fitted(fit0)
res0 <- Y - Yhat0

basis <- build_rotation_basis(design0 = design0, nsamples = nsamples)

null_dist_rot <- numeric(nrot)
max_diff_vec <- numeric(nrot)
corr_orig <- cor(t(Y))

for (b in seq_len(nrot)) {
  Rmat <- sample_rotation_matrix(basis)

  Y_rot <- Yhat0 + res0 %*% Rmat

  rot <- compute_es_from_data(Y_rot, design, set_idx)
  null_dist_rot[b] <- rot$es

  corr_rot <- cor(t(Y_rot))
  max_diff_vec[b] <- max(abs(corr_rot - corr_orig))
}

# Compute empirical p-values
p_perm <- empirical_p_right(ES_obs, null_dist_perm)
p_rot <- empirical_p_right(ES_obs, null_dist_rot)

## -----------------------------
## Diagnostics + plots
## -----------------------------
cat("Dimensions of Y:", dim(Y)[1], "x", dim(Y)[2], "\n")
cat("Observed ES:", ES_obs, "\n")
cat("Empirical p-value (permutation):", p_perm, "\n")
cat("Empirical p-value (rotation):", p_rot, "\n")
cat("Max |delta(corr)| summary:\n")
print(summary(max_diff_vec))

pdf("es_dist.pdf", width = 8, height = 4)
par(mfrow = c(1, 2))
hist(null_dist_rot, breaks = 50, main = "Null distribution (rotation)", xlab = "GSEA statistic")
abline(v = ES_obs, col = "red", lwd = 2)
legend("topright", legend = sprintf("p = %.4g", p_rot), bty = "n")
hist(null_dist_perm, breaks = 50, main = "Null distribution (permutation)", xlab = "GSEA statistic")
abline(v = ES_obs, col = "red", lwd = 2)
legend("topright", legend = sprintf("p = %.4g", p_perm), bty = "n")
dev.off()

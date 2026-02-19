set.seed(1)
if (!requireNamespace("limma", quietly = TRUE)) {
  stop("Package 'limma' is required. Install it with BiocManager::install('limma').")
}
library(limma)

gsea_set_stat <- function(val, set_idx) {
  # Rank genes by statistic (descending)
  ord <- order(val, decreasing = TRUE)
  hits <- ord %in% set_idx
  
  Nh <- sum(hits)
  N  <- length(val)
  Nm <- N - Nh
  if (Nh == 0 || Nm == 0) return(0)
  
  # Running sums for hits and misses
  Phit  <- cumsum(hits / Nh)
  Pmiss <- cumsum((!hits) / Nm)
  
  # Enrichment score (max deviation)
  ES <- max(Phit - Pmiss)
  ES
}

## 1. Simulate a gene expression matrix
ngenes   <- 10000      # number of genes
nsamples <- 8        # number of samples

group    <- factor(rep(c(0, 1), each = nsamples/2))
design   <- model.matrix(~ group)   # intercept + group

# Gene-level coefficients
beta0 <- rnorm(ngenes, 0, 1)       # baseline
beta1 <- numeric(ngenes)
beta1[1:20] <- 1.5                 # first 20 genes are DE

# Simulate noise
E <- matrix(rnorm(ngenes * nsamples), nrow = ngenes, ncol = nsamples)

# Linear mean structure: Y = beta0 + beta1 * group + error
Y <- beta0 %*% t(rep(1, nsamples)) +
     beta1 %*% t(as.numeric(group)) +
     E

dim(Y)      # 100 x 8


## 2. Fit linear model (per-gene)
fit <- lmFit(Y, design)
## get p-values for group effect
P <- eBayes(fit)$p.value[, "group1"]
ES_obs <- gsea_set_stat(-log10(P), set_idx = 1:20)

# create null distribution by permuting gene labels
null_dist_perm <- numeric(10000)
for (i in 1:10000) {
  perm_idx <- sample(ngenes)
  null_dist_perm[i] <- gsea_set_stat(-log10(P[perm_idx]), set_idx = 1:20)
}

# plot histogram of null distribution and observed statistic
hist(null_dist_perm, breaks = 50, main = "Null distribution (permutation)", 
xlab = "GSEA statistic")
abline(v = ES_obs, col = "red", lwd = 2)




# For rotation-based null generation, use residuals from the null model.
# Rotating residuals from the full model keeps the tested t-statistic unchanged.
design0 <- model.matrix(~ 1, data = data.frame(group = group))
fit0 <- lmFit(Y, design0)
Yhat0 <- fitted(fit0)
res0 <- Y - Yhat0

## Residual space projector under the null model
X0      <- design0
XtXinv0 <- solve(t(X0) %*% X0)
H0      <- X0 %*% XtXinv0 %*% t(X0)
Rspace  <- diag(nsamples) - H0




## 3. Construct and apply a random rotation in the residual space

# Get an orthonormal basis for the residual space via eigen-decomposition
e <- eigen(Rspace, symmetric = TRUE)
# eigenvalues ~ 1 correspond to residual space (rank = nsamples - rank(X))
basis <- e$vectors[, e$values > 1e-6, drop = FALSE]
k <- ncol(basis)   # dimension of residual space

## 2. Perform 100 rotations and record max |Δcor| ----------------

nrot <- 100
null_dist_rot <- numeric(nrot)
max_diff_vec <- numeric(nrot)
corr_orig <- cor(t(Y))

for (b in 1:nrot) {
  # random k×k orthogonal matrix
  A    <- matrix(rnorm(k * k), k, k)
  Qrot <- qr.Q(qr(A))
  
  # nsamples×nsamples rotation in residual space
  Rmat <- basis %*% Qrot %*% t(basis)
  
  # rotated residuals and data
  res_rot <- res0 %*% Rmat
  Y_rot   <- Yhat0 + res_rot

  fit_rot <- lmFit(Y_rot, design)
  P_rot <- eBayes(fit_rot)$p.value[, "group1"]
  null_dist_rot[b] <- gsea_set_stat(-log10(P_rot), set_idx = 1:20)

  corr_rot <- cor(t(Y_rot))
  max_diff_vec[b] <- max(abs(corr_rot - corr_orig))
}

pdf("rot.pdf", width = 6, height = 4)
hist(null_dist_rot, breaks = 50, main = "Null distribution (rotation)")
hist(null_dist_perm, breaks = 50, main = "Null distribution (permutation)")
dev.off()

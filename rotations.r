# install limma
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("limma")



# simulate data
set.seed(123)

nGenes <- 20
nSamples <- 6

# Random expression values
expr <- matrix(rnorm(nGenes * nSamples, mean=5, sd=1),
               nrow=nGenes, ncol=nSamples)
rownames(expr) <- paste0("Gene", 1:nGenes)
colnames(expr) <- paste0("Sample", 1:nSamples)

# Add a small signal to first 3 genes in treatment
expr[1:3, 4:6] <- expr[1:3, 4:6] + 2

# Create design matrix for 2 groups (control vs treatment)
agroup <- factor(c("Control","Control","Control","Treatment","Treatment","Treatment"))
design <- model.matrix(~group)
design

# gene set
geneSet <- c("Gene1","Gene2","Gene3")


# Access residuals
# ---- 2. Fit linear model ----
fit <- lmFit(expr, design)
residMat <- residuals(fit, expr)   # genes x samples

# Observed gene set statistic (mean difference between groups)
obs_stat <- rowMeans(expr[geneSet, group=="Treatment"]) - rowMeans(expr[geneSet, group=="Control"])
obs_stat <- mean(obs_stat)
obs_stat

# ---- 3. Create null distribution via rotations ----
nrot <- 1000                     # number of rotations
null_stats <- numeric(nrot)      # to store rotated gene set statistics

for (i in 1:nrot) {
  # Generate random orthogonal matrix (samples x samples)
  R <- qr.Q(qr(matrix(rnorm(nSamples^2), nSamples, nSamples)))
  
  # Rotate residuals
  rotated_resid <- residMat %*% R
  
  # Reconstruct pseudo-expression under null
  y_null <- fit$coefficients %*% t(design) + rotated_resid
  
  # Compute gene set statistic for this rotation
  stat <- rowMeans(y_null[geneSet, group=="Treatment"]) -
          rowMeans(y_null[geneSet, group=="Control"])
  null_stats[i] <- mean(stat)
}

# ---- 4. Visualize null distribution ----
hist(null_stats, breaks=30, main="Null distribution via rotations",
     xlab="Gene set statistic", col="lightblue")
abline(v=obs_stat, col="red", lwd=2)
legend("topright", legend="Observed statistic", col="red", lwd=2)

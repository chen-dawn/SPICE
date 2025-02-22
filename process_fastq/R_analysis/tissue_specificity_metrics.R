# Load necessary libraries
library(dplyr)

# Define the scoring metrics as functions
counts <- function(x, threshold) {
  t <- ifelse(x > threshold, 1, 0)
  return((length(x) - sum(t)) / (length(x) - 1))
}

tau <- function(x) {
  x_hat <- x / max(x)
  return(sum(1 - x_hat) / (length(x) - 1))
}

gini <- function(x) {
  x_sorted <- sort(x)
  n <- length(x)
  gini_val <- sum((2 * 1:n - n - 1) * x_sorted) / (n * sum(x_sorted))
  return(gini_val * n / (n - 1))
}

simpson <- function(x) {
  p <- x / sum(x)
  simpson_val <- sum(p^2)
  return((simpson_val - 1 / length(x)) / (1 - 1 / length(x)))
}

shannon_specificity <- function(x) {
  p <- x / sum(x)
  H <- -sum(p * log2(p))
  HS <- log2(length(x)) - H
  return(HS / log2(length(x)))
}

roku_specificity <- function(x) {
  M <- median(x)
  S <- median(abs(x - M))
  u <- (x - M) / (5 * S + 1e-4)
  w <- ifelse(abs(u) <= 1, (1 - u^2)^2, 0)
  t <- sum(x * w) / sum(w)
  x_prime <- abs(x - t)
  p <- x_prime / sum(x_prime)
  H <- -sum(p * log2(p))
  ROKU <- log2(length(x)) - H
  return(ROKU / log2(length(x)))
}

spm_dpm <- function(x) {
  SPM <- x^2 / (x * sqrt(sum(x^2)))
  mean_SPM <- mean(SPM)
  sigma <- sqrt(sum((SPM - mean_SPM)^2) / (length(SPM) - 1))
  return(sigma * sqrt(length(SPM)))
}

jss_dpm <- function(x) {
  p <- x / sum(x)
  q <- rep(1 / length(x), length(x))
  m <- (p + q) / 2
  H <- -sum(m * log2(m))
  JSS <- 1 - sqrt((sum(p * log2(p)) + sum(q * log2(q))) / 2 - H)
  mean_JSS <- mean(JSS)
  sigma <- sqrt(sum((JSS - mean_JSS)^2) / (length(JSS) - 1))
  return(sigma * sqrt(length(JSS)))
}

tsi <- function(x) {
  return(x / sum(x))
}

z_score <- function(x) {
  mean_x <- mean(x)
  sigma <- sd(x)
  z_scores <- (x - mean_x) / sigma
  return((z_scores + (length(x) - 1) / sqrt(length(x))) / (2 * (length(x) - 1) / sqrt(length(x))))
}

spm <- function(x) {
  return(x^2 / (x * sqrt(sum(x^2))))
}

jss <- function(x) {
  p <- x / sum(x)
  q <- rep(1 / length(x), length(x))
  m <- (p + q) / 2
  H <- -sum(m * log2(m))
  return(1 - sqrt((sum(p * log2(p)) + sum(q * log2(q))) / 2 - H))
}
# Define the vectors
vectors <- list(
  vector1 = c(0, 0, 0, 0, 0.1) + 1,
  vector2 = c(0, 0, 0, 0, 0.1),
  vector3 = c(0.01, 0.02, 0.03, 0.01, 0.01) + 1,
  vector4 = c(0.1, 0.1, 0, 0, 0.9) + 1,
  vector5 = c(0.3, 0.4, 0.5, 0.6, 0.7) + 1,
  vector6 = c(0, 1, 1, 1, 1) + 1,
  vector7 = c(0.5, 1, 1, 1, 1) + 1
)

# Threshold for Counts
threshold <- 0.15

# Initialize a dataframe to store the results
results <- data.frame(
  Metric = c("Counts", "Tau", "Gini", "Simpson", "Shannon Specificity", 
             "ROKU Specificity", "SPM DPM", "JSS DPM", "TSI", 
             "Z-Score", "SPM", "JSS"),
  Vector1 = numeric(12),
  Vector2 = numeric(12),
  Vector3 = numeric(12),
  Vector4 = numeric(12),
  Vector5 = numeric(12),
  Vector6 = numeric(12),
  Vector7 = numeric(12)
)

# Calculate the values for each metric and each vector
for (i in 1:length(vectors)) {
  vec <- vectors[[i]]
  results[1, i + 1] <- counts(vec, threshold)
  results[2, i + 1] <- tau(vec) *2
  results[3, i + 1] <- gini(vec)
  results[4, i + 1] <- simpson(vec)
  results[5, i + 1] <- shannon_specificity(vec)
  results[6, i + 1] <- roku_specificity(vec)
  results[7, i + 1] <- spm_dpm(vec)
  results[8, i + 1] <- jss_dpm(vec)
}

View(results)

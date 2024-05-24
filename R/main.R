# main.R for diffCircadian

# Load the required libraries
require(minpack.lm)
require(parallel)
require(parallelly)
require(lmtest)
require(tidyverse)

# Load the required functions
fitSinCurve <- function(xx, observed, parStart = list(A=3, phase=0, offset=0)) {
  getPred <- function(parS, xx) {
    parS$A * sin((1/12*pi)*(xx+parS$phase)) + parS$offset
  }
  residFun <- function(parS, observed, xx) {
    getPred(parS, xx) - observed
  }
  nls.out <- nls.lm(par=parStart, fn=residFun, observed=observed, xx=xx)
  apar <- nls.out$par

  A0 <- apar$A

  asign <- sign(A0)

  A <- A0 * asign

  phase <- (apar$phase + ifelse(asign==1, 0,12)) %% 24

  offset <- apar$offset

  peak <- (12 * sign(A0) - 6 - phase) %% 24

  if(peak > 18) peak = peak - 24

  SSE <- sum(nls.out$fvec^2)
  SST <- sum((observed-mean(observed))^2)
  R2 <- 1 - SSE/SST
  res <- list(A=A, phase=phase, offset=offset, peak=peak, R2=R2)
  return(res)
}

# Since there is only one observation per sample...
generate_observed_parameters <- function(zeitgeber_times, observed) {
  out <- fitSinCurve(xx=zeitgeber_times, observed=observed)
  observed_parameters <- c(out$A, out$phase, out$offset, out$peak, out$R2)
  return(observed_parameters)
}

ciradianDrawing <- function(tod, expr, apar, labels, specInfo=NULL) {
  geneName <- apar$geneName
  peak <- round(apar$peak)
  if(peak==18) peak <- -6
  amain <- paste("Circadian pattern of", geneName, "expression")
  times <- seq(-6, 18, 0.1)
  pred <- getPred(apar, times)
  curve_fun <- function(x) {
    apar$A * sin((1/12*pi)*(x+apar$phase)) + apar$offset
  }
  plot(tod, expr, xlab="Zeitgeber time (hours)", ylab="Expression level", main=amain, xlim=c(-6, 18), ylim=c(min(expr), max(expr)))
  lines(times, curve_fun(times), col="red")
  box(which="plot", lty="solid", lwd=3)
}

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# Get p-values through lrtest with the null hypothesis that the model is a linear regression
get_lrtest <- function(residuals, zeitgeber_times) {
  lrtest_results <- list()
  lrtest_results <- lrtest(lm(range01(residuals) ~ zeitgeber_times), lm(asin(range01(residuals)) ~ zeitgeber_times))
  lrtest_results_pvalue <- pchisq(lrtest_results$Chisq[2], df=lrtest_results$`#Df`[2], lower.tail=FALSE)
  return(lrtest_results_pvalue)
}

# Load the data
data <- read_csv("Data/data_impute_zeit.csv")
# Impute missing values with the geometric mean
data <- data %>% mutate_at(vars(-Cage), funs(replace(., is.na(.), exp(mean(log(.[!is.na(.)]))))))
pheno_data <- read_csv("Data/pheno.csv")
zeitgeber_times <- data$Cage
snames <- colnames(data)[2:ncol(data)]
A_samples <- pheno_data$Cage[pheno_data$Dose=="A"]
A_samples <- gsub("X", "", A_samples)
B_samples <- pheno_data$Cage[pheno_data$Dose=="B"]
B_samples <- gsub("X", "", B_samples)
C_samples <- pheno_data$Cage[pheno_data$Dose=="C"]
C_samples <- gsub("X", "", C_samples)
data_A <- data[, c(A_samples)]
data_B <- data[, c(B_samples)]
data_C <- data[, c(C_samples)]


# Fit the model
observed_parameters_A <- apply(data_A, 2, generate_observed_parameters, zeitgeber_times=zeitgeber_times)
observed_parameters_A <- as.data.frame(t(observed_parameters_A))
colnames(observed_parameters_A) <- c("A", "phase", "offset", "peak", "R2")

# Fit the model
observed_parameters_B <- apply(data_B, 2, generate_observed_parameters, zeitgeber_times=zeitgeber_times)
observed_parameters_B <- as.data.frame(t(observed_parameters_B))
colnames(observed_parameters_B) <- c("A", "phase", "offset", "peak", "R2")

# Fit the model
observed_parameters_C <- apply(data_C, 2, generate_observed_parameters, zeitgeber_times=zeitgeber_times)
observed_parameters_C <- as.data.frame(t(observed_parameters_C))
colnames(observed_parameters_C) <- c("A", "phase", "offset", "peak", "R2")

# Get p-values
pvalues_A <- apply(data_A, 2, get_lrtest, zeitgeber_times=zeitgeber_times)
pvalues_B <- apply(data_B, 2, get_lrtest, zeitgeber_times=zeitgeber_times)
pvalues_C <- apply(data_C, 2, get_lrtest, zeitgeber_times=zeitgeber_times)

observed_parameters_A <- cbind(observed_parameters_A, pvalues_A)
observed_parameters_B <- cbind(observed_parameters_B, pvalues_B)
observed_parameters_C <- cbind(observed_parameters_C, pvalues_C)

observed_parameters_A <- cbind(Cage=colnames(data_A), observed_parameters_A)
observed_parameters_B <- cbind(Cage=colnames(data_B), observed_parameters_B)
observed_parameters_C <- cbind(Cage=colnames(data_C), observed_parameters_C)

# Save the results
write_csv(observed_parameters_A, "Results/observed_parameters_A.csv")
write_csv(observed_parameters_B, "Results/observed_parameters_B.csv")
write_csv(observed_parameters_C, "Results/observed_parameters_C.csv")

# Get p-values for difference in A, phase, offset, peak and R2
get_p_value_var <- function(data_1, data_2, var) {
  p_value <- t.test(data_1[, var], data_2[, var])$p.value
  return(p_value)
}

for (var in c("A", "phase", "offset", "peak", "R2")) {
  p_value_A_B <- get_p_value_var(observed_parameters_A, observed_parameters_B, var)
  p_value_A_C <- get_p_value_var(observed_parameters_A, observed_parameters_C, var)
  p_value_B_C <- get_p_value_var(observed_parameters_B, observed_parameters_C, var)
  print(paste("P-value for difference in", var, "between A and B:", p_value_A_B))
  print(paste("P-value for difference in", var, "between A and C:", p_value_A_C))
  print(paste("P-value for difference in", var, "between B and C:", p_value_B_C))
}


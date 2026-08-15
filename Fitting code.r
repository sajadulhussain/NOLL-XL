


## Required packages
library(GenSA)
library(numDeriv)
library(goftest)

## ----------------------------------------------------------
## Dataset
## ----------------------------------------------------------

x <- c(
  1.1, 4.1, 1.4, 1.8, 1.3, 1.5, 1.7, 1.2, 1.9, 1.4,
  1.8, 3.0, 1.6, 1.7, 2.2, 2.3, 1.7, 1.6, 2.7, 2.0
)

n <- length(x)


## ----------------------------------------------------------
## NOLL--XL PDF
## ----------------------------------------------------------

dNOLLXL <- function(x, alpha, beta, theta, log = FALSE) {

  if (alpha <= 0 || beta <= 0 || theta <= 0)
    return(rep(ifelse(log, -Inf, 0), length(x)))

  term1 <- 1 + theta * x / (1 + theta)^2

  e_theta <- exp(-theta * x)
  e_beta  <- exp(-beta * theta * x)

  A <- 1 - term1 * e_theta

  ## Numerical protection
  A <- pmax(A, 1e-12)
  term1 <- pmax(term1, 1e-12)
  e_beta <- pmax(e_beta, 1e-300)

  numerator <-
    (theta^2 * (2 + theta + x) * e_beta *
       term1^(beta - 1) *
       A^(alpha - 1) *
       (1 + theta)^(-2)) *
    (alpha + (beta - alpha) *
       (1 - term1 * exp(-theta * x)))

  denominator <-
    (A^alpha + e_beta * term1^beta)^2

  f <- numerator / denominator

  f[!is.finite(f) | f <= 0] <- 1e-16

  if (log)
    return(log(f))

  f
}


## ----------------------------------------------------------
## Negative log-likelihood
## ----------------------------------------------------------

nll_NOLLXL <- function(par, data) {

  alpha <- par[1]
  beta  <- par[2]
  theta <- par[3]

  if (alpha <= 0 || beta <= 0 || theta <= 0)
    return(1e100)

  logf <- log(
    dNOLLXL(
      data,
      alpha = alpha,
      beta = beta,
      theta = theta,
      log = FALSE
    )
  )

  if (any(!is.finite(logf)))
    return(1e100)

  -sum(logf)
}


## ----------------------------------------------------------
## Maximum likelihood estimation
## ----------------------------------------------------------

set.seed(123)

fit <- GenSA(
  par = c(alpha = 1, beta = 1, theta = 1),
  fn = nll_NOLLXL,
  data = x,
  lower = c(alpha = 1e-5,
            beta  = 1e-5,
            theta = 1e-5),
  upper = c(alpha = 10,
            beta  = 10,
            theta = 10),
  control = list(
    max.call = 8000,
    smooth = FALSE
  )
)

## MLEs
alpha_hat <- fit$par[1]
beta_hat  <- fit$par[2]
theta_hat <- fit$par[3]

## Log-likelihood
negloglik <- fit$value
loglik <- -negloglik

cat("MLEs:\n")
cat("alpha =", alpha_hat, "\n")
cat("beta  =", beta_hat, "\n")
cat("theta =", theta_hat, "\n")
cat("-logL =", negloglik, "\n")
cat("LogL  =", loglik, "\n")


## ----------------------------------------------------------
## Standard errors from observed Fisher information
## ----------------------------------------------------------

H <- hessian(
  func = nll_NOLLXL,
  x = fit$par,
  data = x
)

V <- tryCatch(
  solve(H),
  error = function(e) NULL
)

if (!is.null(V) && all(is.finite(V))) {

  SE <- sqrt(abs(diag(V)))

  SE_alpha <- SE[1]
  SE_beta  <- SE[2]
  SE_theta <- SE[3]

} else {

  SE_alpha <- NA
  SE_beta  <- NA
  SE_theta <- NA
}

cat("\nStandard Errors:\n")
cat("SE(alpha) =", SE_alpha, "\n")
cat("SE(beta)  =", SE_beta, "\n")
cat("SE(theta) =", SE_theta, "\n")


## ----------------------------------------------------------
## AIC and BIC
## ----------------------------------------------------------

k <- 3

AIC_NOLLXL <- -2 * loglik + 2 * k
BIC_NOLLXL <- -2 * loglik + k * log(n)

cat("\nInformation Criteria:\n")
cat("AIC =", AIC_NOLLXL, "\n")
cat("BIC =", BIC_NOLLXL, "\n")


## ----------------------------------------------------------
## CDF
##
## For the NOLL--XL construction:
## F(x) = A(x)^alpha /
##        [A(x)^alpha + B(x)^beta]
##
## where
## A(x) = G(x)
## B(x) = 1-G(x)
## ----------------------------------------------------------

pNOLLXL <- function(x, alpha, beta, theta) {

  G <- 1 -
    (1 + theta * x / (1 + theta)^2) *
    exp(-theta * x)

  G <- pmin(pmax(G, 1e-12), 1 - 1e-12)

  F <- G^alpha /
    (G^alpha + (1 - G)^beta)

  F <- pmin(pmax(F, 0), 1)

  F
}


## ----------------------------------------------------------
## Kolmogorov--Smirnov test
## ----------------------------------------------------------

ks <- suppressWarnings(
  ks.test(
    x,
    function(q)
      pNOLLXL(
        q,
        alpha = alpha_hat,
        beta = beta_hat,
        theta = theta_hat
      )
  )
)

KS <- as.numeric(ks$statistic)
KS_pvalue <- ks$p.value

cat("\nKS test:\n")
cat("KS =", KS, "\n")
cat("p-value =", KS_pvalue, "\n")


## ----------------------------------------------------------
## Cramer--von Mises test
## ----------------------------------------------------------

cvm <- suppressWarnings(
  cvm.test(
    x,
    null = function(q)
      pNOLLXL(
        q,
        alpha = alpha_hat,
        beta = beta_hat,
        theta = theta_hat
      )
  )
)

CM <- as.numeric(cvm$statistic)
CM_pvalue <- cvm$p.value

cat("\nCramer--von Mises test:\n")
cat("CM =", CM, "\n")
cat("p-value =", CM_pvalue, "\n")


## ----------------------------------------------------------
## Final results
## ----------------------------------------------------------

results <- data.frame(
  Model = "NOLL-XL",
  Alpha = alpha_hat,
  SE_Alpha = SE_alpha,
  Beta = beta_hat,
  SE_Beta = SE_beta,
  Theta = theta_hat,
  SE_Theta = SE_theta,
  NegLogLik = negloglik,
  LogLik = loglik,
  AIC = AIC_NOLLXL,
  BIC = BIC_NOLLXL,
  KS = KS,
  KS_pvalue = KS_pvalue,
  CM = CM,
  CM_pvalue = CM_pvalue
)

print(results, digits = 8, row.names = FALSE)
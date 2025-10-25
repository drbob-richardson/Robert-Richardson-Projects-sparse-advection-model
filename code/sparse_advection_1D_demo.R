#' Sparse Implied–Advection Utilities (1D demo)
#'
#' Self-contained helpers for the 1D implied-advection + NNGP example.
#' Requires: Matrix, geoR, tmvtnorm
#'

# ---- packages (fail gracefully) ----
req <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. Please install it.", pkg), call. = FALSE)
  }
}
req("Matrix"); req("geoR"); req("tmvtnorm")

# Shorthands
Diagonal <- Matrix::Diagonal
cholSparse <- function(A) Matrix::Cholesky(A, LDL = FALSE, perm = TRUE)

# ---- Nearest-neighbor structure ----
#' Build lower-triangular NNGP neighbor pairs (Vecchia order)
#' @param D distance matrix (k x k)
#' @param m neighbors per site (m << k)
#' @return two-column integer matrix [row i = (site, neighbor)]
make_neighbors <- function(D, m) {
  k <- nrow(D)
  if (m < 1 || m >= k) stop("m must be in [1, k-1]")
  Nbrs <- cbind(2L, 1L)  # site 2 neighbors site 1
  if (k >= 3) {
    for (i in 3:k) {
      # only previous sites are eligible neighbors
      prev <- 1:(i - 1)
      ord  <- order(D[i, prev])[1:min(m, length(prev))]
      Nbrs <- rbind(Nbrs, cbind(i, prev[ord]))
    }
  }
  Nbrs
}

# ---- NNGP precision constructor ----
#' Build NNGP precision C^{-1} = B^T F^{-1} B (sparse)
#' @param D distance matrix (k x k)
#' @param Nstruct neighbor pairs from make_neighbors()
#' @param cov_fun function(dist_matrix) -> covariance matrix for given indices
#' @param k number of locations
#' @return sparse precision matrix (k x k)
make_Cinv <- function(D, Nstruct, cov_fun, k) {
  # Sparse identity B and diagonal F
  B <- Matrix::Diagonal(k)
  Finv_diag <- rep(1, k)
  
  for (i in 2:k) {
    nbrs <- Nstruct[Nstruct[, 1] == i, 2]
    CNs  <- cov_fun(D[nbrs, nbrs, drop = FALSE])
    CsNs <- cov_fun(D[i, nbrs, drop = FALSE, drop = TRUE])
    # Solve CNs * x = CsNs  (stable for small m)
    coef <- as.numeric(solve(CNs, CsNs))
    B[i, nbrs] <- -coef
    Finv_diag[i] <- 1 / (cov_fun(D[i, i, drop = FALSE]) - sum(CsNs * coef))
  }
  Finv <- Diagonal(x = Finv_diag)
  Matrix::t(B) %*% Finv %*% B
}

# ---- Backward-difference advection–diffusion: G^{-1} (1D) ----
#' Construct sparse inverse evolution G^{-1} in 1D (implied advection)
#' @param mu vector of advection (length k)
#' @param sig vector of diffusion (length k, positive)
#' @return sparse k x k matrix Ginv
make_Ginv_1d <- function(mu, sig) {
  k <- length(mu)
  if (length(sig) != k) stop("mu and sig must have same length")
  G <- Matrix::Diagonal(k, x = 1 + sig)
  # boundaries
  G[1, 1]    <- G[1, 1] - sig[1]/2 + mu[1]/2
  G[k, k]    <- G[k, k] - sig[k]/2 - mu[k]/2
  # sub/super diagonals
  if (k >= 2) {
    G[cbind(2:k, 1:(k - 1))] <- -sig[2:k]/2 + mu[2:k]/2
  }
  if (k >= 3) {
    G[cbind(1:(k - 1), 2:k)] <- -sig[2:(k - 1)]/2 - mu[2:(k - 1)]/2
  }
  Matrix::Matrix(G, sparse = TRUE)
}

# ---- Sparse KF + FFBS (with optional taper) ----
#' Sparse extended-KF mean/precision pass + FFBS draw
#' @param Y (k x T) data
#' @param m0 prior mean (k)
#' @param C0inv prior precision (k x k, sparse)
#' @param Ft list length T of (k x k) sparse obs matrices
#' @param Gtinv list length T of (k x k) sparse inverse-evolution matrices
#' @param v obs variance scalar
#' @param delta discount factor in (0,1)
#' @param taper_mask optional logical sparse matrix (pattern for tapering)
#' @return list(Ctinv=..., Rtinv=..., draw= X_{0:T})
sparse_kf_ffbs <- function(Y, m0, C0inv, Ft, Gtinv, v, delta, taper_mask = NULL) {
  k <- nrow(Y); Tlen <- ncol(Y)
  Rtinv <- vector("list", Tlen)
  Ctinv <- vector("list", Tlen)
  mt    <- Matrix::Matrix(0, k, Tlen, sparse = TRUE)
  at    <- mt
  
  # t=1
  Rtinv[[1]] <- delta * (Gtinv[[1]] %*% C0inv %*% Matrix::t(Gtinv[[1]]))
  Ctinv[[1]] <- Rtinv[[1]] + Matrix::t(Ft[[1]]) %*% Ft[[1]] / v
  
  at[, 1] <- solve(Gtinv[[1]], m0, sparse = TRUE)
  # Innovation-precision trick to avoid big solves
  S1 <- Rtinv[[1]] + Matrix::t(Ft[[1]]) %*% Ft[[1]] / v
  rhs <- (Y[, 1] - Ft[[1]] %*% at[, 1]) / v
  mt[, 1] <- at[, 1] + solve(Rtinv[[1]], Matrix::t(Ft[[1]]) %*% (rhs - (Matrix::t(Ft[[1]]) %*% solve(S1, Matrix::t(Ft[[1]]), sparse = TRUE) %*% rhs) / v), sparse = TRUE)
  
  # taper (optional)
  if (!is.null(taper_mask)) {
    L <- Matrix::Cholesky(Ctinv[[1]], LDL = FALSE, perm = TRUE)
    Lt <- Matrix::as(L, "sparseMatrix")
    Lt@x[!as.vector(taper_mask@x)] <- 0
    Ctinv[[1]] <- Lt %*% Matrix::t(Lt)
  }
  
  # t=2..T
  for (t in 2:Tlen) {
    Rtinv[[t]] <- delta * (Gtinv[[t]] %*% Ctinv[[t - 1]] %*% Matrix::t(Gtinv[[t]]))
    Ctinv[[t]] <- Rtinv[[t]] + Matrix::t(Ft[[t]]) %*% Ft[[t]] / v
    at[, t]    <- solve(Gtinv[[t]], mt[, t - 1], sparse = TRUE)
    
    St <- Rtinv[[t]] + Matrix::t(Ft[[t]]) %*% Ft[[t]] / v
    rhs <- (Y[, t] - Ft[[t]] %*% at[, t]) / v
    mt[, t] <- at[, t] + solve(Rtinv[[t]], Matrix::t(Ft[[t]]) %*% (rhs - (Matrix::t(Ft[[t]]) %*% solve(St, Matrix::t(Ft[[t]]), sparse = TRUE) %*% rhs) / v), sparse = TRUE)
    
    if (!is.null(taper_mask)) {
      L <- Matrix::Cholesky(Ctinv[[t]], LDL = FALSE, perm = TRUE)
      Lt <- Matrix::as(L, "sparseMatrix")
      Lt@x[!as.vector(taper_mask@x)] <- 0
      Ctinv[[t]] <- Lt %*% Matrix::t(Lt)
    }
  }
  
  # FFBS draw
  rands <- matrix(stats::rnorm(k * (Tlen + 1)), k, Tlen + 1)
  draw  <- Matrix::Matrix(0, k, Tlen + 1, sparse = TRUE)
  
  # final precision approx for smoother noise
  cholT <- cholSparse(Ctinv[[Tlen]])
  draw[, Tlen + 1] <- mt[, Tlen] + Matrix::solve(cholT, rands[, Tlen + 1])
  
  if (Tlen >= 2) {
    for (t in seq(Tlen, 2, by = -1)) {
      cholPrev <- cholSparse(Ctinv[[t - 1]])
      draw[, t] <- mt[, t - 1] +
        delta * (Gtinv[[t]] %*% (draw[, t + 1] - at[, t])) +
        Matrix::solve(cholPrev, rands[, t])
    }
  }
  # t = 1 uses prior
  chol0 <- cholSparse(C0inv)
  draw[, 1] <- m0 + delta * (Gtinv[[1]] %*% (draw[, 2] - at[, 1])) + Matrix::solve(chol0, rands[, 1])
  
  list(Ctinv = Ctinv, Rtinv = Rtinv, draw = draw)
}

# ---- Sufficient-stat makers for mu and sigma (1D) ----
# NOTE: explicitly pass Y; fix colon slicing; avoid hidden n/T
make_ct_mu <- function(mu, sig, samp, Y) {
  k      <- length(mu)
  Tp1    <- ncol(samp)
  Cmak   <- -Matrix::Diagonal(k, x = 1 + sig)
  if (k >= 2) Cmak[cbind(2:k, 1:(k - 1))] <-  sig[2:k]/2
  if (k >= 3) Cmak[cbind(1:(k - 1), 2:k)] <-  sig[2:(k - 1)]/2
  cts <- samp[, 2:(Tp1 - 1)] + Cmak %*% samp[, 2:Tp1]
  cts[1, ] <- cts[1, ] + sig[1]/2 * Y[1, ]
  cts[k, ] <- cts[k, ] + sig[k]/2 * Y[k, ]
  cts
}

make_Dt_mu <- function(mu, sig, samp, Y) {
  k    <- length(mu)
  Tp1  <- ncol(samp)
  Dmak <- Matrix::Matrix(0, k, k, sparse = TRUE)
  if (k >= 2) Dmak[cbind(2:k, 1:(k - 1))] <-  1/2
  if (k >= 2) Dmak[cbind(1:(k - 1), 2:k)] <- -1/2
  Dts <- Dmak %*% samp[, 2:Tp1]
  Dts[1, ] <- Dts[1, ] + 1/2 * Y[1, ]
  Dts[k, ] <- Dts[k, ] - 1/2 * Y[k, ]
  Dts
}

make_ct_sig <- function(mu, sig, samp, Y) {
  k      <- length(mu)
  Tp1    <- ncol(samp)
  Cmak   <- -Matrix::Diagonal(k)
  if (k >= 2) Cmak[cbind(2:k, 1:(k - 1))] <- -mu[2:k]/2
  if (k >= 3) Cmak[cbind(1:(k - 1), 2:k)] <-  mu[2:(k - 1)]/2
  cts <- samp[, 2:(Tp1 - 1)] + Cmak %*% samp[, 2:Tp1]
  cts[1, ] <- cts[1, ] - mu[1]/2 * Y[1, ]
  cts[k, ] <- cts[k, ] + mu[k]/2 * Y[k, ]
  cts
}

make_Dt_sig <- function(mu, sig, samp, Y) {
  k    <- length(mu)
  Tp1  <- ncol(samp)
  Dmak <- Matrix::Diagonal(k)
  if (k >= 2) Dmak[cbind(2:k, 1:(k - 1))] <- -1/2
  if (k >= 2) Dmak[cbind(1:(k - 1), 2:k)] <- -1/2
  Dts <- Dmak %*% samp[, 2:Tp1]
  Dts[1, ] <- Dts[1, ] - 1/2 * Y[1, ]
  Dts[k, ] <- Dts[k, ] - 1/2 * Y[k, ]
  Dts
}

# ---- Conjugate/posterior samplers for mu, sig (vectorized sitewise) ----
# mu ~ N(Sig0 * b, Sig0), no truncation
samp_mu <- function(Ctinv_list, C0inv, Dts, delta, bt, mu0, Sig0inv) {
  k <- nrow(Dts); Tlen <- ncol(Dts)
  Vpostinv <- Sig0inv
  bpost    <- Sig0inv %*% mu0
  for (t in 1:Tlen) {
    Dmat <- Diagonal(x = as.numeric(Dts[, t]))
    Prec <- if (t == 1) C0inv else Ctinv_list[[t - 1]]
    Vpostinv <- Vpostinv + (Dmat %*% Prec %*% Dmat) * (delta / (1 - delta))
    bpost    <- bpost    + (Dmat %*% Prec %*% Dmat) %*% (bt[, t] * delta / (1 - delta))
  }
  mu <- Matrix::solve(Vpostinv, bpost, sparse = TRUE)
  drop(mu)
}

# sig >= 0 via truncated MVN (fallback to prior mode if NA)
samp_sig <- function(Ctinv_list, C0inv, Dts, delta, bt, mu0, Sig0inv, sig_fallback) {
  k <- nrow(Dts); Tlen <- ncol(Dts)
  Vpostinv <- Sig0inv
  bpost    <- Sig0inv %*% mu0
  for (t in 1:Tlen) {
    Dmat <- Diagonal(x = as.numeric(Dts[, t]))
    Prec <- if (t == 1) C0inv else Ctinv_list[[t - 1]]
    Vpostinv <- Vpostinv + (Dmat %*% Prec %*% Dmat) * (delta / (1 - delta))
    bpost    <- bpost    + (Dmat %*% Prec %*% Dmat) %*% (bt[, t] * delta / (1 - delta))
  }
  mn <- Matrix::solve(Vpostinv, bpost, sparse = TRUE)
  out <- tmvtnorm::rtmvnorm.sparseMatrix(
    n = 1, mean = as.numeric(mn), H = Vpostinv,
    lower = rep(0, k), upper = rep(Inf, k)
  )
  if (any(is.na(out))) out <- sig_fallback
  as.numeric(out)
}

# ---- Example driver (small demo) ----
#' Run a small synthetic 1D demo
#' @param ns number of spatial sites
#' @param Tlen number of time points
#' @param m_nngp neighbors
#' @param seed RNG seed
run_demo_1d <- function(ns = 200, Tlen = 100, m_nngp = 15, seed = 305) {
  set.seed(seed)
  s  <- seq_len(ns)
  D  <- as.matrix(stats::dist(cbind(1, s)))
  
  matern_cov <- function(x, range = 5, kappa = 2.5, sigma2 = 1) {
    geoR::matern(x, phi = range, kappa = kappa, sigma2 = sigma2)
  }
  
  # latent advection/diffusion priors
  C_lat  <- matern_cov(D, range = 25, kappa = 2.5, sigma2 = 1)
  cholC  <- chol(Matrix::Matrix(C_lat, sparse = FALSE))
  mu     <- as.numeric(cholC %*% stats::rnorm(ns)) / 10
  mu[2:3] <- c(-1.26, -1.29)
  sig    <- sqrt(abs(exp(cholC %*% stats::rnorm(ns)))) / 5
  sig[1:3] <- c(0.333, 0.359, 0.387)
  
  # simulate state & data
  M <- stats::dnorm(D, mean = mu, sd = sig * 5)
  M <- M / rowSums(M)
  diag(M) <- diag(M) - 0.05
  C_eps   <- matern_cov(D, range = 5, kappa = 2.5, sigma2 = 1)
  cholE   <- chol(Matrix::Matrix(C_eps, sparse = FALSE))
  
  Y <- matrix(0, ns, Tlen)
  Y[, 1] <- stats::rnorm(ns) / 1000 + cholE %*% stats::rnorm(ns)
  for (t in 2:Tlen) {
    Y[, t] <- M %*% Y[, t - 1] + (1/3) * (cholE %*% stats::rnorm(ns)) + (1/5) * stats::rnorm(ns)
  }
  
  # Build model pieces
  Ginv <- make_Ginv_1d(mu = rep(0, ns), sig = rep(1, ns) / 4)
  Gtinv <- rep(list(Ginv), Tlen)
  Ft    <- rep(list(Matrix::Diagonal(ns)), Tlen)
  
  # Prior precision via NNGP
  N5    <- make_neighbors(D, m = 5)
  mycov <- function(x) matern_cov(x, range = 3, kappa = 1.5, sigma2 = 1)
  C0inv <- make_Cinv(D, N5, mycov, ns)
  
  # Hyper-precisions for mu, sig
  N15     <- make_neighbors(D, m = m_nngp)
  mycov_h <- function(x) matern_cov(x, range = 5, kappa = 2.5, sigma2 = 1)
  Sig0inv <- 100 * make_Cinv(D, N15, mycov_h, ns)
  
  # Settings
  delta <- 0.5
  v     <- 0.85
  
  # KF + FFBS
  kf <- sparse_kf_ffbs(Y, m0 = rep(0, ns), C0inv = C0inv, Ft = Ft, Gtinv = Gtinv, v = v, delta = delta)
  
  # Sufficient stats for sig
  cts_s <- make_ct_sig(mu, sig, samp = kf$draw, Y)
  Dts_s <- make_Dt_sig(mu, sig, samp = kf$draw, Y)
  bt_s  <- cts_s / Dts_s
  sig_new <- samp_sig(kf$Ctinv, C0inv, Dts_s, delta, bt_s, mu0 = rep(0, ns), Sig0inv = 10 * Sig0inv, sig_fallback = sig)
  
  # Sufficient stats for mu
  cts_m <- make_ct_mu(mu, sig_new, samp = kf$draw, Y)
  Dts_m <- make_Dt_mu(mu, sig_new, samp = kf$draw, Y)
  bt_m  <- cts_m / Dts_m
  mu_new <- samp_mu(kf$Ctinv, C0inv, Dts_m, delta, bt_m, mu0 = rep(0, ns), Sig0inv = Sig0inv)
  
  list(Y = Y, draw = kf$draw, mu = mu_new, sig = sig_new)
}

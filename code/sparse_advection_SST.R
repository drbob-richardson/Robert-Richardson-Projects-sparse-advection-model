#' SST 2D Sparse Implied–Advection Utilities
#'
#' 2D preprocessing + model helpers.
#' Requires: Matrix, geoR, fields (optional plotting), tmvtnorm, ncdf4
#'

# ---- packages ----
req <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. Please install it.", pkg), call. = FALSE)
  }
}
req("Matrix"); req("geoR"); req("tmvtnorm"); req("ncdf4")
# fields is only needed for quilt.plot demos:
fields_ok <- requireNamespace("fields", quietly = TRUE)

Diagonal <- Matrix::Diagonal

# ---- NetCDF helpers ----

#' Read a 3D SST block (lon x lat x time) from NetCDF
#' @param nc_path path to NetCDF file with variable (e.g., "sst")
#' @param var_name variable name in NetCDF (default "sst" or first 3D var)
#' @param lon_idx integer vector of longitude indices to keep
#' @param lat_idx integer vector of latitude indices to keep
#' @return list(data=array(lon x lat x time), lon_idx, lat_idx)
read_nc_block <- function(nc_path, var_name = NULL, lon_idx, lat_idx) {
  nc <- ncdf4::nc_open(nc_path)
  on.exit(ncdf4::nc_close(nc))
  
  # pick a default 3D var if not provided
  if (is.null(var_name)) {
    var_name <- names(nc$var)[sapply(nc$var, function(v) length(v$dim) == 3)][1]
  }
  if (is.na(var_name)) stop("No 3D variable found; please set var_name.")
  
  v <- nc$var[[var_name]]
  if (length(v$dim) != 3) stop("Selected variable is not 3D.")
  
  # dims assumed (lon, lat, time) — adjust if needed
  lon_len <- v$dim[[1]]$len
  lat_len <- v$dim[[2]]$len
  time_len <- v$dim[[3]]$len
  
  lon_idx <- lon_idx[lon_idx >= 1 & lon_idx <= lon_len]
  lat_idx <- lat_idx[lat_idx >= 1 & lat_idx <= lat_len]
  
  # read block
  start <- c(min(lon_idx), min(lat_idx), 1)
  count <- c(length(lon_idx), length(lat_idx), time_len)
  arr <- ncdf4::ncvar_get(nc, varid = var_name, start = start, count = count)
  
  list(data = arr, lon_idx = lon_idx, lat_idx = lat_idx)
}

#' Read and vectorize a land/sea mask with same lon/lat window
#' @param nc_path path to mask NetCDF
#' @param var_name mask var (1=ocean/keep, 0=land/drop) — adapt to your file
#' @param lon_idx,lat_idx same window as the SST block
#' @return logical vector length (lon*lat) indicating kept locations
read_mask_block <- function(nc_path, var_name, lon_idx, lat_idx) {
  nc <- ncdf4::nc_open(nc_path)
  on.exit(ncdf4::nc_close(nc))
  v <- nc$var[[var_name]]
  if (length(v$dim) != 2) stop("Mask variable should be 2D (lon x lat).")
  
  start <- c(min(lon_idx), min(lat_idx))
  count <- c(length(lon_idx), length(lat_idx))
  mask2d <- ncdf4::ncvar_get(nc, varid = var_name, start = start, count = count)
  
  # You may need to flip/transpose depending on your file
  as.vector(t(mask2d)) == 1
}

# ---- Anomalies by month ----

#' Compute monthly anomalies (remove month-wise climatology)
#' @param mat matrix (locations x time)
#' @return matrix anomalies (same dim)
monthly_anomalies <- function(mat) {
  L <- nrow(mat); TT <- ncol(mat)
  if (TT %% 12 != 0) warning("Time length not multiple of 12; using available cycles.")
  anoms <- mat
  for (m in 1:12) {
    cols <- seq(m, TT, by = 12)
    mm   <- rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE)
    anoms[, cols] <- mat[, cols, drop = FALSE] - mm
  }
  anoms
}

# ---- Grid / latlon utils ----

#' Build (lon, lat) grid matrix for a lon-by-lat block
#' @param lon_idx,lat_idx integer index vectors
#' @return matrix with 2 columns: (lon, lat) per flattened location
grid_latlon <- function(lon_idx, lat_idx) {
  # match your original indexing: lon=120:300, lat=60:120 (then remapped to 30:-30)
  # Here we just store the integer indices; map to degrees outside if needed
  as.matrix(expand.grid(lon = lon_idx, lat = rev(lat_idx)))[, 1:2]
}

#' Apply a region mask (logical predicate over lon/lat)
#' @param latlon matrix (n x 2), columns (lon, lat)
#' @param predicate function(lon, lat) -> logical
#' @return logical vector keep
apply_region_mask <- function(latlon, predicate) {
  predicate(latlon[, 1], latlon[, 2])
}

# ---- Boundary detection (4-neighborhood on the integer grid) ----

#' Find boundary neighbor placeholders (for values outside the kept set)
#' @param latlon kept grid (lon, lat) integers
#' @return matrix of boundary (lon, lat) integer pairs
find_boundaries <- function(latlon) {
  L <- nrow(latlon)
  uniq_lon <- sort(unique(latlon[, 1]))
  bound <- matrix(NA_integer_, 0, 2)
  for (lon in uniq_lon) {
    init <- latlon[latlon[, 1] == lon, 2]
    left <- latlon[latlon[, 1] == (lon - 1), 2]
    right <- latlon[latlon[, 1] == (lon + 1), 2]
    
    above <- init + 1
    below <- init - 1
    
    if (sum(!(init %in% left))  > 0) bound <- rbind(bound, cbind(lon - 1L, init[!(init %in% left)]))
    if (sum(!(init %in% right)) > 0) bound <- rbind(bound, cbind(lon + 1L, init[!(init %in% right)]))
    if (sum(!(init %in% above)) > 0) bound <- rbind(bound, cbind(lon,       init[!(init %in% above)] - 1L))
    if (sum(!(init %in% below)) > 0) bound <- rbind(bound, cbind(lon,       init[!(init %in% below)] + 1L))
  }
  unique(bound)
}

# ---- 9-point neighbor index table Glocs ----

#' Build Glocs (self + 8 nbh) indices by (lon, lat) adjacency on the kept set
#' @param latlon kept (lon, lat), in integer grid coordinates
#' @return integer matrix (n x 9): columns are [self, E, W, N, S, NE, NW, SE, SW]
build_Glocs <- function(latlon) {
  n <- nrow(latlon)
  Glocs <- matrix(0L, n, 9)
  Glocs[, 1] <- seq_len(n)
  
  # helper: find index of a target (lon,lat) if it exists, else fallback to self
  find_or_self <- function(i, dlon, dlat) {
    target <- which(latlon[, 1] == latlon[i, 1] + dlon & latlon[, 2] == latlon[i, 2] + dlat)
    if (length(target) > 0) target[1] else i
  }
  for (i in 1:n) {
    Glocs[i, 2] <- find_or_self(i, +1,  0) # E
    Glocs[i, 3] <- find_or_self(i, -1,  0) # W
    Glocs[i, 4] <- find_or_self(i,  0, +1) # N
    Glocs[i, 5] <- find_or_self(i,  0, -1) # S
    Glocs[i, 6] <- find_or_self(i, +1, +1) # NE
    Glocs[i, 7] <- find_or_self(i, -1, +1) # NW
    Glocs[i, 8] <- find_or_self(i, +1, -1) # SE
    Glocs[i, 9] <- find_or_self(i, -1, -1) # SW
  }
  Glocs
}

# ---- 2D G^{-1} builder (implied advection) ----

#' Construct sparse inverse evolution G^{-1} for 2D (implied advection)
#' @param mu1,mu2 advection fields (length n) in lon/lat directions
#' @param sig1,sig2 diffusion components (length n), positive
#' @param sig12 cross term (length n)
#' @param Glocs 9-col neighbor index table (from build_Glocs)
#' @return sparse n x n matrix
make_Ginv_2d <- function(mu1, mu2, sig1, sig2, sig12, Glocs) {
  k <- length(mu1)
  G <- Matrix::Matrix(0, k, k, sparse = TRUE)
  G[cbind(1:k, Glocs[, 1])] <- 1 + sig1 + sig2
  G[cbind(1:k, Glocs[, 2])] <- G[cbind(1:k, Glocs[, 2])] - mu1/2 - sig1/2
  G[cbind(1:k, Glocs[, 3])] <- G[cbind(1:k, Glocs[, 3])] + mu1/2 - sig1/2
  G[cbind(1:k, Glocs[, 4])] <- G[cbind(1:k, Glocs[, 4])] - mu2/2 - sig2/2
  G[cbind(1:k, Glocs[, 5])] <- G[cbind(1:k, Glocs[, 5])] + mu2/2 - sig2/2
  G[cbind(1:k, Glocs[, 6])] <- G[cbind(1:k, Glocs[, 6])] - sig12/2
  G[cbind(1:k, Glocs[, 7])] <- G[cbind(1:k, Glocs[, 7])] + sig12/2
  G[cbind(1:k, Glocs[, 8])] <- G[cbind(1:k, Glocs[, 8])] + sig12/2
  G[cbind(1:k, Glocs[, 9])] <- G[cbind(1:k, Glocs[, 9])] - sig12/2
  G
}

# ---- bt helpers (fixing 1:(T+1) style indexing) ----

make_bt_mu1 <- function(mu1, mu2, sig1, sig2, sig12, samp, Glocs) {
  Tlen <- ncol(samp) - 1L
  cts <-   samp[, 1:Tlen] +
    (-1 - sig1 - sig2) * samp[Glocs[, 1], 2:(Tlen + 1)] +
    0.5 * sig1 * samp[Glocs[, 2], 2:(Tlen + 1)] +
    0.5 * sig1 * samp[Glocs[, 3], 2:(Tlen + 1)] +
    (mu2/2 + sig2/2) * samp[Glocs[, 4], 2:(Tlen + 1)] +
    (-mu2/2 + sig2/2) * samp[Glocs[, 5], 2:(Tlen + 1)] +
    0.5 * sig12 * samp[Glocs[, 6], 2:(Tlen + 1)] -
    0.5 * sig12 * samp[Glocs[, 7], 2:(Tlen + 1)] -
    0.5 * sig12 * samp[Glocs[, 8], 2:(Tlen + 1)] +
    0.5 * sig12 * samp[Glocs[, 9], 2:(Tlen + 1)]
  Dts <- -0.5 * samp[Glocs[, 2], 2:(Tlen + 1)] + 0.5 * samp[Glocs[, 3], 2:(Tlen + 1)]
  list(bt = cts / Dts, Dts = Dts)
}

make_bt_mu2 <- function(mu1, mu2, sig1, sig2, sig12, samp, Glocs) {
  Tlen <- ncol(samp) - 1L
  cts <-   samp[, 1:Tlen] +
    (-1 - sig1 - sig2) * samp[Glocs[, 1], 2:(Tlen + 1)] +
    (sig1/2 + mu1/2) * samp[Glocs[, 2], 2:(Tlen + 1)] +
    (0.5 * sig1 - mu1/2) * samp[Glocs[, 3], 2:(Tlen + 1)] +
    (sig2/2) * samp[Glocs[, 4], 2:(Tlen + 1)] +
    (sig2/2) * samp[Glocs[, 5], 2:(Tlen + 1)] +
    0.5 * sig12 * samp[Glocs[, 6], 2:(Tlen + 1)] -
    0.5 * sig12 * samp[Glocs[, 7], 2:(Tlen + 1)] -
    0.5 * sig12 * samp[Glocs[, 8], 2:(Tlen + 1)] +
    0.5 * sig12 * samp[Glocs[, 9], 2:(Tlen + 1)]
  Dts <- -0.5 * samp[Glocs[, 4], 2:(Tlen + 1)] + 0.5 * samp[Glocs[, 5], 2:(Tlen + 1)]
  list(bt = cts / Dts, Dts = Dts)
}

make_bt_sig1 <- function(mu1, mu2, sig1, sig2, sig12, samp, Glocs) {
  Tlen <- ncol(samp) - 1L
  cts <-   samp[, 1:Tlen] +
    (-1 - sig2)  * samp[Glocs[, 1], 2:(Tlen + 1)] +
    (mu1/2)      * samp[Glocs[, 2], 2:(Tlen + 1)] +
    (-mu1/2)     * samp[Glocs[, 3], 2:(Tlen + 1)] +
    (sig2/2+mu2/2) * samp[Glocs[, 4], 2:(Tlen + 1)] +
    (sig2/2-mu2/2) * samp[Glocs[, 5], 2:(Tlen + 1)] +
    0.5 * sig12  * samp[Glocs[, 6], 2:(Tlen + 1)] -
    0.5 * sig12  * samp[Glocs[, 7], 2:(Tlen + 1)] -
    0.5 * sig12  * samp[Glocs[, 8], 2:(Tlen + 1)] +
    0.5 * sig12  * samp[Glocs[, 9], 2:(Tlen + 1)]
  Dts <- samp[Glocs[, 1], 2:(Tlen + 1)] - 0.5 * samp[Glocs[, 2], 2:(Tlen + 1)] - 0.5 * samp[Glocs[, 3], 2:(Tlen + 1)]
  list(bt = cts / Dts, Dts = Dts)
}

make_bt_sig2 <- function(mu1, mu2, sig1, sig2, sig12, samp, Glocs) {
  Tlen <- ncol(samp) - 1L
  cts <-   samp[, 1:Tlen] +
    (-1 - sig1)  * samp[Glocs[, 1], 2:(Tlen + 1)] +
    (sig1/2+mu1/2) * samp[Glocs[, 2], 2:(Tlen + 1)] +
    (0.5*sig1-mu1/2) * samp[Glocs[, 3], 2:(Tlen + 1)] +
    (mu2/2)      * samp[Glocs[, 4], 2:(Tlen + 1)] +
    (-mu2/2)     * samp[Glocs[, 5], 2:(Tlen + 1)] +
    0.5 * sig12  * samp[Glocs[, 6], 2:(Tlen + 1)] -
    0.5 * sig12  * samp[Glocs[, 7], 2:(Tlen + 1)] -
    0.5 * sig12  * samp[Glocs[, 8], 2:(Tlen + 1)] +
    0.5 * sig12  * samp[Glocs[, 9], 2:(Tlen + 1)]
  Dts <- samp[Glocs[, 1], 2:(Tlen + 1)] - 0.5 * samp[Glocs[, 4], 2:(Tlen + 1)] - 0.5 * samp[Glocs[, 5], 2:(Tlen + 1)]
  list(bt = cts / Dts, Dts = Dts)
}

make_bt_sig12 <- function(mu1, mu2, sig1, sig2, sig12, samp, Glocs) {
  Tlen <- ncol(samp) - 1L
  cts <-   samp[, 1:Tlen] +
    (-1 - sig1 - sig2) * samp[Glocs[, 1], 2:(Tlen + 1)] +
    (sig1/2 + mu1/2)   * samp[Glocs[, 2], 2:(Tlen + 1)] +
    (0.5*sig1 - mu1/2) * samp[Glocs[, 3], 2:(Tlen + 1)] +
    (sig2/2 + mu2/2)   * samp[Glocs[, 4], 2:(Tlen + 1)] +
    (sig2/2 - mu2/2)   * samp[Glocs[, 5], 2:(Tlen + 1)]
  Dts <-  -0.5 * sig12 * samp[Glocs[, 6], 2:(Tlen + 1)] +
    0.5 * sig12 * samp[Glocs[, 7], 2:(Tlen + 1)] +
    0.5 * sig12 * samp[Glocs[, 8], 2:(Tlen + 1)] -
    0.5 * sig12 * samp[Glocs[, 9], 2:(Tlen + 1)]
  list(bt = cts / Dts, Dts = Dts)
}

# ---- Conjugate / truncated Gaussian parameter updates ----

samp_mu_vec <- function(Ctinv_list, C0inv, Dts, delta, bt, mu0, Sig0inv) {
  k <- nrow(Dts); Tlen <- ncol(Dts)
  Vpostinv <- Sig0inv
  bpost    <- Sig0inv %*% mu0
  for (t in 1:Tlen) {
    Dmat <- Diagonal(x = as.numeric(Dts[, t]))
    Prec <- if (t == 1) C0inv else Ctinv_list[[t - 1]]
    scal <- delta / (1 - delta)
    Vpostinv <- Vpostinv + (Dmat %*% Prec %*% Dmat) * scal
    bpost    <- bpost    + (Dmat %*% Prec %*% Dmat) %*% (bt[, t] * scal)
  }
  drop(Matrix::solve(Vpostinv, bpost, sparse = TRUE))
}

samp_sig_pos <- function(Ctinv_list, C0inv, Dts, delta, bt, mu0, Sig0inv, fallback) {
  k <- nrow(Dts); Tlen <- ncol(Dts)
  Vpostinv <- Sig0inv
  bpost    <- Sig0inv %*% mu0
  for (t in 1:Tlen) {
    Dmat <- Diagonal(x = as.numeric(Dts[, t]))
    Prec <- if (t == 1) C0inv else Ctinv_list[[t - 1]]
    scal <- delta / (1 - delta)
    Vpostinv <- Vpostinv + (Dmat %*% Prec %*% Dmat) * scal
    bpost    <- bpost    + (Dmat %*% Prec %*% Dmat) %*% (bt[, t] * scal)
  }
  mn <- Matrix::solve(Vpostinv, bpost, sparse = TRUE)
  out <- tmvtnorm::rtmvnorm.sparseMatrix(
    n = 1, mean = as.numeric(mn), H = Vpostinv,
    lower = rep(0, k), upper = rep(Inf, k)
  )
  if (any(is.na(out))) out <- fallback
  as.numeric(out)
}

samp_sig12_bounded <- function(Ctinv_list, C0inv, Dts, delta, bt, mu0, Sig0inv, sig1, sig2, fallback) {
  k <- nrow(Dts); Tlen <- ncol(Dts)
  Vpostinv <- Sig0inv
  bpost    <- Sig0inv %*% mu0
  for (t in 1:Tlen) {
    Dmat <- Diagonal(x = as.numeric(Dts[, t]))
    Prec <- if (t == 1) C0inv else Ctinv_list[[t - 1]]
    scal <- delta / (1 - delta)
    Vpostinv <- Vpostinv + (Dmat %*% Prec %*% Dmat) * scal
    bpost    <- bpost    + (Dmat %*% Prec %*% Dmat) %*% (bt[, t] * scal)
  }
  mn <- Matrix::solve(Vpostinv, bpost, sparse = TRUE)
  bound <- sqrt(pmax(0, sig1 * sig2))
  out <- tmvtnorm::rtmvnorm.sparseMatrix(
    n = 1, mean = as.numeric(mn), H = Vpostinv,
    lower = -bound, upper = bound
  )
  if (any(is.na(out))) out <- fallback
  as.numeric(out)
}

# ---- Example driver (preprocessing + structure) ----

#' Example: Load SST, compute anomalies, mask region, and build 2D structures
#' @param sst_nc path to SST NetCDF
#' @param sst_var variable name (default auto)
#' @param mask_nc path to land/sea mask NetCDF
#' @param mask_var variable name in mask file
#' @param lon_idx,lat_idx integer indices to read (e.g., 120:300, 60:120)
#' @param time_range optional integer vector of time indices to keep
#' @param region_pred predicate(lon,lat)->logical to keep (default: equatorial window)
#' @return list with latlon, anoms, dat, Glocs, boundaries, mask_keep
sst_2d_prepare <- function(
    sst_nc, sst_var = NULL,
    mask_nc, mask_var,
    lon_idx, lat_idx,
    time_range = NULL,
    region_pred = function(lon, lat) (lon >= 200 & lon <= 260 & lat >= -10 & lat <= 10)
) {
  sst <- read_nc_block(sst_nc, var_name = sst_var, lon_idx = lon_idx, lat_idx = lat_idx)
  arr <- sst$data  # lon x lat x time
  
  # vectorize to (locations x time) in (lon-major, then lat)
  L  <- length(sst$lon_idx) * length(sst$lat_idx)
  TT <- dim(arr)[3]
  mat <- matrix(0, L, TT)
  k <- 1L
  for (i in seq_len(dim(arr)[1])) {
    for (j in seq_len(dim(arr)[2])) {
      mat[k, ] <- arr[i, j, ]
      k <- k + 1L
    }
  }
  
  anoms <- monthly_anomalies(mat)
  latlon <- grid_latlon(sst$lon_idx, sst$lat_idx)
  
  # mask by land/sea
  keep_ocean <- read_mask_block(mask_nc, var_name = mask_var, lon_idx = sst$lon_idx, lat_idx = sst$lat_idx)
  latlon <- latlon[keep_ocean, , drop = FALSE]
  anoms  <- anoms[keep_ocean, , drop = FALSE]
  dat    <- mat   [keep_ocean, , drop = FALSE]
  
  if (!is.null(time_range)) {
    anoms <- anoms[, time_range, drop = FALSE]
    dat   <- dat  [, time_range, drop = FALSE]
  }
  
  # region mask (replace with your original shapes if desired)
  keep_region <- apply_region_mask(latlon, region_pred)
  latlon <- latlon[keep_region, , drop = FALSE]
  anoms  <- anoms[keep_region, , drop = FALSE]
  dat    <- dat  [keep_region, , drop = FALSE]
  
  # optional quick plot
  if (fields_ok) {
    fields::quilt.plot(latlon[, 1], latlon[, 2], dat[, 1])
  }
  
  # boundary placeholders (if you use them downstream)
  boundaries <- find_boundaries(latlon)
  
  # neighbor index table for 9-point stencil
  Glocs <- build_Glocs(latlon)
  
  list(latlon = latlon, anoms = anoms, dat = dat,
       Glocs = Glocs, boundaries = boundaries,
       mask_keep = which(keep_ocean & keep_region))
}

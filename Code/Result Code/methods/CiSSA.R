library(Rssa)
library(signal)
library(gsignal)

# Вспомогательные функиции 
mse <- function (f_true, f_reconstructed) 
{
  mean((f_true - f_reconstructed)^2)
}


dftmtx <- function (n) 
{
  y <- stats::mvfft(diag(1, n))
  y
}


diag_averaging <- function (A) 
{
  B <- Re(A[nrow(A):1, ])
  as.numeric(lapply(split(B, -(row(B) - col(B))), mean))
}

shift_vector <- function (vec) 
{
  last_element <- tail(vec, 1)
  vec <- vec[-length(vec)]
  shifted_vec <- c(last_element, vec)
  return(shifted_vec)
}



extend <- function (x, H) 
{
  N <- length(x)
  p <- floor(N/3)
  dx <- diff(x)
  A <- aryule(dx, p)$a
  y <- x
  dy <- diff(y)
  er <- signal::filter(A, 1, dy)
  dy <- signal::filter(1, A, c(er, rep(0, H)))
  y <- y[1] + c(0, cumsum(dy))
  y <- rev(y)
  dy <- diff(y)
  er <- signal::filter(A, 1, dy)
  dy <- signal::filter(1, A, c(er, rep(0, H)))
  y <- y[1] + c(0, cumsum(dy))
  xe <- rev(y)
  xe
}



# Основная функиця для Circular SSA
circulant_SSA <- function (ts, L = NULL, extend_flag = FALSE) 
{
  time_series <- ts
  N <- length(time_series)
  if (is.null(L)) {
    L <- (N + 1)%/%2
  }
  if (extend_flag == FALSE) {
    H <- 0
    time_series <- ts
  }
  else {
    H <- L
    time_series <- extend(ts, H)
  }
  X <- hankel(time_series, L)
  if (L%%2) {
    nf2 <- (L + 1)/2 - 1
  }
  else {
    nf2 <- L/2 - 1
  }
  nft <- nf2 + abs((L%%2) - 2)
  autocov <- numeric(L)
  for (m in 0:(L - 1)) {
    autocov[[m + 1]] <- sum(time_series[1:(N - m)] * time_series[(1 + 
                                                                    m):N])/(N - m)
  }
  circ_first_row <- numeric(L)
  for (m in 0:(L - 1)) {
    circ_first_row[[m + 1]] <- (L - m)/L * autocov[[m + 1]] + 
      (m)/L * autocov[[L - m]]
  }
  S_C <- matrix(circ_first_row, nrow = 1)
  shifted_vector <- circ_first_row
  for (i in 2:(L)) {
    shifted_vector <- shift_vector(shifted_vector)
    S_C <- rbind(as.vector(shifted_vector), S_C)
  }
  U <- dftmtx(L)/sqrt(L)
  U[, 1] <- Re(U[, 1])
  for (k in 1:nf2) {
    u_k <- U[, k + 1]
    U[, k + 1] <- sqrt(2) * Re(u_k)
    U[, L + 2 - (k + 1)] <- sqrt(2) * Im(u_k)
  }
  if (L%%2 != 0) {
    U[, nft] <- Re(U[, nft])
  }
  psd <- abs(diag(t(U) %*% S_C %*% U))
  W <- t(U) %*% X
  R <- matrix(0, nrow = N + 2 * H, ncol = L)
  for (k in 1:L) {
    R[, k] <- diag_averaging(U[, k] %*% t(W[k, ]))
  }
  Z <- matrix(0, nrow = N + 2 * H, ncol = nft)
  Z[, 1] <- R[, 1]
  imp <- numeric(nft)
  lambda_sm <- sum(psd)
  imp[1] <- psd[1]/lambda_sm
  for (k in 1:nf2) {
    Z[, k + 1] <- R[, k + 1] + R[, L + 2 - (k + 1)]
    imp[k + 1] <- (psd[k + 1] + psd[L + 2 - (k + 1)])/lambda_sm
  }
  if (L%%2 != 0) {
    Z[, nft] <- R[, nft]
    imp[nft] <- psd[nft]/lambda_sm
  }
  list(t_series = Z[(H + 1):(N + H), ], importance = imp, freq = (0:(length(imp) - 
                                                                       1))/L)
}


# группировка по частотам
grouping_cissa <- function (cissa_res, groups) 
{
  freq <- cissa_res$freq
  t_series <- cissa_res$t_series
  residuals <- 0
  result <- setNames(as.list(rep(0, length(groups))), names(groups))
  result_freqs <- list()
  for (i in 1:length(cissa_res$freq)) {
    flag <- FALSE
    for (name in names(groups)) {
      if (flag == TRUE) {
        break
      }
      if (groups[[name]][1] <= freq[i] & freq[i] <= groups[[name]][2]) {
        flag <- TRUE
        result[[name]] <- result[[name]] + t_series[, 
                                                    i]
        result_freqs[[name]] <- c(result_freqs[[name]], 
                                  freq[i])
      }
    }
    if (flag == FALSE) {
      residuals <- residuals + t_series[, i]
    }
  }
  result[["residuals"]] <- residuals
  return(list(t_series = result, freqs_by_group = result_freqs))
  result
}


# Разложение Фурье
reconstruct_fft <- function (x_init, y_init, extend_flag = FALSE) 
{
  x <- x_init
  y <- y_init
  N <- length(y_init)
  H <- 0
  if (extend_flag == TRUE) {
    H <- length(y)%/%2
    y <- extend(y, H)
    x <- 0:(length(y) - 1)
  }
  frequencies <- (0:(length(x) - 1))/length(x)
  fft_y <- fft(y)
  amplitudes <- Mod(fft_y)
  phases <- Arg(fft_y)
  reconstructed <- matrix(0, length(amplitudes), length(x))
  n <- length(amplitudes)
  L <- n
  for (i in 1:(length(amplitudes))) {
    reconstructed[i, ] <- amplitudes[i] * cos(2 * pi * frequencies[i] * 
                                                (x) + phases[i])/n
  }
  group_by_elementary_freq_foureir <- function(res_averaged) {
    nf2 <- 0
    if (L%%2) {
      nf2 <- (L + 1)/2 - 1
    }
    else {
      nf2 <- L/2 - 1
    }
    nft <- nf2 + abs((L%%2) - 2)
    Z <- matrix(0, ncol = nft, nrow = n)
    Z[, 1] <- res_averaged[1, ]
    for (k in 1:nf2) {
      Z[, k + 1] <- res_averaged[k + 1, ] + res_averaged[L + 
                                                           2 - (k + 1), ]
    }
    if (L%%2 != 0) {
      Z[, nft] <- res_averaged[nft, ]
    }
    return(list(t_series = Z[(H + 1):(N + H), ], freq = (0:(dim(Z)[2] - 
                                                              1))/L))
  }
  rs <- group_by_elementary_freq_foureir(reconstructed)
  return(list(t_series = rs$t_series, freq = rs$freq))
}
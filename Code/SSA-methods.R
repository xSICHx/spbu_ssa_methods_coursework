# GSSA

library(Rssa)


afc <- function(filter, omega) {
  k <- seq_along(filter) - 1 #значения от 0 до length(filter) - 1
  h <- function(o) sum(rev(filter) * exp(-k*1i * o))
  abs(sapply(omega, h))
}


gssa <- function(ts, L, alpha = 0){
  N <- length(ts)
  tmp <- 1:(L)
  w <- sqrt(sin(pi*tmp/(L+1))^2)**alpha
  X <- diag(w) %*% hankel(ts, L)
  K <- dim(X)[2] 
  X_t <- t(X)
  S <- X %*% X_t
  
  tmp <- eigen(S)
  lambda <- tmp$values
  lambda <- lambda[lambda > 0]
  d <- length(lambda)
  U <- tmp$vectors[, 1:d]
  
  
  get_V_i <- function(i){
    X_t %*% U[, i] / sqrt(lambda[i])
  }
  
  get_X_i <- function(i){
    sqrt(lambda[i]) * as.matrix(U[, i]) %*% t(as.matrix(get_V_i(i)))
  }
  
  X_mats <- lapply(1:d, get_X_i)
  
  create_A_alpha <- function(){
    A <- c()
    for (i in 1:N){
      if (i <= L-1){
        A <- c(A, sum(w[1:i]))
      }
      else if (L <= i & i <= K){
        A <- c(A, sum(w))
      }
      else{
        A <- c(A, sum(w[(i-K+1):L]))
      }
    }
    A
  }
  
  A_alpha <- create_A_alpha()
  
  # A_alpha |> print()
  
  
  mean_w <- function(x){
    sum(x)/A_alpha[length(x)]
  }
  
  diag_averaging_w <- function(A){
    B <- A[nrow(A):1, ] |> Re()
    lapply(split(B, -(row(B) - col(B)) ), mean_w) |> as.numeric()
  }
  
  reconstr_elementary <- function(){
    lapply(X_mats, diag_averaging_w)
  }
  
  r_elem <- reconstr_elementary()
  # r_elem |> print()
  
  return(
    list(
      ts = r_elem,
      lambda = lambda
    )
  )
  
}




get_middle_point_filters <- function(ts, L, alpha = 0){
  N <- length(ts)
  tmp <- 1:(L)
  w <- sqrt(sin(pi*tmp/(L+1))^2)^alpha
  X <- diag(w) %*% hankel(ts, L)
  K <- dim(X)[2] 
  X_t <- t(X)
  S <- X %*% X_t
  
  tmp <- eigen(S)
  lambda <- tmp$values
  lambda <- lambda[lambda > 0]
  d <- length(lambda)
  U <- tmp$vectors[, 1:d]
  
  weight <- sum(w)
  
  get_filter_i <- function(i){
    filt <- c()
    u <- U[, i]
    
    get_sum_j <- function(j){
      sm <- 0
      for ( k in 1:(L-abs(j)) ){
        sm <- sm + u[k]*u[k+abs(j)]*w[k]
      }
      sm / weight
    }
    
    lapply( (-(L-1)):(L-1), get_sum_j ) |> unlist()
  }
  
  lapply(1:d, get_filter_i)
  
  
}









# CISSA
library(Rssa)
library(signal)
library(gsignal)


dftmtx <- function(n) {
  y <- stats::mvfft(diag(1, n))
  y
}

diag_averaging <- function(A){
  B <- A[nrow(A):1, ] |> Re()
  lapply(split(B, -(row(B) - col(B)) ), mean) |> as.numeric()
}

shift_vector <- function(vec) {
  last_element <- tail(vec, 1)
  vec <- vec[-length(vec)]
  shifted_vec <- c(last_element, vec)
  return(shifted_vec)
}

extend <- function(x, H){
  # Вычисление коэффициентов AR модели для дифференцированного ряда
  N <- length(x)
  p <- floor(N / 3)
  dx <- diff(x)
  # A <- ar(dx, order.max = p, method = "yule-walker")$ar
  A <- aryule(dx, p)$a
  
  # Правое расширение
  y <- x
  dy <- diff(y)
  er <- signal::filter(A, 1, dy)
  dy <- signal::filter(1, A, c(er, rep(0, H)))
  y <- y[1] + c(0, cumsum(dy))
  
  # Левое расширение
  y <- rev(y)
  dy <- diff(y)
  er <- signal::filter(A,1,dy)
  dy <- signal::filter(1,A,c(er, rep(0, H)))
  y <- y[1] + c(0, cumsum(dy))
  
  # Расширенный ряд
  xe <- rev(y)
  
  # Вывод результатов
  xe 
}


circulant_SSA <- function(ts, L = NULL, extend_flag = FALSE){
  time_series <- ts
  # Construct trajectory matrix
  N <- length(time_series)
  if (is.null(L)){
    L <- (N + 1)%/%2
  }
  # Проверка на расширения ряда
  if (extend_flag == FALSE){
    H <- 0
    time_series <- ts
  }
  else{
    H <- L
    time_series <- extend(ts, H)
  }
  
  X <- hankel(time_series, L)
  
  # Number of symmetric frequency pairs around 1/2
  if (L %% 2) {
    nf2 <- (L + 1) / 2 - 1
  } else {
    nf2 <- L / 2 - 1
  }
  
  # Number of frequencies <= 1/2
  nft <- nf2 + abs((L %% 2) - 2)
  
  # Decomposition
  # Estimate autocovariance     OK
  autocov <- numeric(L)
  for (m in 0:(L-1)){
    autocov[[m+1]] <- sum(time_series[1:(N-m)] * time_series[(1+m):N]) / (N-m)
  }
  
  # First row of circulant matrix
  circ_first_row <- numeric(L)
  for (m in 0:(L-1)){
    circ_first_row[[m+1]] <- (L-m)/L * autocov[[m+1]] + (m)/L * autocov[[L-m]]
  }
  
  # Build circulant matrix
  S_C <- matrix(circ_first_row, nrow = 1)
  shifted_vector <- circ_first_row
  for (i in 2:(L)) {
    shifted_vector <- shift_vector(shifted_vector)
    # S_C <- rbind(S_C, as.vector(shifted_vector))
    S_C <- rbind(as.vector(shifted_vector), S_C)
  }
  
  # Eigenvectors of circulant matrix (unitary base)
  U <- dftmtx(L)/sqrt(L)
  
  # Real eigenvectors (orthonormal base)
  U[, 1] <- Re(U[, 1])
  for (k in 1:nf2) {
    u_k <- U[, k + 1]
    U[, k + 1] <- sqrt(2) * Re(u_k)
    U[, L + 2 - (k + 1)] <- sqrt(2) * Im(u_k)
  }
  if (L %% 2 != 0) {
    U[, nft] <- Re(U[, nft])
  }
  
  # Eigenvalues of circulant matrix: estimated power spectral density
  psd <- abs(diag(t(U) %*% S_C %*% U))
  
  # Principal components
  W <- t(U) %*% X
  # Reconstruction
  # Elementary reconstructed series
  R <- matrix(0, nrow = N+2*H, ncol = L)
  for (k in 1:L) {
    R[, k] <- U[ ,k] %*% t(W[k, ]) |> diag_averaging()
  }
  
  # Grouping by frequency
  # Elementary reconstructed series by frequency
  Z <- matrix(0, nrow = N+2*H, ncol = nft)
  Z[, 1] <- R[, 1]
  # Importance of component
  imp <- numeric(nft)
  lambda_sm <- sum(psd)
  imp[1] <- psd[1]/lambda_sm
  for (k in 1:nf2) {
    Z[, k + 1] <- R[, k + 1] + R[, L + 2 - (k + 1)]
    imp[k+1] <- (psd[k+1] + psd[ L + 2 - (k + 1)])/lambda_sm
  }
  if (L %% 2 != 0) {
    Z[, nft] <- R[, nft]
    imp[nft] <- psd[nft] / lambda_sm
  }
  
  list(t_series = Z[(H+1):(N+H),],
       importance = imp,
       freq = (0:(length(imp) -1))/L
  )
}


# groups - list of frequencies
grouping_cissa <- function(cissa_res, groups){
  freq <- cissa_res$freq
  t_series <- cissa_res$t_series
  
  residuals <- 0
  result <- setNames(as.list(rep(0, length(groups))), names(groups))
  for (i in 1:length(cissa_res$freq)){
    flag <- FALSE
    for (name in names(groups)){
      if (groups[[name]][1] <= freq[i] & freq[i] <= groups[[name]][2]){
        flag <- TRUE
        result[[name]] <- result[[name]] + t_series[, i]
      }
    }
    
    if (flag == FALSE){
      residuals <- residuals + t_series[, i]
    }
  }
  
  result[["residuals"]] <- residuals
  result
}
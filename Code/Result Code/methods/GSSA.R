library(Rssa)

# Основной алгоритм gssa
gssa <- function(ts, # значения временного ряда
                 L, # длина окна
                 alpha = 0 # параметр для весов
                 ){
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


# Получить фильтр для средней точки 
get_middle_point_filters <- function(ts, # временной ряд
                                     L,  # длина окна
                                     alpha = 0 # параметр для весов 
                                     ){
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


# Функция фильтра для любой точки
get_all_point_filters <- function(ts, # значения временного ряда
                                  L,  # длина окна
                                  alpha = 0 # параметр для весов 
                                  ){
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
    
    W <- matrix(0, nrow = K, ncol = length(ts))
    W_w <- matrix(0, nrow = K, ncol = length(ts))
    for (i in 1:K){
      W[i, i:(i+L-1)] <- u * w
      W_w[i, i:(i+L-1)] <- u 
    }
    W_w <- t(W_w)
    
    D <- c()
    cm_sm <- 0
    for (i in 1:L){
      cm_sm <- cm_sm + w[i]
      D <- c(D, cm_sm)
    }
    for (i in (L+1):K){
      D<- c(D, cm_sm)
    }
    for (i in (K+1):length(ts)){
      cm_sm <- cm_sm - w[i-K]
      D <- c(D, cm_sm)
    }
    
    D <- diag(D)
    D <- solve(D)
    
    # print((D %*% W_w %*% W) %*% ts)
    
    (D %*% W_w %*% W)
    
  }
  lapply(1:d, get_filter_i)
}





# Вспомогательные функции 

# АЧХ линейного фильтра для каждой точки
afc <- function(filter, omega) {
  k <- seq_along(filter) - 1 #значения от 0 до length(filter) - 1
  h <- function(o) sum(rev(filter) * exp(-k*1i * o))
  abs(sapply(omega, h))
}

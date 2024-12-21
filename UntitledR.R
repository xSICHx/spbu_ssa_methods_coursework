
afc <- function(filter, omega) {
  k <- seq_along(filter) - 1 #значения от 0 до length(filter) - 1
  h <- function(o) sum(rev(filter) * exp(-k*1i * o))
  abs(sapply(omega, h))
}

### WOSSA


library(Rssa)






n <- 96*2-1
x <- 0:(n-1)
L <- 96
y1 <- sin(2*pi/12 * x)
y2 <- cos(2*pi/3 * x)/2
y3 <- exp(x/100) + 1
y <- y1 + y2 + y3
alpha <- 0.5
filt <- sqrt(sin(pi * 1:(L) / (L + 1))^2)^alpha


s <- ssa(y, L)
r <- reconstruct(s, groups = list(
  e = c(1),
  s1 = 2:3,
  s2 = 4:5
))


sw <- ssa(y, L, row.oblique = filt)
rw <- reconstruct(sw, groups = list(
  e = c(1),
  s1 = 2:3,
  s2 = 4:5
))

mse <- function(f_true, f_reconstructed){
  mean((f_true - f_reconstructed)^2) 
}

mse(r$e, rw$e)
mse(r$e, y3)
mse(rw$e, y3)
# mse()

# plot(rw$e)
# lines(y3, col="red")


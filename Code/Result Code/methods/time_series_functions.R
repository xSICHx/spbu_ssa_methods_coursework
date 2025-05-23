f_const <- function (x, C = 0) 
{
  rep(C, length(x))
}


f_cos <- function (x, A = 1, omega = 1/4, phi = 0) 
{
  f_exp_mod_harm_series(x, A, alpha = 0, omega = omega, phi = phi)
}


f_exp <- function (x, A = 1, alpha = 1) 
{
  A * exp(alpha * x)
}


f_exp_cos <- function (x, A = 1, alpha = 1, omega = 1/4, phi = 0) 
{
  f_exp_mod_harm_series(x, A, alpha, omega, phi)
}


f_exp_mod_harm_series <- function (x, A = 1, alpha = 1, omega = 1/4, phi = 0) 
{
  A * exp(alpha * x) * cos(2 * pi * omega * x + phi)
}


f_linear <- function (x, a = 1, b = 0) 
{
  a * x + b
}


f_sin <- function (x, A = 1, omega = 1/4, phi = 3 * pi/2) 
{
  f_exp_mod_harm_series(x, A, alpha = 0, omega = omega, phi = phi)
}


f_sum <- function (x) 
{
  f_cos(x, omega = omega_cs) + f_exp_mod_harm_series(x, a = a, 
                                                     omega = omega_exp) + f_sin(x, omega = omega_sn)
}

generate_ts <- function (func, n = 1000, ...) 
{
  ts(func(1:n, ...))
}
cal_d <- function(trt_m, trt_sd, trt_n, ctl_m, ctl_sd, ctl_n){
  d <- (trt_m - ctl_m)/sqrt(((trt_n-1)*trt_sd^2 + (ctl_n-1)*ctl_sd^2)/(trt_n + ctl_n -2))
  v_d <- (trt_n + ctl_n)/(trt_n*ctl_n) + d^2 /(2*(trt_n + ctl_n))
  cbind(d, v_d)
}

cal_t <- function(n_trt, m_trt, sd_trt, m_ctl, sd_ctl, n_ctl){
  sp <- sqrt(((n_trt-1)*sd_trt^2 + (n_ctl-1)*sd_ctl^2)/(n_trt + n_ctl -2))
  t <- (m_trt - m_ctl)/(sp*sqrt((1/n_ctl)+ (1/n_ctl)))
  t
}

cal_rce <- function(d, n1, n2) {
  a <- (n1+n2)^2/(n1*n2)
  r_ce <- d/sqrt(d^2 + a^2)
  r_ce
}

# large sample variance of r_ce
cal_v_rce <- function(d, n1, n2){
  a <- (n1+n2)^2/(n1*n2)
  v_rce <- (a^4/(d^2 + a^2)^3)*((n1+n2)/(n1*n2) + (d^2)/(2*(n1+n2)))
  v_rce
}


## Calculate t to d - we may not need this. 
cal_t_d <- function(t, trt_n, ctl_n){
  d <- t*sqrt(1/trt_n + 1/ctl_n) 
  d
}

# d to r_b
dtob <- function(d, n1, n2){
  h <- (n1 + n2 -2)/n1 + (n1+n2-2)/n2
  r_pb <- d/sqrt(d^2 + h)
  p <- n1/(n1+n2)
  c <- qnorm(p)
  r_b <- (sqrt(p*(1-p))/dnorm(c))*r_pb
  return(r_b)
}

# t to r_b 
ttob <- function(t, n1, n2){
  m <- n1+n2-2
  r_pb <- t/sqrt(t^2 + m)
  p <- n1/(n1+n2)
  c <- qnorm(p)
  r_b <- (sqrt(p*(1-p))/dnorm(c))*r_pb
  return(r_b)
}

### Sampling variance of r_b - Soper's approxiamte method.

cal_v_rb <- function(r_b, n1, n2) {
  p <- n1/(n1+n2)
  c <- qnorm(p)
  v_rb <- (1/(n1+n2-1))*((sqrt(p*(1-p))/dnorm(c) - r_b^2)^2)
  return(v_rb)
}

# function from John Fox's polychor function
# This uses ML estimator. 

pcor <- function(x) {
  require(polycor)
  mat <- matrix(x, 2, 2)
  pc <- polychor(x = mat, ML = TRUE, std.err=TRUE)
  c(pcor = pc$rho, v_pcor = pc$var[1,1])
}

# function from ?
calc_tet <- function(a00, a11, a10, a01) {
  or <- a00 * a11 / (a01 * a10)
  v_or <- 1/a00 + 1/a01 + 1/a10 + 1/a11
  delta <- 1 + sqrt(a00 * a11 / (a10 * a01))
  r_tet <- cos(pi / delta) # from Pearson, 1900
  n <- a00 + a11 + a10 + a01
  p0 <- (a00 + a01) / n
  p1 <- (a00 + a10) / n
  z0 <- qnorm(p0)
  z1 <- qnorm(p1)
  h0 <- dnorm(z0)
  h1 <- dnorm(z1)
  v_tet <- (p0 * (1 - p0) * p1 * (1 - p1)) / (n * (h0 * h1)^2)
  return(cbind(or, v_or, r_tet, v_tet))
}

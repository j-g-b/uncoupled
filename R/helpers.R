#'
#'
#'
quantize_real_line <- function(V, sigma, n){
  alpha0 <- -(V + sigma)*log(n)
  N <- ceiling(2*(n^(1/4))*log(n))
  alpha <- seq(alpha0, alpha0+(N*(V+sigma))/(n^(1/4)), length.out = N+1)
  return(alpha[alpha < V & alpha > -V])
}
#'
#'
#'
update_alpha_masses <- function(alpha, mus, sigma){
  alpha_masses <- rep(0, length(alpha))
  for(i in 1:length(alpha)){
    if(i == 1){
      alpha_masses[i] <- pnorm(alpha[i + 1], mean = alpha, sd = sigma) %>%
        magrittr::multiply_by(mus) %>%
        sum()
    } else if(i == length(alpha)){
      alpha_masses[i] <- (1 - pnorm(alpha[i], mean = alpha, sd = sigma)) %>%
        magrittr::multiply_by(mus) %>%
        sum()
    } else {
      alpha_masses[i] <- (pnorm(alpha[i + 1], mean = alpha, sd = sigma) - pnorm(alpha[i], mean = alpha, sd = sigma)) %>%
        magrittr::multiply_by(mus) %>%
        sum()
    }
  }
  return(alpha_masses)
}
#'
#'
#'
simplex_project <- function(alpha){
  u <- sort(alpha, decreasing = T)
  rho <- plyr::aaply(1:length(u), 1, function(j){
    u[j] + (1/j)*(1 - sum(u[1:j]))
  }) %>%
    magrittr::is_greater_than(0) %>%
    which() %>%
    max()
  lambda <- (1/rho)*(1 - sum(u[1:rho]))
  return(pmax(alpha + lambda, 0))
}
#'
#'
#'
initialize_mus <- function(alpha, y, sigma){
  mus <- rep(0, length(alpha))
  for(i in 1:length(alpha)){
    mus[i] <- dnorm(y, mean = alpha[i]) %>% sum()
  }
  return(mus / sum(mus))
}
#'
#'
#'
quantile_func <- function(x, alpha, mus, w = NULL){
  n <- length(x)
  if(is.null(w)){
    w <- seq(0, 1, length.out = n)[-1]
  }
  g_x <- rep(NA, n)
  F_x <- cumsum(mus[order(alpha)])
  F_x[length(F_x)] <- 1
  alpha_sort <- sort(alpha)
  for(i in 1:length(g_x)){
    g_x[i] <- alpha_sort[min(which(F_x >= w[i]))]
  }
  return(g_x)
}
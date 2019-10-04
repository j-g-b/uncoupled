#'
#'
#'
min_wasserstein <- function(x, y, weights = NULL, V = NULL, cost_fun = NULL, sigma = NULL, maxit = 2000){
  #
  require(magrittr)
  require(CVXR)
  require(Matrix)
  #
  n <- length(y)
  #
  if(is.null(weights)){
    weights <- rep(1/n, n)
  }
  #
  alpha <- uncoupled::quantize_real_line(V, sigma, n)
  mus <- uncoupled::initialize_mus(alpha, y, sigma)
  alpha_masses <- uncoupled::update_alpha_masses(alpha, mus, sigma)
  #
  costs <- sapply(y, function(x){(x - alpha)^2}) %>% c()
  #
  alpha_sum_mat <- Matrix::sparseMatrix(i = rep(1:length(alpha), each = length(y)),
                                        j = sapply(1:length(alpha), function(i){i + length(alpha)*(0:(length(y)-1))}),
                                        x = 1,
                                        dims = c(length(alpha), length(y)*length(alpha)))
  y_sum_mat <- Matrix::sparseMatrix(i = rep(1:length(y), each = length(alpha)),
                                    j = sapply(1:length(y), function(i){((i-1)*length(alpha) + 1):(i*length(alpha))}) %>% c(),
                                    x = 1,
                                    dims = c(length(y), length(y)*length(alpha)))
  #
  grad_mus <- matrix(1, nrow = length(alpha), ncol = length(alpha))
  C <- CVXR::Variable(length(alpha)*n)
  objective <- CVXR::Minimize(t(costs) %*% C)
  iter <- 1
  best_iter <- 1
  best_value <- 1000
  best_mu <- mus
  s_prev <- 0.01*rep(1, length(mus))
  #
  cat("Optimizing with projected subgradient descent\n")
  #
  while(iter < maxit){
    # Compute Wasserstein-2 distance between discrete measures `alpha_masses` and `weights`
    constraints <- list(C >= 0, alpha_sum_mat %*% C == alpha_masses, y_sum_mat %*% C == weights)
    problem <- CVXR::Problem(objective, constraints)
    result <- CVXR::psolve(problem)
    # Update best objective value, `mus` value, and iter index
    if(result$value < best_value){
      best_value <- result$value
      best_mu <- mus
      best_iter <- iter
    } else {
      if((iter - best_iter) > 100 & iter > 500){
        break
      }
    }
    # Extract negative dual weights, the subgradient of the objective w.r.t. `alpha_masses`
    grad_alpha <- -1*result$getDualValue(constraints[[2]]) %>% c()
    # Calculate subgradient w.r.t. `mus`
    grad_mus <- matrix(0, nrow = length(alpha), ncol = length(alpha))
    for(i in 1:length(alpha)){
      if(i == 1){
        grad_mus[i, ] <- pnorm(alpha[i+1], mean = alpha, sd = sigma)
      } else if(i == length(alpha)){
        grad_mus[i, ] <- 1 - pnorm(alpha[i], mean = alpha, sd = sigma)
      } else {
        grad_mus[i, ] <- (pnorm(alpha[i+1], mean = alpha, sd = sigma) -
                            pnorm(alpha[i], mean = alpha, sd = sigma))
      }
    }
    grad_mus <- apply(grad_mus, 2, function(x){x*grad_alpha}) %>% colSums()
    # Update `mus` and `alpha_masses` with projected subgradient descent using `CFM` step size
    # Camerini P.M., Fratta L., Maffioli F. (1975) On improving relaxation methods by modified gradient techniques. In: Balinski M.L., Wolfe P. (eds) Nondifferentiable Optimization. Mathematical Programming Studies, vol 3. Springer, Berlin, Heidelberg
    gamma <- 1.5
    lam <- (maxit / 100) / j
    beta <- max(0, -1*gamma*sum(s_prev*grad_mus)/sum(s_prev^2))
    s_curr <- grad_mus + beta*s_prev
    s_prev <- s_curr
    step_size <- (result$value - best_value + lam) / sum(s_curr^2)
    mus <- uncoupled::simplex_project(mus - step_size*s_curr)
    alpha_masses <- uncoupled::update_alpha_masses(alpha, mus, sigma)
    # Update `iter`
    iter <- iter + 1
    # Print to console
    if(iter%%10 == 0){
      cat("\rvalue: ", sprintf("%f", result$value), "  iter: ", iter)
    }
  }
  cat("\n")
  #
  q_x <- uncoupled::quantile_func(x, alpha, best_mu, w = cumsum(weights))
  #
  return(q_x)
}
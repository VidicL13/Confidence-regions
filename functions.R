T_hotelling <- function(mu, col_means, critical, norm_data, n, p){
  n * (t(col_means - mu) %*% solve(cov(norm_data)) %*% (col_means - mu)) * ((n - p)/(p * (n - 1))) - critical
}

T_hotelling_radial <- function(r, phi, critical, norm_data, n, p){
  n * r^2 * (t(c(cos(phi), sin(phi))) %*% solve(cov(norm_data)) %*% (c(cos(phi), sin(phi)))) * ((n - p)/(p * (n - 1))) - critical
}

T_likelihood <- function(mu_x, mu_y, hat_rho, hat_sigma_x, hat_sigma_y, n, hat_mu_x, hat_mu_y, critical){
  (n/(1 - hat_rho^2)) * (((hat_mu_x - mu_x)/hat_sigma_x)^2 + ((hat_mu_y - mu_y)/hat_sigma_y)^2 - 2 * hat_rho * ((hat_mu_x - mu_x)/hat_sigma_x) * ((hat_mu_y - mu_y)/hat_sigma_y)) - critical
}

T_likelihood_phi <- function(r, phi, rho, sigma_x, sigma_y, n, critical){
  ((n) / (1 - rho**2)) * 
    (((r * cos(phi)) / (sigma_x))**2 +
       ((r * sin(phi)) / (sigma_y))**2 - 
       (2 * rho * ((r * cos(phi)) / (sigma_x)) * ((r * sin(phi)) / (sigma_y)))
    ) - critical
}
# 1D function
fun_likelihood <- function(x, critical, hat_mu, norm_data, n) {
  -2 * (sum(dnorm(
    norm_data,
    mean = x,
    sd = sqrt((sum((
      norm_data - x
    ) ^ 2)) / (n)),
    log = TRUE
  )) - sum(dnorm(
    norm_data,
    mean = hat_mu,
    sd = sqrt((sum((norm_data - hat_mu) ^ 2
    )) / (n)),
    log = TRUE
  ))) - critical
}

#### 1D ####
# Get all the data for 1d cases
cr_1d <- function(row){
  run <- row['run']
  n <- row['n']
  mu <- row['mu']
  sigma <- row['sigma']
  # Create a sample
  norm_data <- rnorm(n = n, mean = mu, sd = sigma)
  
  # normal procedure
  tic()
  s_n <- sigma / sqrt(n)
  mean_n <- mean(norm_data)
  xmin <- mean_n + qnorm(.975, mean = 0, sd = 1, lower.tail = FALSE) * s_n
  xmax <- mean_n + qnorm(.975, mean = 0, sd = 1, lower.tail = TRUE) * s_n
  elapsed <- toc(quiet = TRUE)
  elapsed <- elapsed$toc - elapsed$tic
  names(elapsed) <- NULL
  CI_n <- data.frame(xmin = xmin,
                     xmax = xmax,
                     method = 'Normal',
                     comp_time = elapsed)
  
  # Hotteling procedure in 1D is equal to student's t/test
  tic()
  s_h <- sd(norm_data) / sqrt(n)
  mean_n <- mean(norm_data)
  xmin <- mean_n + qt(.975, df = n-1, lower.tail = FALSE) * s_h
  xmax <- mean_n + qt(.975, df = n-1, lower.tail = TRUE) * s_h
  elapsed <- toc(quiet = TRUE)
  elapsed <- elapsed$toc - elapsed$tic
  names(elapsed) <- NULL
  CI_h <- data.frame(xmin = xmin,
                     xmax = xmax,
                     method = 'Hotelling',
                     comp_time = elapsed)
  
  # Log Likelihood ratio 
  tic()
  s_l <- sd(norm_data)
  mean_n <- mean(norm_data)
  xmin <- uniroot(fun_likelihood, interval = c(mean_n, mean_n-5*sigma), critical = qchisq(.95, df = 1), hat_mu = mean_n, norm_data = norm_data, n = n)$root
  xmax <- uniroot(fun_likelihood, interval = c(mean_n, mean_n+5*sigma), critical = qchisq(.95, df = 1), hat_mu = mean_n, norm_data = norm_data, n = n)$root
  elapsed <- toc(quiet = TRUE)
  elapsed <- elapsed$toc - elapsed$tic
  names(elapsed) <- NULL
  CI_w <- data.frame(xmin = xmin,
                     xmax = xmax,
                     method = 'Wilks',
                     comp_time = elapsed)
  
  CI <- bind_rows(CI_n, CI_h, CI_w)
  CI <- CI %>% mutate(ci_len = xmax - xmin,
                      is_in = (xmin <= mu & mu <= xmax) %>% as.numeric(),
                      run = run,
                      n = n,
                      sigma = sigma,
                      mu = mu)
  return(CI)
}





#### Hotelling ####
Hotelling_cr <- function(norm_data, p, n, len_out = 100, epsilon, mu, points = FALSE) {
  tic()
  col_means <- colMeans(norm_data)
  col_sd <- apply(norm_data, 2, sd)
  critical <- qf(0.95, df1 = p, df2 = n - p)
  i <- 0
  dimx <- 0
  dimy <- 0
  # search for the border until you have at least 3 points 
  while (dimx <= 4 | dimy <= 4) {
    # initial grid
    if (i == 0) {
      grid <- expand.grid(
        mux = seq(
          col_means[1] - 4 * col_sd[1]/sqrt(n),
          col_means[1] + 4 * col_sd[1]/sqrt(n),
          length.out = len_out
        ),
        muy = seq(
          col_means[2] - 4 * col_sd[2]/sqrt(n),
          col_means[2] + 4 * col_sd[2]/sqrt(n),
          length.out = len_out
        )
      )
    }
    # define grid based on T values
    else{
      min_t <- vred %>% arrange(desc(abs(T_h2))) %>% tail(10)
      minx <- min(min_t$mux) - 3 * col_sd[1]/sqrt(n)
      maxx <- max(min_t$mux) + 3 * col_sd[1]/sqrt(n)
      miny <- min(min_t$muy) - 3 * col_sd[2]/sqrt(n)
      maxy <- max(min_t$muy) + 3 * col_sd[2]/sqrt(n)
      grid <- expand.grid(
        mux = seq(minx,
                  maxx,
                  length.out = len_out),
        muy = seq(miny,
                  maxy,
                  length.out = len_out)
      )
    }
    # calculate T values
    vred <- grid %>%
      rowwise() %>%
      mutate(T_h2 = T_hotelling(
        mu = c(mux, muy),
        col_means = c(col_means),
        critical = critical,
        norm_data = norm_data,
        n = n,
        p = p
      )) %>%
      ungroup()
    
    # gather all that are in epsilon proximity
    vred_out <- vred %>%
      filter(abs(T_h2) <= epsilon) %>%
      dplyr::select(mux, muy)
    dimx <- length(unique(vred_out$mux))
    dimy <- length(unique(vred_out$muy))
    i <- i + 1
    if(i==10){
      return()
    }
  }
  # return points to plot
  if(points){
    vred_out <- vred_out %>% mutate(method = 'Hotelling grid')
    return(vred_out)
  }
  # Calculate metrics
  convhull_1 <- convhulln(vred_out, options = 'FA')
  elapsed_1 <- toc(quiet = TRUE)
  elapsed_1 <- elapsed_1$toc - elapsed_1$tic
  names(elapsed_1) <- NULL
  res1 <- data.frame(
    method = 'Hotelling grid',
    vol = convhull_1$vol,
    comp_time = elapsed_1,
    is_in = inhulln(convhull_1, p = matrix(mu, ncol = 2)),
    num_border_points = nrow(vred_out),
    num_recalc = i)
  return(res1)
}


#### Likelihood ####
Likelihood_cr <- function(norm_data, n, epsilon, mu, len_out = 100, points = FALSE){
  tic()
  critical <- qchisq(0.95, df = 2)
  col_means <- colMeans(norm_data)
  col_sd <- apply(norm_data, 2, sd)
  hat_rho <- cor(norm_data)[1,2]
  
  i <- 0
  dimx <- 0
  dimy <- 0
  # search for the border until you have at least 3 points 
  while (dimx <= 4 | dimy <= 4) {
    # initial grid
    if (i == 0) {
      grid <- expand.grid(
        mux = seq(
          col_means[1] - 4 * col_sd[1]/sqrt(n),
          col_means[1] + 4 * col_sd[1]/sqrt(n),
          length.out = len_out
        ),
        muy = seq(
          col_means[2] - 4 * col_sd[2]/sqrt(n),
          col_means[2] + 4 * col_sd[2]/sqrt(n),
          length.out = len_out
        )
      )
    }
    # define grid based on T values
    else{
      min_t <- vred %>% arrange(desc(abs(T_like))) %>% tail(10)
      minx <- min(min_t$mux) - 3 * col_sd[1]/sqrt(n)
      maxx <- max(min_t$mux) + 3 * col_sd[1]/sqrt(n)
      miny <- min(min_t$muy) - 3 * col_sd[2]/sqrt(n)
      maxy <- max(min_t$muy) + 3 * col_sd[2]/sqrt(n)
      grid <- expand.grid(
        mux = seq(minx,
                  maxx,
                  length.out = len_out),
        muy = seq(miny,
                  maxy,
                  length.out = len_out)
      )
    }
    vred <- grid %>% 
    rowwise() %>%
    mutate(T_like = T_likelihood(
      mu_x = mux,
      mu_y = muy,
      hat_rho = hat_rho,
      hat_mu_x = col_means[1],
      hat_mu_y = col_means[2],
      n = n,
      critical = critical,
      hat_sigma_x = col_sd[1],
      hat_sigma_y = col_sd[2]
    )) %>% 
    ungroup()
    # gather all that are in epsilon proximity
    vred_out <- vred %>%
      filter(abs(T_like) <= epsilon) %>%
      dplyr::select(mux, muy)
    dimx <- length(unique(vred_out$mux))
    dimy <- length(unique(vred_out$muy))
    i <- i + 1
    if(i==10){
      return()
    }
  }
  # return points to plot
  if(points){
    vred_out <- vred_out %>% mutate(method = 'Likelihood grid')
    return(vred_out)
  }
  # Calculate metrics
  convhull_1 <- convhulln(vred_out, options = 'FA')
  elapsed_1 <- toc(quiet = TRUE)
  elapsed_1 <- elapsed_1$toc - elapsed_1$tic
  names(elapsed_1) <- NULL
  res1 <- data.frame(
    method = 'Likelihood grid',
    vol = convhull_1$vol,
    comp_time = elapsed_1,
    is_in = inhulln(convhull_1, p = matrix(mu, ncol = 2)),
    num_border_points = nrow(vred_out),
    num_recalc = i)
  return(res1)
}

#### Hotelling radial ####
Hotelling_radial_cr <- function(n, p, norm_data, len_out=180, mu, points = FALSE){
  tic()
  critical <- qf(0.95, df1 = p, df2 = n - p)
  col_means <- colMeans(norm_data)
  
  vred <- data.frame(phi = seq(0, 2 * pi, length.out = len_out)) %>%
    rowwise() %>%
    mutate(r = uniroot(
      f = T_hotelling_radial,
      phi = phi,
      critical = critical,
      norm_data = norm_data,
      n = n,
      p = p,
      interval = c(0, 400)
    )$root) %>%
    ungroup() %>%
    transmute(mux = col_means[1] + r * cos(phi), muy = col_means[2] + r * sin(phi))
  # return points to plot
  if(points){
    vred_out <- vred %>% mutate(method = 'Hotelling radial')
    return(vred_out)
  }
  # compute metrics
  convhull_1 <- convhulln(vred, options = 'FA')
  elapsed_1 <- toc(quiet = TRUE)
  elapsed_1 <- elapsed_1$toc - elapsed_1$tic
  names(elapsed_1) <- NULL
  res1 <- data.frame(
    method = 'Hotelling radial',
    vol = convhull_1$vol,
    comp_time = elapsed_1,
    is_in = inhulln(convhull_1, p = matrix(mu, ncol = 2)),
    num_border_points = nrow(vred),
    num_recalc = 1)
  return(res1)
}

#### Likelihood radial ####
Likelihood_radial_cr <- function(n, norm_data, len_out = 180, mu, points = FALSE){
  tic()
  col_sd <- apply(norm_data, 2, sd)
  col_means <- colMeans(norm_data)
  hat_rho <- cor(norm_data)[1,2]
  critical <- qchisq(0.95, df = 2)
  vred <- data.frame(phi = seq(0, 2 * pi, length.out = len_out)) %>%
    rowwise() %>%
    mutate(r = uniroot(
      f = T_likelihood_phi,
      phi = phi,
      rho = hat_rho,
      sigma_x = col_sd[1],
      sigma_y = col_sd[2],
      n = n,
      critical = critical,
      interval = c(0,400)
    )$root) %>% 
    ungroup() %>% 
    transmute(mux = col_means[1] + r * cos(phi), muy = col_means[2] + r * sin(phi))
  # return points to plot
  if(points){
    vred_out <- vred %>% mutate(method = 'Likelihood radial')
    return(vred_out)
  }
  # compute metrics
  convhull_1 <- convhulln(vred, options = 'FA')
  elapsed_1 <- toc(quiet = TRUE)
  elapsed_1 <- elapsed_1$toc - elapsed_1$tic
  names(elapsed_1) <- NULL
  res1 <- data.frame(
    method = 'Likelihood radial',
    vol = convhull_1$vol,
    comp_time = elapsed_1,
    is_in = inhulln(convhull_1, p = matrix(mu, ncol = 2)),
    num_border_points = nrow(vred),
    num_recalc = 1)
  return(res1)
}

#### cr all ####
cr_2d <- function(row, points=FALSE){
  mu_x <- row['mu_x']
  mu_y <- row['mu_y']
  rho <- row['rho']
  sigma_x <- row['sigma_x']
  sigma_y <- row['sigma_y']
  n <- row['n']
  epsilon <- row['epsilon']
  run <- row['run']
  
  # get the data
  mu <- c(mu_x, mu_y)
  Sigma <- matrix(c(sigma_x^2, rho * sigma_x * sigma_y, rho * sigma_x * sigma_y, sigma_y ^2), nrow = 2)
  norm_data <- mvrnorm(n = n, mu = mu, Sigma = Sigma) %>% as.data.frame()
  
  p <- ncol(norm_data)
  
  res1 <- Hotelling_cr(norm_data = norm_data, p = p, n = n, epsilon = epsilon, mu = mu, points = points)
  res2 <- Hotelling_radial_cr(n = n, p = p, norm_data = norm_data, mu = mu, points = points)
  res3 <- Likelihood_cr(norm_data = norm_data, n = n, epsilon = epsilon, mu = mu, points = points)
  res4 <- Likelihood_radial_cr(n = n, norm_data = norm_data, mu = mu, points = points)
  
  CR <- bind_rows(res1, res2, res3, res4)
  CR <- CR %>% mutate(mu_x = mu_x,
                      mu_y = mu_y,
                      rho = rho,
                      sigma_x = sigma_x,
                      sigma_y = sigma_y,
                      n = n,
                      epsilon = epsilon,
                      run = run
  )
  val <- getTxtProgressBar(pb)
  setTxtProgressBar(pb, val+1)
  return(CR)
}


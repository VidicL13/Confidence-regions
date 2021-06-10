# Monte-Carlo simulation where we calculate: 
#   - volume, 
#   - computation time, 
#   - if the confidence region contains the theoretical value

mu <- c(0,0)
rho <- 0.5
sigma_x <- 1
sigma_y <- 1
n <- 10
n_runs <- 50
save_file <- TRUE
sigma <- matrix(c(sigma_x^2, rho * sigma_x * sigma_y, rho * sigma_x * sigma_y, sigma_y ^2), nrow = 2)
source('./libraries.r')
source('./functions.r')

results <- data.frame(
  method = character(),
  volume = numeric(),
  comupte_time = numeric(),
  is_inhul = logical(),
  run = integer()
)
pb <- txtProgressBar(min = 0, max = n_runs, style = 3)
for(i in 1:n_runs){
  binorm_data <- mvrnorm(n = n, mu = mu, Sigma = sigma) %>% as.data.frame()
  p <- ncol(binorm_data)
  
  # Hotelling method
  col_means <- colMeans(binorm_data)
  
  col_sd <- apply(binorm_data, 2, sd)
  critical <- qf(0.95, df1 = p, df2 = n - p)
  tic()
  vred <- expand.grid(mux = seq(-2*col_sd[1], 2*col_sd[1], length.out = 100),
                      muy = seq(-2*col_sd[2], 2*col_sd[2], length.out = 100)) %>% 
    rowwise() %>% 
    mutate(T_h2 = T_hotelling(
      mu = c(mux, muy),
      col_means = c(col_means),
      critical = critical,
      binorm_data = binorm_data,
      n = n,
      p = p
    )) %>% 
    ungroup() %>% 
    mutate(nula = abs(T_h2)<=0.1)
  
  convhull_1 <- convhulln(vred %>% filter(nula) %>% dplyr::select(mux, muy), options = 'FA')
  
  # Get the time required to compute the CR & SCI
  elapsed_1 <- toc(quiet = TRUE)
  elapsed_1 <- elapsed_1$toc - elapsed_1$tic
  names(elapsed_1) <- NULL
  
  res1 <- data.frame(
    method = 'Hotelling grid',
    volume = convhull_1$vol,
    comupte_time = elapsed_1,
    is_inhul = inhulln(convhull_1, p = matrix(mu, ncol = 2)),
    run = i
  )
  
  
  # Radial transform version of Hotelling method
  tic()
  vred2 <- data.frame(phi = seq(0, 2 * pi, length.out = 180)) %>%
    rowwise() %>%
    mutate(r = uniroot(
      f = T_hotelling_radial,
      phi = phi,
      critical = critical,
      binorm_data = binorm_data,
      n = n,
      p = p,
      interval = c(0, 20)
    )$root) %>%
    ungroup() %>%
    transmute(mux = col_means[1] + r * cos(phi), muy = col_means[2] + r * sin(phi))
  
  convhull_2 <- convhulln(vred2, options = 'FA')
  
  # Get the time required to compute the CR & SCI
  elapsed_2 <- toc(quiet = TRUE)
  elapsed_2 <- elapsed_2$toc - elapsed_2$tic
  names(elapsed_2) <- NULL
  
  res2 <- data.frame(
    method = 'Hotelling radial',
    volume = convhull_2$vol,
    comupte_time = elapsed_2,
    is_inhul = inhulln(convhull_2, p = matrix(mu, ncol = 2)),
    run = i
  )
  
  
  # Profile likelihood ratio radial version
  
  tic()
  vred3 <- confRegionBoundary(
    num_points = 180,
    alpha_value = 0.05,
    Sigma = sigma,
    mu_x = col_means[1],
    mu_y = col_means[2],
    n = n
  )
  
  convhull_3 <- convhulln(vred3$p, options = 'FA')
  
  # Get the time required to compute the CR & SCI
  elapsed_3 <- toc(quiet = TRUE)
  elapsed_3 <- elapsed_3$toc - elapsed_3$tic
  names(elapsed_3) <- NULL
  
  res3 <- data.frame(
    method = 'Likelihood ratio radial',
    volume = convhull_3$vol,
    comupte_time = elapsed_3,
    is_inhul = inhulln(convhull_3, p = matrix(mu, ncol = 2)),
    run = i
  )
  setTxtProgressBar(pb, i)
  results <- bind_rows(results, res1, res2, res3)
  
}

if(save_file){
  write.csv2(results, file = paste0('./MonteCarlo_', i, '_runs_', Sys.Date(), '.csv'), fileEncoding = 'UTF-8')
}


source('./libraries.R')

#### Functions ----
T_hotelling_radial <- function(r, phi, critical, norm_data, n, p){
  n * r^2 * (t(c(cos(phi), sin(phi))) %*% solve(cov(norm_data)) %*% (c(cos(phi), sin(phi)))) * ((n - p)/(p * (n - 1))) - critical
}

T_likelihood_phi <- function(r, phi, rho, sigma_x, sigma_y, n, critical){
  ((n) / (1 - rho**2)) * 
    (((r * cos(phi)) / (sigma_x))**2 +
       ((r * sin(phi)) / (sigma_y))**2 - 
       (2 * rho * ((r * cos(phi)) / (sigma_x)) * ((r * sin(phi)) / (sigma_y)))
    ) - critical
}

mvnormSample_2d <- function(mu_x = 0, mu_y = 0, n = 20, rho = 0.5, sigma_x = 1, sigma_y = 1){
  mu <- c(mu_x , mu_y)
  names(mu) <- c('mu_x', 'mu_y')
  Sigma <- matrix(c(sigma_x^2, 
                    sigma_x * sigma_y * rho, 
                    sigma_x * sigma_y * rho, 
                    sigma_y^2), 
                  2, 2, byrow = TRUE)
  sample <- mvrnorm(n = n, mu = mu, Sigma = Sigma) %>% as.data.frame()
  names(sample) <- c('x', 'y')
  return(sample)
}

Hotelling_radial_cr <- function(data, len_out=180, signf_lvl = 0.05){
  n <- nrow(data)
  p <- ncol(data)
  critical <- qf(1 - signf_lvl, df1 = p, df2 = n - p)
  col_means <- colMeans(data)
  
  vred <- data.frame(phi = seq(0, 2 * pi, length.out = len_out)) %>%
    rowwise() %>%
    mutate(r = uniroot(
      f = T_hotelling_radial,
      phi = phi,
      critical = critical,
      norm_data = data,
      n = n,
      p = p,
      interval = c(0, 400)
    )$root) %>%
    ungroup() %>%
    transmute(x = col_means[1] + r * cos(phi), y = col_means[2] + r * sin(phi))
  return(vred)
}

Likelihood_radial_cr <- function(data, len_out = 180, signf_lvl = 0.05){
  n = nrow(data)
  col_sd <- apply(data, 2, sd)
  col_means <- colMeans(data)
  hat_rho <- cor(data)[1,2]
  critical <- qchisq(1 - signf_lvl, df = 2)
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
    transmute(x = col_means[1] + r * cos(phi), y = col_means[2] + r * sin(phi))
  # return points to plot
  return(vred)
}

# 3D sample
mvnormSample_3d <- function(n = 20){
  mu <- c(0,0,0)
  names(mu) <- c('mu_x', 'mu_y', 'mu_z')
  Sigma <- matrix(c(1, 0.5, 0.5,
                    0.5, 1, 0.5,
                    0.5, 0.5, 1), 
                  3, 3, byrow = TRUE)
  sample <- mvrnorm(n = n, mu = mu, Sigma = Sigma) %>% as.data.frame()
  names(sample) <- c('x', 'y', 'z')
  return(sample)
}

# Hotelling 3D 
Hot_pivot <-  function(r, phi_1, phi_2, n, p, data, critical){
  n * r^2 * 
     t(c(
       cos(phi_1),
       sin(phi_1) * cos(phi_2),
       sin(phi_1) * sin(phi_2)
     )) %*%
     solve(cov(data)) %*% 
     c(
       cos(phi_1),
       sin(phi_1) * cos(phi_2),
       sin(phi_1) * sin(phi_2)
     ) *
     ((n - p)/(p * (n - 1))) - critical
}

Hot_cr_3d <- function(data, len_out=180, signf_lvl = 0.05){
  len_2 <- floor(sqrt(len_out/2))
  len_1 <- 2 * len_2
  n <- nrow(data)
  p <- ncol(data)
  critical <- qf(1-signf_lvl, df1 = p, df2 = n-p)
  col_means <- colMeans(data)
  phi_1 = seq(0, pi, length.out = len_1)
  phi_2 = seq(0, 2 * pi, length.out = len_2)
  vred <- expand_grid(phi_1, phi_2) %>%
    rowwise() %>%
    mutate(r = uniroot(
      f = Hot_pivot,
      phi_1 = phi_1,
      phi_2 = phi_2,
      critical = critical,
      data = data,
      n = n,
      p = p,
      interval = c(0, 400)
    )$root) %>%
    ungroup() %>%
    transmute(x = col_means[1] + r * cos(phi_1), y = col_means[2] + r * sin(phi_1) * cos(phi_2), z = col_means[3] + r * sin(phi_1) * sin(phi_2))
  return(vred)
}

# Likelihood  3D
ell_0 <- function(x, mu_hat, k, Sigma_inv, Sigma_det){
  -(1/2) * (
    t(x - mu_hat) %*%
      Sigma_inv %*%
      (x - mu_hat)
  ) -
    (k/2) * log(2* pi) - (1/2) * log(Sigma_det)
}

ell <- function(r, x, phi_1, phi_2, mu_hat, k, Sigma_inv, Sigma_det){
  -(1/2) * (
    t(x - mu_hat - r * 
        c(
          cos(phi_1),
          sin(phi_1) * cos(phi_2),
          sin(phi_1) * sin(phi_2)
        )
      ) %*%
      Sigma_inv %*%
      (x - mu_hat - r * 
         c(
           cos(phi_1),
           sin(phi_1) * cos(phi_2),
           sin(phi_1) * sin(phi_2)
         )
      )
    ) -
    (k/2) * log(2* pi) - (1/2) * log(Sigma_det)
}

Lik_pivot <- function(r, phi_1, phi_2, mu_hat, k, critical, Sigma_inv, Sigma_det, data){
  
  ells <- data %>%
    rowwise() %>% 
    transmute(
      var_ell = ell(
        r = r,
        x = c(x, y, z),
        phi_1 = phi_1,
        phi_2 = phi_2,
        mu_hat = mu_hat,
        k = k,
        Sigma_inv = Sigma_inv,
        Sigma_det = Sigma_det
      ),
      var_ell_0 = ell_0(
        x = c(x, y, z),
        mu_hat = mu_hat,
        k = k,
        Sigma_inv = Sigma_inv,
        Sigma_det = Sigma_det
      )
    ) %>% 
    ungroup() %>%
    summarise(sum_ell = sum(var_ell),
              sum_ell_0 = sum(var_ell_0))
  value <- -2 * (ells$sum_ell - ells$sum_ell_0) - critical
  return(value)
}

Lik_cr_3d <- function(data, len_out=180, signf_lvl=0.05){
  k <- ncol(data)
  n <- nrow(data)
  cov_mtrx <- cov(data)
  Sigma_inv <- solve(cov_mtrx)
  Sigma_det <- det(cov_mtrx)
  critical <- qchisq(1-signf_lvl, df = k)
  mu_hat <- colMeans(data)
  len_2 <- floor(sqrt(len_out/2))
  len_1 <- 2 * len_2
  
  phi_1 = seq(0, pi, length.out = len_1)
  phi_2 = seq(0, 2 * pi, length.out = len_2)
  vred <- expand_grid(phi_1, phi_2) %>%
    rowwise() %>%
    mutate(r = uniroot(
      f = Lik_pivot,
      phi_1 = phi_1,
      phi_2 = phi_2,
      mu_hat = mu_hat,
      k = k,
      critical = critical,
      Sigma_inv = Sigma_inv,
      Sigma_det = Sigma_det,
      data = data,
      interval = c(0, 400)
    )$root) %>%
    ungroup() %>%
    transmute(x = mu_hat[1] + r * cos(phi_1), y = mu_hat[2] + r * sin(phi_1) * cos(phi_2), z = mu_hat[3] + r * sin(phi_1) * sin(phi_2))
  return(vred)
  
}

#### 2D example ----
# Generate 2D data


data_2d <- mvnormSample_2d(n = n)
mu_hat <- colMeans(data_2d)
cr_2d_hot <- Hotelling_radial_cr(data = data_2d, len_out = 180) %>% slice(chull(x, y))
cr_2d_lik <- Likelihood_radial_cr(data = data_2d, len_out = 180) %>% slice(chull(x, y))

data_2d %>% ggplot(aes(x = x, y = y)) +
  geom_point() + 
  geom_point(aes(x = mu_hat[1], y = mu_hat[2]), size = 3, alpha = .5, color = 'red') +
  geom_polygon(data = cr_2d_hot, alpha = 0.3, fill='darkgreen') +
  geom_polygon(data = cr_2d_lik, alpha = 0.3, fill='yellow')


#### 3D example ----
cr_plot <- function(data, colour){
  chul_3d <- convhulln(data, options = 'FA')
  mesh_3d <- mesh3d(chul_3d$p)
  mesh_3d$it <- chul_3d$hull %>% t()
  mesh_3d$ib <- NULL
  wire3d(mesh_3d, col = colour, alpha = 0.4)
  shade3d(mesh_3d, col = colour, alpha = 0.1, override = T)
}

run_3d <- function(n_pts, n_bnry_pts){
  
  data_3d <- mvnormSample_3d(n = n_pts)
  mean_3d <- colMeans(data_3d) %>% t()
  
  # confidence regions
  data_3d_cr_hot <-
    Hot_cr_3d(data = data_3d,
              len_out = n_bnry_pts,
              signf_lvl = 0.05)
  data_3d_cr_lik <-
    Lik_cr_3d(data = data_3d,
              len_out = n_bnry_pts,
              signf_lvl = 0.05)
  
  # plot points, mu and estimates
  points3d(data_3d, alpha = 0.5)
  axes3d(c('x', 'y', 'z'))
  points3d(data_3d,
           size = 20,
           alpha = 0.2,
           col = 'grey')
  points3d(mean_3d, size = 3, col = 'red')
  points3d(mean_3d,
           size = 20,
           alpha = 0.3,
           col = 'red')
  
  # original
  points3d(data.frame(x = 0, y = 0, z = 0), size = 3, col = 'orange')
  points3d(
    data.frame(x = 0, y = 0, z = 0),
    size = 20,
    alpha = 0.3,
    col = 'orange'
  )
  
  # CR
  cr_plot(data = data_3d_cr_hot, colour = 'green')
  cr_plot(data = data_3d_cr_lik, colour = 'blue')
}
run_3d(n_pts = 300, n_bnry_pts = 180)
clear3d()




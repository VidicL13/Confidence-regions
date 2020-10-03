
#' Profile log likelihood ratio test statistic
#' 
#' @param n Number of observations in original sample
#' @param phi Polar coordinate
#' @param r Distance from central point
#' @param rho Correlation coefficient
#' @param sigma_x Standard deviation of x-coordinate
#' @param sigma_y Standard deviation of y-coordinate
#' 
T_phi <- function(r, phi, rho, sigma_x, sigma_y, n){
  ((n) / (1 - rho**2)) * 
    (((r * cos(phi)) / (sigma_x))**2 +
       ((r * sin(phi)) / (sigma_y))**2 - 
       (2 * rho * ((r * cos(phi)) / (sigma_x)) * ((r * sin(phi)) / (sigma_y)))
    )
}


#' Compute points on the confidence region boundary
#' 
#' @description The function computes the (1 - alpha) confidence region for bivariate normal distribution
#' based on Chi squared distribution and profile log likelihood ratio
#' 
#' @param alpha_value Significance level of confidence region.
#' @param mu_x Mean value of x-coordinate
#' @param mu_y Mean value of y-coordinate
#' @param n Number of observations in original sample
#' @param num_points Number of points to be determined.
#' @param rho Correlation coefficient
#' @param sigma_x Standard deviation of x-coordinate
#' @param sigma_y Standard deviation of y-coordinate
#' @param sample Number used to mark the sample number
#' @return list of 2: time and dataframe with cartesian coordinates of x & y
confRegionBoundary <- function(num_points, alpha_value, sigma_x, sigma_y, mu_x, mu_y, rho, n, sample = NULL){
  # set critical region and angles 
  Chisq_critical <- qchisq(1- alpha_value, 2)
  phi_seq <- seq(0, 2*pi, length.out = num_points)
  # start timer
  tic()
  # calculate the points
  region_points <- lapply(phi_seq, function(phi){
    # prepare the function to be optimized
    opti <- function(x){
      T_phi(r = x, phi = phi, rho = rho, sigma_x = sigma_x, sigma_y = sigma_y, n = n) - Chisq_critical
    }
    # find optimum value
    r <- uniroot(f = opti, interval = c(0,max(sigma_x, sigma_y)))$root
    # reparameterization
    result <- data.frame('x' = mu_x + r * cos(phi), 'y' = mu_y + r * sin(phi))
    return(result)
  }) %>% bind_rows()
  
  # end timer and calculate computation time
  time <- toc(quiet = TRUE)
  time <- time$toc - time$tic
  if(!is.null(sample)){
    region_points$sample <- c(as.character(sample))
  }
  results <- list(time = time, p = region_points)
  return(results)
}

#' Generate bivariate normal sample
#' 
#' @param mu_x Mean value of x-coordinate
#' @param mu_y Mean value of y-coordinate
#' @param n Number of observations in original sample
#' @param rho Correlation coefficient
#' @param sigma_x Standard deviation of x-coordinate
#' @param sigma_y Standard deviation of y-coordinate
mvnormSample <- function(mu_x, mu_y, n, rho, sigma_x, sigma_y){
  mu <- c(mu_x , mu_y)
  names(mu) <- c('mu_x', 'mu_y')
  Sigma <- matrix(c(sigma_x^2, 
                    sigma_x * sigma_y * rho, 
                    sigma_x * sigma_y * rho, 
                    sigma_y^2), 
                  2, 2, byrow = TRUE)
  sample <- mvrnorm(n = n, mu = mu, Sigma = Sigma) %>% as.data.frame()
  names(sample) <- c('x', 'y')
  return(sample = sample)
}

#' Generate MLE's
#' 
#' @param n Number of observations in original sample
#' @param sample Sample dataset
#' @return A list with: sample and all MLEs
mvnormMLE <- function(n, sample){
  
  # Compute MLE for mean values
  mu_x_MLE <- sum(sample$x) / n
  mu_y_MLE <- sum(sample$y) / n
  
  # Compute MLE for Sigma matrix
  Sigma_MLE <- apply(sample, 1, function(row){
    (row - c(mu_x_MLE, mu_y_MLE)) %*% t(row - c(mu_x_MLE, mu_y_MLE))
  }) %>% rowSums() / n
  
  # Compute MLE for standard deviation
  sigma_x_MLE <- sqrt(Sigma_MLE[1])
  sigma_y_MLE <- sqrt(Sigma_MLE[4])
  
  # Compute MLE for correlation coefficient
  rho_MLE <- Sigma_MLE[3] / (sigma_x_MLE * sigma_y_MLE)
  return(list(sample = sample,
              Sigma_MLE = Sigma_MLE,
              mu_x_MLE = mu_x_MLE,
              mu_y_MLE = mu_y_MLE,
              sigma_x_MLE = sigma_x_MLE,
              sigma_y_MLE = sigma_y_MLE,
              rho_MLE =  rho_MLE))
}





























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
#' @return time:  time required to compute the points on the boundary of confidence region
#'         label: label given for the sample 
#'         p:     data.frame with cartesian coordinates for confidence region boundary
confRegionBoundary <- function(num_points, alpha_value, sigma_x, sigma_y, mu_x, mu_y, rho, n, label = NULL){
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
    r <- uniroot(f = opti, interval = c(0,max(sigma_x, sigma_y, 20)))$root
    # reparameterization
    result <- data.frame('x' = mu_x + r * cos(phi), 'y' = mu_y + r * sin(phi))
    return(result)
  }) %>% bind_rows()
  
  # end timer and calculate computation time
  time <- toc(quiet = TRUE)
  time <- time$toc - time$tic
  names(time) <- NULL
  results <- list(time = time, p = region_points, label = label)
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
#' @param sample Sample dataset
#' @return Sigma_MLE:     2x2 correlation matrix
#'         mu_x_MLE:      MLE of mean value of x-coordinate of the sample
#'         mu_y_MLE:      MLE of mean value of y-coordinate of the sample
#'         sigma_x_MLE:   MLE of standard deviation of x-coordinate of the sample
#'         sigma_y_MLE:   MLE of standard deviation of y-coordinate of the sample
#'         rho_MLE:       MLE of correlation coefficient from the sample
mvnormMLE <- function(sample){
  n = nrow(sample)
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
  return(list(Sigma_MLE = Sigma_MLE,
              mu_x_MLE = mu_x_MLE,
              mu_y_MLE = mu_y_MLE,
              sigma_x_MLE = sigma_x_MLE,
              sigma_y_MLE = sigma_y_MLE,
              rho_MLE =  rho_MLE))
}

#' Compute metrics via 'method' for given sample
#' 
#' @param alpha_value Number: Significance level of confidence region.
#' @param method Character: Method to be used in \cset. See details for available methods
#' @param mu_x Mean value of x-coordinate
#' @param mu_y Mean value of y-coordinate
#' @param optimal_region Points of 'optimal' confidence region
#' @param sample Sample of points on which we want to calculate confidence region
#' @details Available methods for confidence regions are: boot.kern for the nonparametric bootstrap method
#' using kernel density estimation described in Pallmann & Jaki (2017); emp.bayes for the empirical 
#' Bayes region described in Casella & Hwang (1983); hotelling for the Hotelling-type region
#' described in Wang et al (1999); limacon.asy for the limacon-shaped mimimum expected volume
#' region described in Brown et al (1995); limacon.fin for the finite-sample variant of the minimum
#' expected volume region described in Berger & Hsu (1996); standard.cor for the standard region
#' incorporating correlation between parameters described in Chew (1966); standard.ind for the
#' standard region ignoring correlation between parameters; tost for the two one-sided test (TOST)
#' intervals described in Schuirmann (1987); tseng for the mimimum expected interval length region
#' described in Tseng (2002); tseng.brown for the pseudo-empirical Bayes region described in Tseng
#' & Brown (1997). Available methods for confidence intervals are: expanded for the two one-sided 
#' est (TOST) procedure (Schuirmann 1987) using the expanded intervals described e.g., in Bofinger (1992) and Hsu
#' et al. (1994); fix.seq for the fixed sequence intervals described in Maurer et al (1995) and Hsu &
#' Berger (1999); tost for the two one-sided test (TOST) intervals described in Schuirmann (1987).
#' See also an overview and comparison of all methods in Pallmann & Jaki (2017).
#' @return results:      a data.frame object with desired metrics
#'         CR_convhull:  convex hull for confidence region
#'         SCI_convhull: convex hull for simultaneous confidence intervals
metrics <- function(sample, 
                    method, 
                    mu_x, 
                    mu_y, 
                    alpha_value,
                    optimal_region,
                    sample_num){
  
  ####---------------------------------------------------------------------------------------------
  #
  # Calculate confidence set (CR & SCI) for our sample
  tic()
  confidence_set <- cset(dat = sample, method = method, alpha = alpha_value)
  elapsed <- toc(quiet = TRUE)
  
  # Get the time required to compute the CR & SCI
  elapsed <- elapsed$toc - elapsed$tic
  names(elapsed) <- NULL
  
  # check if confidence intervals have no limits
  if(any(confidence_set$ci == Inf)){
    return(NULL)
  }
  n <- nrow(sample)
  #### Simultaneous Confidence Intervals ----------------------------------------------------------
  #
  # get corner points for simultaneous confidence intervals
  SCI <- expand.grid(confidence_set$ci[1,], 
                     confidence_set$ci[2,],
                     KEEP.OUT.ATTRS = FALSE) %>% as.data.frame()
  names(SCI) <- c('x', 'y')
  
  # Compute the intersection & convex hulls
  SCI_intersection <- intersectn(SCI, optimal_region)
  SCI_is_inhull <- inhulln(SCI_intersection$ch1, matrix(c(mu_x, mu_y), 1, 2))
  
  # Gather the results
  SCI_results <- list(
    'computation time' = elapsed,
    'volume' = SCI_intersection$ch1$vol, 
    'intersection volume' = SCI_intersection$ch$vol,
    'theoretical volume' = SCI_intersection$ch2$vol, 
    'method' = paste0(method, '_SCI'),
    'is_inhull' = SCI_is_inhull,
    'convex hull' = SCI_intersection$ch1,
    'sample number' = sample_num,
    'sample size' = n
  )
  
  #### Confidence Regions -------------------------------------------------------------------------
  #
  # Some methods don't produce confidence regions, return only SCI
  if(is.null(confidence_set$cr)){
    return(list(SCI = SCI_results))
  }
  
  # Get border points of confidence regions
  CR <- confidence_set$cr %>% as.data.frame()
  names(CR) <- c('x', 'y')
  
  # Compute the intersection & convex hulls
  CR_intersection <- intersectn(CR, optimal_region)
  CR_is_inhull <- inhulln(CR_intersection$ch1, matrix(c(mu_x, mu_y), 1, 2))
  
  # Gather the results
  CR_results <- list(
    'computation time' = elapsed,
    'volume' = CR_intersection$ch1$vol, 
    'intersection volume' = CR_intersection$ch$vol,
    'theoretical volume' = CR_intersection$ch2$vol, 
    'method' = paste0(method, '_CR'),
    'is_inhull' = CR_is_inhull,
    'convex hull' = CR_intersection$ch1,
    'sample number' = sample_num,
    'sample size' = n
  )
  
  return(list(SCI = SCI_results, CR = CR_results))
}

#' Generate sample, compute theoretical CR, MLE approximate CR, CR & SCI from jocre package
#' 
#' @param alpha_value Number: Significance level of confidence region.
#' @param methods list: Methods to be used in \cset. See details for available methods
#' @param mu_x Mean value of x-coordinate
#' @param mu_y Mean value of y-coordinate
#' @param n Number of observations in original sample
#' @param num_points Number of points to be computed for theoretical CR & MLE approximate CR
#' @param rho Correlation coefficient
#' @param sample_num Number: Number or label of the sample in MC simulation
#' @param sigma_x Standard deviation of x-coordinate
#' @param sigma_y Standard deviation of y-coordinate
#' @param sample data.frame: sample
#' 
#' @details Available methods for confidence regions are: boot.kern for the nonparametric bootstrap method
#' using kernel density estimation described in Pallmann & Jaki (2017); emp.bayes for the empirical 
#' Bayes region described in Casella & Hwang (1983); hotelling for the Hotelling-type region
#' described in Wang et al (1999); limacon.asy for the limacon-shaped mimimum expected volume
#' region described in Brown et al (1995); limacon.fin for the finite-sample variant of the minimum
#' expected volume region described in Berger & Hsu (1996); standard.cor for the standard region
#' incorporating correlation between parameters described in Chew (1966); standard.ind for the
#' standard region ignoring correlation between parameters; tost for the two one-sided test (TOST)
#' intervals described in Schuirmann (1987); tseng for the mimimum expected interval length region
#' described in Tseng (2002); tseng.brown for the pseudo-empirical Bayes region described in Tseng
#' & Brown (1997). Available methods for confidence intervals are: expanded for the two one-sided 
#' est (TOST) procedure (Schuirmann 1987) using the expanded intervals described e.g., in Bofinger (1992) and Hsu
#' et al. (1994); fix.seq for the fixed sequence intervals described in Maurer et al (1995) and Hsu &
#' Berger (1999); tost for the two one-sided test (TOST) intervals described in Schuirmann (1987).
#' See also an overview and comparison of all methods in Pallmann & Jaki (2017).
ComputeConfidenceRegions <- function(mu_x,
                                     mu_y,
                                     n, 
                                     rho, 
                                     sigma_x, 
                                     sigma_y,
                                     num_points,
                                     alpha_value,
                                     methods,
                                     sample_num,
                                     sample){
  results <- list()
  
  #### CR Original data ---------------------------------------------------------------------------
  #
  # Confidence region based on original data & proposed reparametrisation
  confRegion_points <- confRegionBoundary(num_points = num_points, 
                                          alpha_value = alpha_value,
                                          sigma_x = sigma_x, 
                                          sigma_y = sigma_y,
                                          mu_x = mu_x, 
                                          mu_y = mu_y,
                                          rho = rho,
                                          n = n,
                                          label = 'original data')
  #### MLE estimators -----------------------------------------------------------------------------
  #
  # Get MLE estimators
  MLE <- mvnormMLE(sample = sample)
  
  #### Metrics Original data ----------------------------------------------------------------------
  #
  # Generate convex hull for CR generated from original data
  original_convhull <- convhulln(confRegion_points$p, options = 'FA')
  
  # Check if MLE is within the CR generated from original data
  original_is_inhull <- inhulln(original_convhull, matrix(c(MLE$mu_x_MLE, MLE$mu_y_MLE), 1, 2))
  
  # Gather the results 
  results <- append(results, list(list(
    'computation time' = confRegion_points$time,
    'volume' = original_convhull$vol, 
    'intersection volume' = original_convhull$vol,
    'theoretical volume' = original_convhull$vol, 
    'method' = confRegion_points$label,
    'is_inhull' = original_is_inhull,
    'convex hull' = original_convhull,
    'sample number' = sample_num,
    'sample size' = n
  )))
  
  #### CR from MLE --------------------------------------------------------------------------------
  #
  # Compute the CR from MLE estimators
  MLEconfRegion_points <- confRegionBoundary(num_points = num_points, 
                                             alpha_value = alpha_value,
                                             sigma_x = MLE$sigma_x_MLE, 
                                             sigma_y = MLE$sigma_y_MLE,
                                             mu_x = MLE$mu_x_MLE, 
                                             mu_y = MLE$mu_y_MLE,
                                             rho = MLE$rho_MLE,
                                             n = n,
                                             label = 'MLE')
  
  # Generate convex hull & intersection with CR from original data
  MLE_intersection <- intersectn(MLEconfRegion_points$p, confRegion_points$p)
  
  # Check if original mu_x & mu_y are within the convex hull of CR
  MLE_is_inhull <- inhulln(MLE_intersection$ch1, matrix(c(mu_x, mu_y), 1, 2))
  
  # Gather the results
  results <- append(results, list(list(
    'computation time' = MLEconfRegion_points$time,
    'volume' = MLE_intersection$ch1$vol, 
    'intersection volume' = MLE_intersection$ch$vol,
    'theoretical volume' = MLE_intersection$ch2$vol, 
    'method' = MLEconfRegion_points$label,
    'is_inhull' = MLE_is_inhull,
    'convex hull' = MLE_intersection$ch1,
    'sample number' = sample_num,
    'sample size' = n
  )))
  
  #### CR & SCI Other methods ---------------------------------------------------------------------
  #
  # Compute confidence regions for other methods from jocre package
  other_methods <- lapply(methods, function(method){
    metrics(sample = sample, 
            method = method, 
            mu_x = mu_x,
            mu_y = mu_y,
            alpha_value = alpha_value, 
            optimal_region = confRegion_points$p,
            sample_num = sample_num)
  })
  
  # gather all the results
  results <- append(results, unlist(other_methods, recursive = FALSE))
  return(results)
}






















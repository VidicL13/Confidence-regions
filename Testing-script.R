# Setup the environment
source('./libraries.r')
source('./functions.r')




# Make a data.frame of all possible combinations
values <- expand.grid(mu_x = c(0,10),
                      mu_y = c(0,10),
                      n = c(10,20,50), 
                      rho = c(0,0.5,0.9), 
                      sigma_x = c(sqrt(10),sqrt(100)), 
                      sigma_y = c(sqrt(10),sqrt(100)),
                      num_points = 180,
                      alpha_value = c(0.01, 0.05, 0.1),
                      num_runs = 100)

# select methods for cset function
methods <- list('boot.kern', 'emp.bayes')
             #    , 'hotelling',
             # 'limacon.asy', 'limacon.fin', 'standard.cor',
             # 'standard.ind', 'tost', 'tseng',
             # 'tseng.brown')

od <- 241
do <- 280

pb <- progress_bar$new(total = nrow(values[od:do,]),
                       width = 100,
                       clear = FALSE,
                       format = paste("Progress: [:bar] :percent | estimated time: :eta | elapsed time: :elapsed", '(', od, '-', do, ')'))

results <- apply(values[od:do,], 1, function(data_row){
  lapply(1:(data_row['num_runs']), function(i){
    # Generate bivariate normal sample
    sample <- mvnormSample(mu_x = data_row['mu_x'], 
                           mu_y = data_row['mu_y'], 
                           n = data_row['n'], 
                           rho = data_row['rho'], 
                           sigma_x = data_row['sigma_x'], 
                           sigma_y = data_row['sigma_y'])
    
    # Compute MLEs
    sample_mle <- mvnormMLE(sample = sample)
    
    #### Get CR-s and all the metrics
    results <- ComputeConfidenceRegions(mu_x = data_row['mu_x'],
                                        mu_y = data_row['mu_y'],
                                        n = data_row['n'], 
                                        rho = data_row['rho'], 
                                        sigma_x = data_row['sigma_x'], 
                                        sigma_y = data_row['sigma_y'],
                                        num_points = data_row['num_points'],
                                        alpha_value = data_row['alpha_value'],
                                        methods = methods,
                                        sample_num = i,
                                        sample = sample)
    pb$tick(1/data_row['num_runs'])
    # prepare the results
    lapply(results, function(item){
      data.frame(
        'method' = item[['method']],
        'sample number' = item[['sample number']],
        'sample size' = item[['sample size']],
        'computation time' = item[['computation time']],
        'volume' = item[['volume']], 
        'intersection volume' = item[['intersection volume']],
        'theoretical volume' = item[['theoretical volume']], 
        'is_inhull' = item[['is_inhull']],
        'mu_x' = data_row['mu_x'],
        'mu_y' = data_row['mu_y'],
        'rho' = data_row['rho'], 
        'sigma_x' = data_row['sigma_x'], 
        'sigma_y' = data_row['sigma_y'],
        'alpha_value' = data_row['alpha_value'],
        'mu_x_MLE' = sample_mle$mu_x_MLE,
        'mu_y_MLE' = sample_mle$mu_y_MLE,
        'rho_MLE' = sample_mle$rho_MLE, 
        'sigma_x_MLE' = sample_mle$sigma_x_MLE, 
        'sigma_y_MLE' = sample_mle$sigma_y_MLE
      )
    }) %>% bind_rows()
  }) %>% bind_rows()
}) %>% bind_rows()
# took 2 hours
write.csv2(results, file = paste0('./Data/data_',od ,'-', do, '.csv'), fileEncoding = 'UTF-8')



# took 2 hours
write.csv2(results, file = './Data/data_81-100.csv', fileEncoding = 'UTF-8')

# took 1 hour
write.csv2(results, file = './Data/data_1-10.csv', fileEncoding = 'UTF-8')

# took 5 hours
write.csv2(results, file = './Data/data_11-70.csv', fileEncoding = 'UTF-8')


# took 47 minutes
write.csv2(results, file = './Data/data_71-80.csv', fileEncoding = 'UTF-8')



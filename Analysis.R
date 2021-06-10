# load libraries
source('./libraries.r')
source('./functions.r')

# gather the data
CR_data <- lapply(list.files('./Data'), 
                  function(x) {read.csv2(file = paste0('./Data/', x), 
                                         fileEncoding = 'UTF-8', 
                                         row.names = 'X')}) %>% bind_rows()
CR_data <- results

# Summary of my data
CR_data %>% summary()

# confidence level of model based on sample size
CR_data %>% 
  group_by(method, alpha_value, sample.size, sample.number) %>% 
  summarise(mean_cov = mean(volume),
            mean_time = mean(computation.time),
            confidence_level = mean(as.numeric(is_inhull))) %>% 
  arrange(desc(mean_cov)) %>% 
  ggplot() + 
  geom_jitter(aes(x = method, y = confidence_level, col = method)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(sample.size ~ alpha_value)

# Average coverage based on alpha and sample size
CR_data %>% 
  group_by(method, alpha_value, sample.size, sample.number) %>% 
  summarise(mean_cov = mean(volume),
            mean_inter = mean(intersection.volume),
            confidence_level = mean(as.numeric(is_inhull))) %>% 
  ggplot() + 
  geom_boxplot(aes(x = method, y = mean_cov, fill = method)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(sample.size ~ alpha_value, scales = 'free')

# Total computation time
CR_data %>% 
  filter(!grepl('.*SCI', method)) %>% 
  group_by(method) %>% 
  summarise(comp_time = sum(computation.time)/(60*60), cases = n())

#### NONINFORMATIVE ----
# compare the confidence level based on size and mean value
CR_data %>% 
  filter(alpha_value == 0.1) %>% 
  group_by(method, alpha_value, sample.size, mu_x, mu_y) %>% 
  summarise(mean_cov = mean(volume),
            mean_inter = mean(intersection.volume),
            confidence_level = mean(as.numeric(is_inhull)),
            num_samples = n()) %>% 
  mutate(mux_muy = paste0(mu_x, '_', mu_y)) %>% 
  ggplot() +
  geom_point(aes(x = mux_muy, y = confidence_level, col = method, shape = mux_muy, size = mean_cov)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(method ~ sample.size)

#### NONINFORMATIVE ----
# compare the confidence level based on size and variance
CR_data %>% 
  filter(alpha_value == 0.1) %>%
  group_by(method, alpha_value, sample.size, sigma_x, sigma_y) %>% 
  summarise(mean_cov = mean(volume),
            mean_inter = mean(intersection.volume),
            confidence_level = mean(as.numeric(is_inhull)),
            num_samples = n()) %>% 
  mutate(sigmas = paste0(round(sigma_x), '_', round(sigma_y))) %>% 
  ggplot() +
  geom_point(aes(x = sigmas, y = confidence_level, col = method, shape = sigmas, size = mean_cov)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(method ~ sample.size)

#### NONINFORMATIVE ----
# compare the confidence level based on size and rho
CR_data %>% 
  group_by(method, alpha_value, sample.size, rho) %>% 
  summarise(confidence_level = mean(as.numeric(is_inhull))) %>% 
  ggplot() +
  geom_point(aes(x = method, y = confidence_level, col = alpha_value %>% as.character())) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap('rho')
  # facet_grid(sample.size ~ alpha_value)

# potential metric to represent quality of CR
CR_data %>% 
  group_by(method, sample.size) %>% 
  summarise(perc_coverage = mean(intersection.volume) / mean(theoretical.volume),
            confidence_level = mean(as.numeric(is_inhull)),
            vol = mean(volume)) %>%
  mutate(metrica = confidence_level / vol) %>% 
  arrange(desc(metrica)) %>% view()

CR_data %>% 
  mutate(mux_muy = paste0(mu_x, '_', mu_y)) %>% 
  group_by(is_inhull, mux_muy) %>% 
  summarise(MLE_x = mean(mu_x_MLE),
            MLE_y = mean(mu_y_MLE),
            n = n(),
            n_mle = sum(method == 'original data'))

# computed confidence level with actual one
CR_data %>% 
  group_by(method, alpha_value, sample.size) %>% 
  summarise(confidence_level = mean(as.numeric(is_inhull))) %>% 
  ggplot() +
  geom_point(aes(x = method, y = confidence_level, col = method), size = 5) + 
  geom_hline(aes(yintercept = 1-alpha_value)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(sample.size ~ alpha_value)
  
  
  
  
  









# load libraries
source('./libraries.r')

# gather the data
CR_data <- lapply(list.files('./Data'), 
                  function(x) {read.csv2(file = paste0('./Data/', x), 
                                         fileEncoding = 'UTF-8', 
                                         row.names = 'X')}) %>% bind_rows()

# confidence level of model based on sample size
CR_data %>% 
  group_by(method, alpha_value, sample.size, sample.number) %>% 
  summarise(mean_cov = mean(volume),
            mean_time = mean(computation.time),
            confidence_level = mean(as.numeric(is_inhull))) %>% 
  arrange(desc(mean_cov)) %>% 
  filter(alpha_value == 0.01) %>% 
  ggplot() + 
  geom_boxplot(aes(x = method, y = confidence_level, fill = method)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap('sample.size')

CR_data %>% left_join(values, by = c("mu_x", "mu_y", 'sample.size' = "n", "rho", "sigma_x", "sigma_y", "alpha_value"))

names(values)

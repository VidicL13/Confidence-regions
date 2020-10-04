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
  ggplot() + 
  geom_boxplot(aes(x = method, y = confidence_level, fill = method)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(sample.size ~ alpha_value)

CR_data %>% 
  group_by(method, alpha_value, sample.size, sample.number) %>% 
  summarise(mean_cov = mean(volume),
            mean_inter = mean(intersection.volume),
            confidence_level = mean(as.numeric(is_inhull))) %>% 
  ggplot() + 
  geom_boxplot(aes(x = method, y = mean_inter, fill = method)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(sample.size ~ alpha_value, scales = 'free')

# Total computation time
CR_data %>% 
  filter(!grepl('.*SCI', method)) %>% 
  group_by(method) %>% 
  summarise(comp_time = sum(computation.time)/(60*60))

args = commandArgs(trailingOnly = TRUE)
distro_dir = args[1]
healthy_libs_str = args[2]
healthy_med_file = args[3]

library(tidyverse)

healthy_libs_distros = paste0(distro_dir, "/", unlist(strsplit(healthy_libs_str, " ")), "_gc_distro.csv")

read_in_gc = function(gc_csv){
  read.csv(gc_csv, header = T)
}

healthy_list = lapply(healthy_libs_distros, read_in_gc)

# Bind
healthy_all = do.call(rbind, healthy_list)

# Summarize
healthy_med =
  healthy_all %>%
  group_by(gc_strata) %>%
  summarise(med_frag_fract = median(fract_frags))

# Save
# saveRDS(healthy_med, file = healthy_med_file)
healthy_med %>%
  write.csv(file = healthy_med_file, row.names = F)

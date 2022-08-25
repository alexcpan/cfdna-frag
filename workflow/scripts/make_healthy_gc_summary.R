#########1#########2#########3#########4#########5#########6#########7#########8
source(snakemake@config[["r_lib_loads"]])

# Read in healthy plasma gc distributions
all_distros = list.files(path = paste0(snakemake@config[["data_dir"]],"/frag"),
                       pattern = "gc_distro")
healthy_libs = snakemake@config[["healthy_plasma"]]

saveRDS(all_distros, file = snakemake@output[[1]])
healthy_distros = paste0(snakemake@config[["data_dir"]],"/frag/",
                         grep(paste(healthy_libs, collapse="|"),
                              all_distros, value = T))


read_in_gc = function(gc_csv){
  read.csv(gc_csv, header = T)
}
healthy_list = lapply(healthy_distros, read_in_gc)

# Bind
healthy_all = do.call(rbind, healthy_list)

# Summarize
healthy_med =
  healthy_all %>%
  group_by(gc_strata) %>%
  summarise(med_frag_fract = median(fract_frags))

# Save
saveRDS(healthy_med, file = snakemake@output[["healthy_med"]])

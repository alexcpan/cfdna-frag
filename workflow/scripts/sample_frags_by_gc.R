#########1#########2#########3#########4#########5#########6#########7#########8
source(snakemake@config[["r_lib_loads"]])

healthy_fract = readRDS(snakemake@input[["healthy_med"]])
frag_file = read.table(snakemake@input[["frag_bed"]], sep = '\t', header = F)

reject_sample = function(frag_bed,healthy_fract){
  names(frag_bed) = c("chr", "start", "end", "gc_raw", "len")
  sampled = frag_bed %>%
    mutate(gc_strata = round(gc_raw, 2)) %>%
    left_join(healthy_fract, by = "gc_strata") %>%
    mutate(include = ifelse(runif(nrow(.),0,1) < med_frag_fract / max(med_frag_fract, na.rm = T), "yes", "no")) %>%
    filter(include == "yes")
  return(sampled)
}

sampled = reject_sample(frag_file, healthy_fract)

write.table(sampled, sep = "\t", col.names = F, row.names = F, quote = F, file = snakemake@output[[1]])

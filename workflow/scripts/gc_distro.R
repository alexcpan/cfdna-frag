#########1#########2#########3#########4#########5#########6#########7#########8

# Source config
source(snakemake@params[[1]])

# Read in modified bed
bed = read.table(snakemake@input[[1]], sep = '\t')
names(bed) = c("chr","start","end","gc_raw","len")

# Generate distribution csv
distro =
  bed %>%
  # Round GC
  mutate(gc_strata = round(gc_raw, 2)) %>%
  # Count frags per strata
  count(gc_strata) %>%
  # Get fraction frags
  mutate(fract_frags = n/sum(n)) %>% mutate(library_id = gsub("_frag.bed", "", gsub("^.*lib", "lib", snakemake@input[[1]]))) %>%
  select(library_id,gc_strata,fract_frags) %>%
  write.csv(file = snakemake@output[[1]], row.names = F)

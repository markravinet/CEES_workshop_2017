## combine SNP statistics ##
rm(list = ls())

library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)

# FREQ
# find files
freq_files <- list.files(path = "./data/", pattern = "*.frq", full.names = T)
# read in data
freq <- lapply(freq_files, read_delim, delim = "\t", skip = 1, 
               col_names = c("chr", "pos", "n_all", "n_chr", "freq1", "freq2"))

# calculate minor allele freq
freq <- lapply(freq, function(x){
  x$maf <- apply(as.matrix(select(x, freq1, freq2)), 1, function(y) min(y))
  x <- mutate(x, id = paste0(chr, "_", pos))
  select(x, -freq1, -freq2)
})

# join together
freq_join <- left_join(freq[[1]], freq[[2]], by = "id", suffix = c("_angsd", "_fb"))
freq_join <- left_join(freq_join, freq[[3]], by = "id")
freq_join <- rename(freq_join, chr_gatk = chr, pos_gatk = pos, n_all_gatk = n_all, n_chr_gatk = n_chr,
                    maf_gatk = maf)
# now add gl mafs
gl_mafs <- read_delim("./data/angsd_gl.mafs", delim = "\t", col_names = c("chr", "pos", "maf_angsd_gl"))
gl_mafs <- mutate(gl_mafs, id = paste0(chr, "_", pos))
freq_join <- select(left_join(freq_join, gl_mafs), -chr, -pos)

# clean up
freq <- select(freq_join, chr = chr_angsd, pos = pos_angsd, contains("n_chr"), contains("maf"), id)
freq <- na.omit(freq)

# DEPTH
# find files
depth_files <- list.files(path = "./data/", pattern = "*.ldepth.mean", full.names = T)
# read in data
depth <- lapply(depth_files, read_delim, delim = "\t", skip = 1, 
               col_names = c("chr", "pos", "mean_depth", "var_depth"))

# set id
depth <- lapply(depth, function(x){
  mutate(x, id = paste0(chr, "_", pos))
})

# join together
depth_join <- left_join(depth[[1]], depth[[2]], by = "id", suffix = c("_angsd", "_fb"))
depth_join <- left_join(depth_join, depth[[3]], by = "id")
# clean up
depth_join <- rename(depth_join, chr_gatk = chr, pos_gatk = pos, mean_depth_gatk = mean_depth)

# clean up
depth <- select(depth_join, chr = chr_angsd, pos = pos_angsd, contains("mean_depth"), -mean_depth_angsd, id)
depth <- na.omit(depth)

# MISSING
# find files
miss_files <- list.files(path = "./data/", pattern = "*.lmiss", full.names = T)
# read in data
miss <- lapply(miss_files, read_delim, delim = "\t", skip = 1, 
                col_names = c("chr", "pos", "n_data", "n_geno_fil", "n_miss", "f_miss"))

# set id
miss <- lapply(miss, function(x){
  mutate(x, id = paste0(chr, "_", pos))
})

# join together
miss_join <- left_join(miss[[1]], miss[[2]], by = "id", suffix = c("_angsd", "_fb"))
miss_join <- left_join(miss_join, miss[[3]], by = "id")
# clean up
miss_join <- rename(miss_join, f_miss_gatk = f_miss)

# clean up
miss <- select(miss_join, chr = chr_angsd, pos = pos_angsd, contains("f_miss"), id)
miss <- na.omit(miss)

# QUAL
# find files
qual_files <- list.files(path = "./data/", pattern = "*.lqual", full.names = T)
# read in data
qual <- lapply(qual_files, read_delim, delim = "\t", skip = 1, 
               col_names = c("chr", "pos", "qual"))

# set id
qual <- lapply(qual, function(x){
  mutate(x, id = paste0(chr, "_", pos))
})

# join together
qual_join <- left_join(qual[[1]], qual[[2]], by = "id", suffix = c("_angsd", "_fb"))
qual_join <- left_join(qual_join, qual[[3]], by = "id")
# clean up
qual_join <- rename(qual_join, qual_gatk = qual)

# clean up
qual <- select(qual_join, chr = chr_angsd, pos = pos_angsd, contains("qual"), id)
qual <- na.omit(qual)

snp_data <- left_join(freq, miss, by = "id")
snp_data <- left_join(snp_data, depth, by = "id")
snp_data <- left_join(snp_data, qual, by = "id")
#  now clean up
snp_data <- select(snp_data, chr = chr.x, pos = pos.x, contains("maf"), contains("f_miss"),
       contains("mean_depth"), contains("qual"))

# combine data 
# snp_data <- tbl_df(cbind(select(freq, chr, pos, contains("maf")),
#                          select(miss, contains("f_miss")),
#                          select(depth, contains("mean_depth")),
#                          select(qual, contains("qual"))))

write_csv(snp_data, "./snp_stats_combined.csv")

### INDIVIDUAL DATA
# MISS
# find files
miss_files <- list.files(path = "./data/", pattern = "*.imiss", full.names = T)
methods <- c("angsd", "fb", "gatk")
# read in data
miss <- lapply(miss_files, read_delim, delim = "\t", skip = 1, 
               col_names = c("ind", "n_data", "n_geno_fil", "n_miss", "f_miss"))
names(miss) <- methods
ind <- select(miss[[3]], ind)
# combine them all
f_miss <- do.call(cbind, lapply(miss, select, f_miss))
colnames(f_miss) <- paste0(colnames(f_miss), "_", methods)

# DEPTH
# find files
depth_files <- list.files(path = "./data/", pattern = "*.idepth", full.names = T)
# read in data
depth <- lapply(depth_files, read_delim, delim = "\t", skip = 1,
               col_names = c("ind", "n_sites", "mean_depth"))

depth <- do.call(cbind, lapply(depth, select, mean_depth))
colnames(depth) <- paste0(colnames(depth), "_", methods)

# HET
# find files
het_files <- list.files(path = "./data/", pattern = "*.het", full.names = T)
# read in data
het <- lapply(het_files, read_delim, delim = "\t", skip = 1,
                col_names = c("ind", "ho", "he", "nsites", "f"))

het <- do.call(cbind, lapply(het, select, f))
colnames(het) <- paste0(colnames(het), "_", methods)

# create a final file
ind_data <- tbl_df(cbind(ind, f_miss, depth, het))
write_csv(ind_data, "./ind_stats_combined.csv")


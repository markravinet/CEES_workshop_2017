rm(list = ls())
library(tidyverse)
# or #
library(ggplot2)
library(readr)
library(dplyr)

# read in data
snps <- read_csv("./snp_stats_combined.csv", col_names = T)
# downsample to just 10%
snps_d <- sample_frac(snps, 0.1)

# compare the ANGSD mafs (hard calls vs genotype likelihoods)
a <- ggplot(snps_d, aes(maf_angsd, maf_angsd_gl))
a + geom_point()

# compare GATK and FreeBayes
a <- ggplot(snps_d, aes(maf_gatk, maf_fb,
                        colour = log(mean_depth_gatk)))
a <- a + geom_point(position = "jitter", alpha = 0.6)
a + geom_abline(lty = 2, colour = "red")
a + geom_smooth(lty = 2, colour = "red")

# GATK vs ANGSD GL estimated MAFs
a <- ggplot(snps_d, aes(maf_gatk, maf_angsd_gl))
a <- a + geom_point(position = "jitter", alpha = 0.1)
a + geom_abline(lty = 2, colour = "red")

# FreeBayes vs ANGSD GL estimated MAFs
a <- ggplot(snps_d, aes(maf_fb, maf_angsd_gl))
a <- a + geom_point(position = "jitter", alpha = 0.1)
a + geom_abline(lty = 2, colour = "red")

# mean depth of FreeBayes vs GATK
a <- ggplot(snps_d, aes(mean_depth_fb, mean_depth_gatk))
a <- a + geom_point(position = "jitter", alpha = 0.3)
a + geom_abline(lty = 2, colour = "red")

# plot missing data
f_miss <- select(snps_d, contains("f_miss"))
f_miss <- gather(f_miss, key = "caller", value = "f_miss")
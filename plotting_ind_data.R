rm(list = ls())

library(tidyverse)

# read in inds data
inds <- read_csv("ind_stats_combined.csv", col_names = T)

# examine the distribution of missing data
f_miss <- gather(select(inds, contains("f_miss")), key = "method", value = "f_miss")
ggplot(f_miss, aes(f_miss)) + geom_histogram() + facet_wrap(~method) 

# plot het vs missing data
ggplot(inds, aes(mean_depth_gatk, f_gatk)) + geom_point()
ggplot(inds, aes(mean_depth_fb, f_fb)) + geom_point()
# this lets you see the relationship between depth and also heterozygosity (measured by FIS)
# here you can see that for FB there is quite a clear positive relatioship ut that this isn't the case
# for GATK
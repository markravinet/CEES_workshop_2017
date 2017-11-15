## Plotting SNP values to set thresholds for filtering
rm(list = ls())

library(ggplot2)
library(dplyr)

# read in data
snps <- read.table("./snp_data.txt", header = F, na.strings = ("."))
snps <- tbl_df(snps)

# rename columns
colnames(snps) <- c("an", "ac", "af", "qual", "qd", "dp", "fs",
                    "mq", "bsqr", "readPosRS") 

# this is no longer necessary as we fixed this with the input
# remove NAs in the mq column
# tidyverse method
# snps <- mutate(snps, mq = as.numeric(as.character(snps$mq)))
# baseR (you choose!)
# snps$mq <- as.numeric(as.character(snps$mq))

# plot some data values to get an idea of what we might want to filter on
# allele number
ggplot(snps, aes(an)) + geom_density(fill = "dodgerblue1", alpha = 0.3) + theme_light()

# number of alleles called (max is 70 as there are 35 inds)

# allele count
ggplot(snps, aes(ac)) + geom_density(fill = "dodgerblue1", alpha = 0.3) + theme_light()
# number of alleles in genotypes for the ALT allele

# allele frequency for ALT allele
ggplot(snps, aes(af)) + geom_density(fill = "dodgerblue1", alpha = 0.3) + theme_light()
# number of alleles in genotypes for the ALT allele

# quality
ggplot(snps, aes(qual)) + geom_density(fill = "dodgerblue1", alpha = 0.3) + theme_light() + xlim(0, 10000)
# quality scores

# quality by depth
ggplot(snps, aes(qd)) + geom_density(fill = "dodgerblue1", alpha = 0.3) + theme_light()

# fisher strand bias
ggplot(snps, aes(fs)) + geom_density(fill = "dodgerblue1", alpha = 0.3) + theme_light() + xlim(0, 30)

# mq
ggplot(snps, aes(mq)) + geom_density(fill = "dodgerblue1", alpha = 0.3) + theme_light() + xlim(0, 100)

# bsqr
ggplot(snps, aes(bsqr)) + geom_density(fill = "dodgerblue1", alpha = 0.3) + theme_light()

# readPos
ggplot(snps, aes(readPosRS)) + geom_density(fill = "dodgerblue1", alpha = 0.3) + theme_light()

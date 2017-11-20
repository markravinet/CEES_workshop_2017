rm(list = ls())
library(tidyverse)

# read in point data
fst <- read_delim("./house_spanish.weir.fst", delim = '\t')
fst <- rename(fst, chr = CHROM, bp = POS, fst = WEIR_AND_COCKERHAM_FST)
# plot using ggplot
a <- ggplot(fst, aes(bp/10^6, fst)) 
a <- a + geom_point(alpha = 0.1)
a <- a + xlab("Position (Mb)") + ylab(expression(italic("F")[ST]))
a + theme_light()

# read in windowed data
fst <- read_delim("./house_spanish.windowed.weir.fst", delim = '\t')
fst <- rename(fst, chr = CHROM, start = BIN_START, 
              end = BIN_END, nvar = N_VARIANTS,
              weighted_fst = WEIGHTED_FST,
              mean_fst = MEAN_FST)
# generate midpoint for plotting
fst <- mutate(fst, mid = (start + end-1)/2)
# draw plot
a <- ggplot(fst, aes(mid/10^6, weighted_fst)) 
a <- a + geom_line()
a <- a + xlab("Position (Mb)") + ylab(expression(italic("F")[ST]))
a + theme_light()




# CEES Workshop 2017 - script guidelines

**You will need to alter some of the paths in these script to get them to run for you. However you should also be able to look at the scripts and get an idea how to make them work for your own data**

All of these scripts will work with the toy sparrow dataset which you can find here:

```
/projects/cees2/in_progress/cees_workshop
```
This is also your the `$DATA_DIR` variable, so if you run `source prep_environment.sh` you can access it that way.

## Getting access to the reads 

Raw reads are available here:

```
/projects/cees2/in_progress/cees_workshop/full_raw
```

Trimmed reads are here:

```
/projects/cees2/in_progress/cees_workshop/full_trimmed
```

The vcfs are present here:

```
/projects/cees2/in_progress/cees_workshop/full_trimmed
```

## Mapping reads

There are two scripts for mapping the reads. One which we wrote in class which you can run on a just a few samples - that is `map_reads_class.sh`. 

There is also another script which I wrote to map the data we used downstream, this is a bit more complicated but you should be able to use it as a template `map_toy_reads.sh` - this includes information on setting readgroups. 

## GATK pipeline

To run the `GATK` pipeline, you will need to run a series of scripts. 

1. First of all you need to realign the bams using `bam_realignment.sh`
2. After this, we generate gVCFs using `haplotype_caller.sh`
3. Finally you need to run the joint genotyping module across all of the gVCFs with `joint_genotyping.sh`
4. To annotate the VCF with filters, you need to use `annotate_gatk_filters.sh

## VCF filtering and statistics

To apply the `GATK` filters and also filter vcfs in general, use the `filtering_and_pruning_variants.sh` script. This will also output the files you need to run some of the R scripts.

We also filtered the complete chr15 vcf with `randomly_sample_vcf.sh`. Statistics for chr15 were generated using `chr15_stats.sh`.

Finally we calculated *F*<sub>ST</sub> values using the script `calculate_fst_with_vcftools.sh`.

## R scripts for plotting

To plot the thresholds from the randomly sampled vcf data, you can use the `plotting_snps_query_data.R` script. I also generated [a report of this data](https://markravinet.github.io/sparrow_SNP_thresholds.html) using [Rmarkdown](http://rmarkdown.rstudio.com/). To see the markdown file I used, look at `sparrow_SNP_thresholds.Rmd`. 

I combined the data from chr15 using `combine_datasets.R`- this script is clunky and not the most straightforward to understand, but it should give you some idea of how I put together such a large dataset. 

As you would expect, `plot_fst.R` plots the *F*<sub>ST</sub> values we calculated using `vcftools`.

One the last day we examined the difference between the output of the variant callers/pipelines. You can recreate this by using the `plotting_diffs_btwn_callers.R` script.

One thing we didn't get do (meaning I forgot) was to look at the differences between individuals. I wrote a quick Rscript to do this for you - so feel free to have a look at it - `plotting_ind_data.R`

## Datasets on Github

Lastly, there are couple of datasets present on the Github that can be used with the R scripts.

These are:

1.`ind_stats_combined.csv` - individual statistics calculated from chr15
2. `snp_stats_combined.csv` - as above but for each site on chr15 called using all methods
3. `snp_data.txt` - the randomly sampled data used to set `GATK` filters.



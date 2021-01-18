# lTWAS

## Tutorial

You may run the following command:

```
python3 ltwas.py location_of_sumstats1 location_of_sumstats2 \
--N1 sample_size_1 \
--N2 sample_size_2 \
--bfile location_of_reference_panel \
--chr chromosome_of_the_region \
--start start_genomic_position_of_the_region\
--end end_genomic_position_of_the_region\
--h1 local_heritability_of_trait1\
--h2 local_heritability_of_trait2\
--out location_of_results
```
### Explanation of Command-Line Arguments

- The first two arguments denote the locations of the first and second summary statistics files. These files may be compressed using gzip, bz2, zip, xz, or not compressed at all. The program will infer the compression method if the files end with .gz, .bz2, .zip, xz, respectively. As previously mentioned, we assume that the files are in the standard format that `ldsc` understands.

- The `N1` and `N2` arguments (optional) denote the sample sizes of the summary statistics files. If they are not provided, they will be inferred from the summary statistics files.

- The `bfile` argument denotes the prefix of the `.bed/.bim/.fam` genotypic data file. Please provide only one file which contain the genomic region you are interested in.

- The `chr` argument denotes the chromosome of the region.

- The `start` and `end` denote the start and the end genomic position of the region, respectively.

- The `h1` and `h2` denote the local heritability of trait1 and trait2 for the region, respectively.	

- The `out` flag denotes the file location for the results to be outputted to.

### Additional Command-Line Arguments

- `--thread`: The thread number used for calculation (optional). If they are not provided, the default is the number of detected CPUs.

### Explanation of Output
The output will be a whitespace-delimited text file, with the rows corresponding to different annotations and the columns as such:

- `rho`: The estimation of local genetic covariance.
- `corr`: The estimation of local genetic correlation.
- `h2_1`: The estimation of local heritability of the first trait by method of moments.
- `h2_2`: The estimation of local heritability of the second trait by method of moments.
- `var`: The variance of the estimation of local genetic covariance.
- `p`: The p value of local genetic covariance.
- `m`: The number of SNPs involved in the estimation of local genetic covariance in the genomic region.
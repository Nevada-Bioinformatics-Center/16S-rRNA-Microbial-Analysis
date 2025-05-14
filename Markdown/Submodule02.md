# Submodule 2: Data Preprocessing and Quality Control


## Overview
In this module, we’ll cover essential steps for preprocessing and quality control of 16S rRNA sequencing data. Effective data preprocessing ensures that our analyses are based on high-quality data, minimizing errors and maximizing accuracy. This module will guide you through importing raw sequencing data, performing quality filtering, and removing potential errors such as chimeric sequences.

We’ll explain the **FastQ** file formats and use **DADA2** for quality assessment and feature extraction. This hands-on practice will import and process real data from the Nevada Wolfpack Study which we will use in the next module for downstream analyses.

<center>
    <img src="../images/2.jpg" alt="2"/>
</center>

## Learning Objectives
+ Preprocess raw 16S rRNA sequencing data effectively.
+ Perform quality control checks on sequencing data.

## Import Sequencing Data

### Demultiplexing

Our starting point is demultiplexed samples sorted in individual files. Demultiplexing is a crucial step in processing high-throughput sequencing data when multiple samples are sequenced together in a single run (a process called multiplexing). Multiplexing allows researchers to save time and resources by sequencing multiple samples at once, but it introduces the need to sort, or "demultiplex," the data afterward.

In multiplexing, unique DNA "barcode" sequences are added to each sample before sequencing. During demultiplexing, these barcodes are used to identify which sequences belong to each original sample. The demultiplexing step occurs immediately after the sequencing machine outputs the raw data.

Demultiplexing software scans the data, recognizes the barcodes, and then separates sequences into individual files for each sample. Once demultiplexed, the data can be analyzed separately for each sample, moving on to further steps such as quality control, trimming, and chimera removal. Demultiplexing is typically performed by the sequencing platform's software and files of individual samples are given to the researcher. 

### Working with Our Samples

Working with microbiome data can be challenging. Traditional analysis methods cluster sequences based on similarity (operational taxonomic units (OTUs)), often missing the true biological diversity in the sample and sometimes lumping together different species. This can reduce the accuracy of the analysis, especially when studying complex microbial communities. Also, errors are often introduced during the sequencing process, such as substitutions (wrong bases) and chimera (sequences formed by two different DNA strands merging). These errors must be properly handled to avoid incorrect species identification. 

<center>
    <img src="../images/dada2.png" alt="dada2" width="600"/>
</center>

We will use [DADA2](https://benjjneb.github.io/dada2/tutorial_1_8.html) to load, process, and quality control our sequencing files. DADA2 is a widely used software package in bioinformatics for analyzing microbial communities from DNA sequencing data and handling errors. It works particularly well in processing 16S rRNA data. 

### How DADA2 Works

DADA2 takes a new approach by attempting to resolve individual sequences down to the level of single-nucleotide differences while correcting for sequencing errors<sup>1</sup>. This approach is called **Amplicon Sequence Variants (ASVs)**, which are highly accurate sequences that represent actual biological diversity in the sample without clustering.

The DADA2 workflow includes several key steps:

 1. **Inspecting Read Quality Profiles:** Before processing, DADA2 allows us to visualize the quality scores across all sequences. Quality scores provide an estimate of the accuracy of each base call in a sequence. By inspecting these scores, we can identify where the quality drops, which helps determine appropriate trimming parameters.
 2. **Filtering and Trimming:** Low-quality sequences and ambiguous bases can introduce errors in downstream analysis. In this step, we set quality thresholds to remove reads that don’t meet quality standards. Trimming removes low-quality ends of reads based on quality scores, helping to retain only high-quality sequences for further analysis.
 3. **Learning Error Rates:** DADA2 builds a model of the specific error rates associated with the sequencing run. Since error patterns vary between sequencing technologies and even between runs, DADA2 analyzes a subset of data to estimate these unique error rates. This error model helps DADA2 accurately correct errors in the sequencing reads.
 4. **Inferring Sample Composition (ASV Inference):** Using the error model, DADA2 corrects sequencing errors to identify unique, true biological sequences known as Amplicon Sequence Variants (ASVs). ASVs provide finer resolution than Operational Taxonomic Units (OTUs) and often represent exact biological sequences rather than clustered approximations.
 5. **Merging Paired Reads:** Since we are using paired-end sequencing, this step merges the forward and reverse reads that cover overlapping regions of the same DNA fragment. The merged reads typically provide higher accuracy by combining information from both ends, but reads that fail to overlap are discarded.
 6. **Constructing a Sequence Table:** After ASVs are identified, DADA2 constructs a sequence table, which is a matrix of samples and their ASV counts. This table serves as the primary data structure for downstream analysis, showing how many times each ASV appears in each sample.
 7. **Removing Chimeras:** Chimeric sequences are artifacts from PCR amplification, formed by combining segments from different DNA templates. DADA2 identifies these artificial sequences by comparing ASVs and removes them to reduce false positives in microbial identification.
 8. **Tracking Reads Through the Pipeline:** Throughout the pipeline, DADA2 tracks the number of reads retained at each step. This read-tracking helps to monitor data quality and determine where reads may have been lost, providing a quality control check on each processing stage.
 9. **Assigning Taxonomy:** After obtaining high-quality ASVs, DADA2 compares them to a reference database (e.g., SILVA, Greengenes) to assign taxonomic labels. This step allows us to identify the microbial taxa present in each sample down to species, genus, or other taxonomic levels, depending on the reference database used.

This process gives us highly accurate and reproducible results. We will start by installing and loading DADA2 and necessary packages.



``` r
#Install the DADA2 and other packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

if (!requireNamespace("dada2", quietly = TRUE)) {
    BiocManager::install("dada2")
}

if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

if (!requireNamespace("patchwork", quietly = TRUE)) {
    devtools::install("patchwork")
}
```


``` r
#Load the packages
library(dada2)
library(ggplot2)
library(patchwork)
```


Our sequencing files are stored in FASTQ format, a common file type for sequencing data. Each entry (or sequence) in a FASTQ file contains four lines: the sequence identifier, the DNA sequence, a separator, and a quality score line. This format is compact and efficient, making it essential for storing large-scale sequencing data.

<center>
    <img src="../images/fastq.png" alt="fastq" width="900"/>
</center>

FASTQ files follow the same format with each sequence entry:
 1. **Sequence Identifier:** Starts with an "@" symbol, followed by a unique identifier for the sequence. It may also contain additional details, such as the instrument used and the read number.
 2. **DNA Sequence:** The actual sequence of DNA bases (A, T, C, G) observed in the sample. The length of this line is equal to the number of bases sequenced.
 3. **Separator Line:** A "+" symbol acts as a separator between the sequence and its quality scores.
 4. **Quality Scores:** A string of characters representing the quality of each base in the sequence. Each character encodes a [Phred quality score](https://en.wikipedia.org/wiki/Phred_quality_score), where higher scores indicate more confidence in the accuracy of the base call. The quality string is the same length as the sequence line.

Now we will download our data and import it with DADA2.


``` r
# Download data
system("wget -r -np -R 'index.html*' --no-check-certificate --cut-dirs=2 -nH https://biox.unr.edu/ftp/16S_Training_Module/")
```


``` r
# Check files
path <- "fastqs"
list.files(path)
```

```
##  [1] "filtered"                       "WP-11O29_S66_L001_R1_001.fastq"
##  [3] "WP-11O29_S66_L001_R2_001.fastq" "WP-1641C_S68_L001_R1_001.fastq"
##  [5] "WP-1641C_S68_L001_R2_001.fastq" "WP-1J0A0_S50_L001_R1_001.fastq"
##  [7] "WP-1J0A0_S50_L001_R2_001.fastq" "WP-3GH4E_S51_L001_R1_001.fastq"
##  [9] "WP-3GH4E_S51_L001_R2_001.fastq" "WP-4SP0U_S57_L001_R1_001.fastq"
## [11] "WP-4SP0U_S57_L001_R2_001.fastq" "WP-6RCV3_S39_L001_R1_001.fastq"
## [13] "WP-6RCV3_S39_L001_R2_001.fastq" "WP-7FAEC_S47_L001_R1_001.fastq"
## [15] "WP-7FAEC_S47_L001_R2_001.fastq" "WP-8SXNH_S53_L001_R1_001.fastq"
## [17] "WP-8SXNH_S53_L001_R2_001.fastq" "WP-8V6E3_S72_L001_R1_001.fastq"
## [19] "WP-8V6E3_S72_L001_R2_001.fastq" "WP-95020_S58_L001_R1_001.fastq"
## [21] "WP-95020_S58_L001_R2_001.fastq" "WP-9ZJJI_S55_L001_R1_001.fastq"
## [23] "WP-9ZJJI_S55_L001_R2_001.fastq" "WP-ATJLX_S75_L001_R1_001.fastq"
## [25] "WP-ATJLX_S75_L001_R2_001.fastq" "WP-DFG11_S74_L001_R1_001.fastq"
## [27] "WP-DFG11_S74_L001_R2_001.fastq" "WP-EC5AP_S42_L001_R1_001.fastq"
## [29] "WP-EC5AP_S42_L001_R2_001.fastq" "WP-G375R_S73_L001_R1_001.fastq"
## [31] "WP-G375R_S73_L001_R2_001.fastq" "WP-G9EBL_S40_L001_R1_001.fastq"
## [33] "WP-G9EBL_S40_L001_R2_001.fastq" "WP-J06MK_S31_L001_R1_001.fastq"
## [35] "WP-J06MK_S31_L001_R2_001.fastq" "WP-J0B3X_S33_L001_R1_001.fastq"
## [37] "WP-J0B3X_S33_L001_R2_001.fastq" "WP-KM8BK_S62_L001_R1_001.fastq"
## [39] "WP-KM8BK_S62_L001_R2_001.fastq" "WP-KXGT6_S37_L001_R1_001.fastq"
## [41] "WP-KXGT6_S37_L001_R2_001.fastq" "WP-MRCZG_S36_L001_R1_001.fastq"
## [43] "WP-MRCZG_S36_L001_R2_001.fastq" "WP-QIGEC_S59_L001_R1_001.fastq"
## [45] "WP-QIGEC_S59_L001_R2_001.fastq" "WP-RAJES_S35_L001_R1_001.fastq"
## [47] "WP-RAJES_S35_L001_R2_001.fastq" "WP-RPV50_S49_L001_R1_001.fastq"
## [49] "WP-RPV50_S49_L001_R2_001.fastq" "WP-SUFNU_S34_L001_R1_001.fastq"
## [51] "WP-SUFNU_S34_L001_R2_001.fastq" "WP-TDY7G_S32_L001_R1_001.fastq"
## [53] "WP-TDY7G_S32_L001_R2_001.fastq" "WP-U4XI6_S44_L001_R1_001.fastq"
## [55] "WP-U4XI6_S44_L001_R2_001.fastq" "WP-UIXW3_S45_L001_R1_001.fastq"
## [57] "WP-UIXW3_S45_L001_R2_001.fastq" "WP-X576L_S56_L001_R1_001.fastq"
## [59] "WP-X576L_S56_L001_R2_001.fastq" "WP-X803O_S69_L001_R1_001.fastq"
## [61] "WP-X803O_S69_L001_R2_001.fastq" "WP-XDI0I_S65_L001_R1_001.fastq"
## [63] "WP-XDI0I_S65_L001_R2_001.fastq" "WP-YF6XT_S61_L001_R1_001.fastq"
## [65] "WP-YF6XT_S61_L001_R2_001.fastq" "WP-YRYYF_S67_L001_R1_001.fastq"
## [67] "WP-YRYYF_S67_L001_R2_001.fastq" "WP-YY8WG_S38_L001_R1_001.fastq"
## [69] "WP-YY8WG_S38_L001_R2_001.fastq" "WP-Z01LM_S43_L001_R1_001.fastq"
## [71] "WP-Z01LM_S43_L001_R2_001.fastq"
```



``` r
# Forward and reverse fastq filenames have format: SAMPLENAME_L001_R1_001.fastq and SAMPLENAME_L001_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_L001_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_L001_R2_001.fastq", full.names = TRUE))
```



``` r
# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```

Check to make sure our sample names are correct.


``` r
# Display sample names and number of samples
cat("Sample names:\n")
```

```
## Sample names:
```

``` r
cat(sample.names, sep = "\n")
```

```
## WP-11O29
## WP-1641C
## WP-1J0A0
## WP-3GH4E
## WP-4SP0U
## WP-6RCV3
## WP-7FAEC
## WP-8SXNH
## WP-8V6E3
## WP-95020
## WP-9ZJJI
## WP-ATJLX
## WP-DFG11
## WP-EC5AP
## WP-G375R
## WP-G9EBL
## WP-J06MK
## WP-J0B3X
## WP-KM8BK
## WP-KXGT6
## WP-MRCZG
## WP-QIGEC
## WP-RAJES
## WP-RPV50
## WP-SUFNU
## WP-TDY7G
## WP-U4XI6
## WP-UIXW3
## WP-X576L
## WP-X803O
## WP-XDI0I
## WP-YF6XT
## WP-YRYYF
## WP-YY8WG
## WP-Z01LM
```

``` r
cat("There are", length(sample.names), "samples")
```

```
## There are 35 samples
```

An .rds file saves R objects that can be quickly loaded back into R using `readRDS()`. While .rds files are optimized for R and can't be viewed directly outside of R, they are efficient for storing individual R objects and preserving their exact structure. This format is especially useful when you need to save and reload between sessions.


``` r
# Save our data
saveRDS(sample.names, file = 'sampleNames.rds')
```

## 2. Data Preprocessing and Quality Control

### Read Quality Score
We will first look at the quality of the individual bases in the reads by pulling the Phred score for each read in a sample. The heatmap produced shows the frequency of each quality score at each base over all the reads for a given sample. The green line represents the median with orange lines representing quartiles. 

*Note:* Some of these steps take several minutes to run. 


``` r
# View quality profiles of forward reads for the first 2 samples
plotQualityProfile(fnFs[1:2])
```

![](../Submodule02_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

*This step will take several minutes to run.*


``` r
# We can create a plot for all of the samples and
p1 = plotQualityProfile(fnFs) +
  ggtitle('Forward Reads')

p2 = plotQualityProfile(fnRs)+
  ggtitle('Reverse Reads')
```



``` r
# Save the plot
ggsave('viz_data/plotQualityProfile.png', plot = p1 + p2, width = 18, height = 12, units = 'in', dpi = 600)
```

We can manually visualize our quality profiles for all of the samples (plotQualityProfile.png) in the viz_data folder. Quality scores typically decrease over sequencing cycles (positions along the read) due to several factors:

1. **Chemistry degradation**: As sequencing progresses, the chemistry becomes less efficient and error rates increase
2. **Signal interference**: Accumulation of incomplete washes, phasing/pre-phasing, and fluorescent noise
3. **Base calling uncertainty**: Higher error rates in later cycles make base calling less confident


``` r
# View quality profiles of reverse reads for the first 2 samples
plotQualityProfile(fnRs[1:2])
```

![](../Submodule02_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

Reverse reads typically show lower quality scores than forward reads. This is normal in paired-end sequencing data for several reasons:

1. **Sequential Nature**: Reverse reads are generated after forward ones, when:
   - Reagents are more depleted
   - More sequencing cycles have occurred
   - More errors have accumulated

2. **DADA2's Handling**: DADA2 is designed to work with this quality variation because:
   - It incorporates quality scores into its error models
   - It can differentiate between sequencing errors and real biological variation
   - Each base is weighted by its quality score in the algorithm

3. **Trimming Strategy**: While DADA2 can handle lower quality scores, we still trim the poorest quality bases because:
   - It improves the algorithm's performance
   - Reduces computational time
   - Increases confidence in the final results
   - Helps remove technical artifacts


<br>

### Understanding Our Sequencing Run Quality

Due to a recent system update on the MiSeq, the sequencing run was inadvertently over-clustered, producing an excess of low-quality reads we see in the plots. Let's break down what happened and why it matters:
 
#### What is Over-Clustering?
During DNA sequencing, molecules attach to the flow cell surface in spots called "clusters". Think of it like planting seeds in a garden:
- Optimal spacing: Each plant has room to grow and be identified clearly
- Over-crowded: Plants grow too close, making it hard to tell them apart
 
#### What Happened in Our Run?
A MiSeq system update led to over-clustering, meaning:
- Too many DNA clusters formed on the flow cell
- Clusters grew too close together
- Their signals interfered with each other
- Like trying to read many overlapping words
 
#### Impact on Our Data
1. **Initial Data:**
   - High number of raw reads
   - But many reads have poor quality scores
  
2. **After Quality Filtering:**
   - Significant reduction in read count
   - Remaining reads are high quality
   - More trustworthy for analysis
 
This is okay for our analysis because quality is more important than quantity, removing unreliable data improves our final results, and we still have sufficient good quality reads for analysis.


### Here is an example of good quality reads: 
High-quality reads typically have consistent Phred scores across their length, with minimal drop-off at the ends.


![Forward Reads Quality Profile](../images/plotFs.png)

![Reverse Reads Quality Profile](../images/plotRs.png)


### Filter and Trim



``` r
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- names(filtRs) <- sample.names
```


DADA2's `filterAndTrim()` function processes raw sequencing reads to ensure high-quality data for analysis. Here's a explanation of each parameter:

#### Key Parameters for V4 16S rRNA Data:

**`truncLen=c(200,130)`**: 
- For V4 region (~250 bp), we expect ~200-250 bp reads
- Forward reads (R1) at 200 bp because quality usually stays good
- Reverse reads (R2) truncated more aggressively at 130 bp due to quality drop
- These numbers are chosen by examining quality plots for your specific run
- Must ensure enough overlap between F and R reads (~20 bp minimum)

**Quality Control Parameters:**
- `maxEE=c(2,5)`: Maximum "expected errors" allowed
  - More stringent for forward reads (2 errors)
  - More lenient for reverse reads (5 errors)
  - Higher numbers = more lenient filtering

- `truncQ=2`: Cuts reads when quality score drops to Q2
  - Removes poor quality ends dynamically
  - Q2 is very low quality, ensures removal of worst bases

**Technical Parameters:**
- `maxN=0`: Removes reads with any ambiguous bases (N)
- `trimLeft=c(0,0)`: No trimming from read starts
- `rm.phix=TRUE`: Removes PhiX spike-in control DNA (if you are using those)

#### How to Choose Parameters:
1. Check quality plots first
2. Look for where quality drops significantly
3. Ensure enough overlap remains for merging
4. Consider your downstream analysis needs
5. Balance stringency with maintaining sufficient read depth


``` r
# Set trimming standards first
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,130),
                     maxN=0, maxEE=c(2,5), trimLeft = c(0,0), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
```



``` r
# View outputs
head(out)
```

```
##                                reads.in reads.out
## WP-11O29_S66_L001_R1_001.fastq   114247     39363
## WP-1641C_S68_L001_R1_001.fastq   180306     63123
## WP-1J0A0_S50_L001_R1_001.fastq   118458     41544
## WP-3GH4E_S51_L001_R1_001.fastq   105493     37698
## WP-4SP0U_S57_L001_R1_001.fastq   245187     83598
## WP-6RCV3_S39_L001_R1_001.fastq   114269     40530
```


### Learn the Error Rates

While quality scores tell us about the confidence in each base call, DADA2 goes further by learning the specific error patterns in your dataset. Here's how it works:

#### Quality Scores vs. Error Models
- **Quality scores** are like confidence ratings for each base
  - Assigned during sequencing
  - Based on how clear the signal was
  - Same score means same confidence across any run

- **DADA2's error models** are like fingerprints of your specific run
  - Learns what errors look like in your data
  - Accounts for systematic patterns specific to your run
  - Different for every dataset

#### How DADA2 Learns
1. First looks at your sequences and their frequencies
2. Makes initial guesses about which sequences are real vs. errors
3. Keeps refining these guesses as it processes more data
4. Creates a customized error profile for your specific dataset

Think of it like this: Quality scores are like spell-checking individual words, while DADA2's error model is like learning someone's specific typing patterns to better catch their common mistakes.

*Note: This is why DADA2 is more sophisticated than simple quality filtering - it understands the unique "personality" of your sequencing run's errors.*

*These steps will take several minutes to run.*


``` r
# Error model for forward reads
errF <- learnErrors(filtFs, multithread=TRUE)
```

```
## 108897000 total bases in 544485 reads from 9 samples will be used for learning the error rates.
```



``` r
# Error model for reverse reads
errR <- learnErrors(filtRs, multithread=TRUE)
```

```
## 104995670 total bases in 807659 reads from 13 samples will be used for learning the error rates.
```



``` r
# Save our outputs
saveRDS(errF,file="errF.rds")
saveRDS(errR,file="errR.rds")
```


<details>
<summary>If you are coming back to this module, reload your saved data here.</summary>
<br>
If you are coming back to this module, you can load your saved data instead of rerunning. Use the following R code to upload saved data.
    

``` r
errF <- readRDS("errF.rds")
errR <- readRDS("errR.rds")
```
</details>


We can now plot the estimated error rates. These plots displays the observed error rates across various quality scores for each base position in the reads (A, C, G, T). It allows us to compare the empirical error rates (from the actual data) to the estimated error rates (based on the model DADA2 has learned represented in a black line).


``` r
# Plot forward reads estimated error rates
plotErrors(errF, nominalQ=TRUE) +
  ggtitle('errF')
```

```
## Warning in scale_y_log10(): log-10 transformation introduced infinite values.
```

```
## Warning: Removed 82 rows containing missing values or values outside the scale range
## (`geom_line()`).
## Removed 82 rows containing missing values or values outside the scale range
## (`geom_line()`).
```

![](../Submodule02_files/figure-html/unnamed-chunk-20-1.png)<!-- -->



``` r
# Plot reverse reads estimated error rates
plotErrors(errR, nominalQ = TRUE) +
  ggtitle('errR')
```

```
## Warning: Removed 82 rows containing missing values or values outside the scale range
## (`geom_line()`).
## Removed 82 rows containing missing values or values outside the scale range
## (`geom_line()`).
```

![](../Submodule02_files/figure-html/unnamed-chunk-21-1.png)<!-- -->


### Sample Inference

We will now use DADA2's learned error rates and identify amplicon sequence variants (ASV) while removing likely errors. 

*This step will take several minutes to run.*


``` r
# Creating DADA object
dadaFs <- dada(filtFs, err=errF, multithread=TRUE) # filtered reads in unique sequences
```

```
## Sample 1 - 39363 reads in 14187 unique sequences.
## Sample 2 - 63123 reads in 26683 unique sequences.
## Sample 3 - 41544 reads in 17784 unique sequences.
## Sample 4 - 37698 reads in 12290 unique sequences.
## Sample 5 - 83598 reads in 32871 unique sequences.
## Sample 6 - 40530 reads in 14527 unique sequences.
## Sample 7 - 99030 reads in 42889 unique sequences.
## Sample 8 - 39699 reads in 13841 unique sequences.
## Sample 9 - 99900 reads in 39434 unique sequences.
## Sample 10 - 110089 reads in 40284 unique sequences.
## Sample 11 - 75378 reads in 28399 unique sequences.
## Sample 12 - 36181 reads in 12523 unique sequences.
## Sample 13 - 41526 reads in 16087 unique sequences.
## Sample 14 - 31948 reads in 11354 unique sequences.
## Sample 15 - 36964 reads in 13680 unique sequences.
## Sample 16 - 19481 reads in 9230 unique sequences.
## Sample 17 - 77711 reads in 29452 unique sequences.
## Sample 18 - 114355 reads in 42651 unique sequences.
## Sample 19 - 49661 reads in 20678 unique sequences.
## Sample 20 - 45830 reads in 17676 unique sequences.
## Sample 21 - 122888 reads in 46755 unique sequences.
## Sample 22 - 103522 reads in 36545 unique sequences.
## Sample 23 - 106278 reads in 44948 unique sequences.
## Sample 24 - 45389 reads in 18754 unique sequences.
## Sample 25 - 131341 reads in 55136 unique sequences.
## Sample 26 - 49243 reads in 20042 unique sequences.
## Sample 27 - 29334 reads in 10554 unique sequences.
## Sample 28 - 104203 reads in 42992 unique sequences.
## Sample 29 - 43671 reads in 16335 unique sequences.
## Sample 30 - 97132 reads in 41151 unique sequences.
## Sample 31 - 69682 reads in 25638 unique sequences.
## Sample 32 - 72768 reads in 28278 unique sequences.
## Sample 33 - 97360 reads in 33458 unique sequences.
## Sample 34 - 56868 reads in 21112 unique sequences.
## Sample 35 - 77521 reads in 28580 unique sequences.
```

``` r
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

```
## Sample 1 - 39363 reads in 15064 unique sequences.
## Sample 2 - 63123 reads in 27307 unique sequences.
## Sample 3 - 41544 reads in 16990 unique sequences.
## Sample 4 - 37698 reads in 12863 unique sequences.
## Sample 5 - 83598 reads in 33883 unique sequences.
## Sample 6 - 40530 reads in 13599 unique sequences.
## Sample 7 - 99030 reads in 35760 unique sequences.
## Sample 8 - 39699 reads in 14950 unique sequences.
## Sample 9 - 99900 reads in 37298 unique sequences.
## Sample 10 - 110089 reads in 37099 unique sequences.
## Sample 11 - 75378 reads in 36755 unique sequences.
## Sample 12 - 36181 reads in 10374 unique sequences.
## Sample 13 - 41526 reads in 13581 unique sequences.
## Sample 14 - 31948 reads in 10654 unique sequences.
## Sample 15 - 36964 reads in 11463 unique sequences.
## Sample 16 - 19481 reads in 8674 unique sequences.
## Sample 17 - 77711 reads in 31860 unique sequences.
## Sample 18 - 114355 reads in 32740 unique sequences.
## Sample 19 - 49661 reads in 20391 unique sequences.
## Sample 20 - 45830 reads in 14569 unique sequences.
## Sample 21 - 122888 reads in 40315 unique sequences.
## Sample 22 - 103522 reads in 37254 unique sequences.
## Sample 23 - 106278 reads in 38028 unique sequences.
## Sample 24 - 45389 reads in 17592 unique sequences.
## Sample 25 - 131341 reads in 43710 unique sequences.
## Sample 26 - 49243 reads in 16616 unique sequences.
## Sample 27 - 29334 reads in 10477 unique sequences.
## Sample 28 - 104203 reads in 37545 unique sequences.
## Sample 29 - 43671 reads in 17089 unique sequences.
## Sample 30 - 97132 reads in 39960 unique sequences.
## Sample 31 - 69682 reads in 25054 unique sequences.
## Sample 32 - 72768 reads in 26365 unique sequences.
## Sample 33 - 97360 reads in 44237 unique sequences.
## Sample 34 - 56868 reads in 18983 unique sequences.
## Sample 35 - 77521 reads in 32367 unique sequences.
```



``` r
# Review DADA object
dadaFs[[1]]
```

```
## dada-class: object describing DADA2 denoising results
## 144 sequence variants were inferred from 14187 input unique sequences.
## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16
```


### Merging Paired Reads

Our data is paired-end sequencing, DNA fragments are sequenced from both ends - the forward read (5' end) and reverse read (3' end) of the template strand. This approach provides higher quality data than single-end sequencing because the overlapping regions between reads allow for error validation and longer final sequences. DADA2 merges these paired reads by identifying matching overlap regions and combining them into single, high-confidence sequences, discarding any pairs that don't align properly.

The `mergePairs()` funcation requires the dada objects we created as well as the filtered reads, because not all of that information is stored in the dada objects we created in the last step. 


``` r
# Combine the forward and reverse reads into one
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

```
## 35785 paired-reads (in 138 unique pairings) successfully merged out of 38686 (in 999 pairings) input.
```

```
## 53717 paired-reads (in 446 unique pairings) successfully merged out of 60741 (in 3081 pairings) input.
```

```
## 36690 paired-reads (in 322 unique pairings) successfully merged out of 40576 (in 2088 pairings) input.
```

```
## 35362 paired-reads (in 134 unique pairings) successfully merged out of 37119 (in 738 pairings) input.
```

```
## 71673 paired-reads (in 622 unique pairings) successfully merged out of 81328 (in 4326 pairings) input.
```

```
## 37855 paired-reads (in 181 unique pairings) successfully merged out of 40114 (in 871 pairings) input.
```

```
## 85388 paired-reads (in 691 unique pairings) successfully merged out of 97325 (in 5002 pairings) input.
```

```
## 36520 paired-reads (in 216 unique pairings) successfully merged out of 39028 (in 1120 pairings) input.
```

```
## 87580 paired-reads (in 632 unique pairings) successfully merged out of 97779 (in 4232 pairings) input.
```

```
## 96754 paired-reads (in 712 unique pairings) successfully merged out of 108420 (in 4084 pairings) input.
```

```
## 64124 paired-reads (in 380 unique pairings) successfully merged out of 73004 (in 2917 pairings) input.
```

```
## 33886 paired-reads (in 190 unique pairings) successfully merged out of 35589 (in 836 pairings) input.
```

```
## 37238 paired-reads (in 301 unique pairings) successfully merged out of 40612 (in 1776 pairings) input.
```

```
## 29704 paired-reads (in 94 unique pairings) successfully merged out of 31397 (in 581 pairings) input.
```

```
## 33822 paired-reads (in 221 unique pairings) successfully merged out of 36347 (in 1036 pairings) input.
```

```
## 17383 paired-reads (in 158 unique pairings) successfully merged out of 18865 (in 745 pairings) input.
```

```
## 67794 paired-reads (in 299 unique pairings) successfully merged out of 76338 (in 2477 pairings) input.
```

```
## 103320 paired-reads (in 491 unique pairings) successfully merged out of 112706 (in 3154 pairings) input.
```

```
## 43015 paired-reads (in 299 unique pairings) successfully merged out of 48412 (in 2035 pairings) input.
```

```
## 41990 paired-reads (in 154 unique pairings) successfully merged out of 45184 (in 1202 pairings) input.
```

```
## 106122 paired-reads (in 790 unique pairings) successfully merged out of 120866 (in 4848 pairings) input.
```

```
## 90427 paired-reads (in 582 unique pairings) successfully merged out of 101653 (in 3410 pairings) input.
```

```
## 93042 paired-reads (in 783 unique pairings) successfully merged out of 104018 (in 4968 pairings) input.
```

```
## 41204 paired-reads (in 228 unique pairings) successfully merged out of 44584 (in 1520 pairings) input.
```

```
## 110386 paired-reads (in 1180 unique pairings) successfully merged out of 128944 (in 7236 pairings) input.
```

```
## 44761 paired-reads (in 237 unique pairings) successfully merged out of 48302 (in 1502 pairings) input.
```

```
## 27224 paired-reads (in 120 unique pairings) successfully merged out of 28749 (in 696 pairings) input.
```

```
## 88260 paired-reads (in 754 unique pairings) successfully merged out of 102142 (in 4912 pairings) input.
```

```
## 39647 paired-reads (in 254 unique pairings) successfully merged out of 42872 (in 1303 pairings) input.
```

```
## 82639 paired-reads (in 576 unique pairings) successfully merged out of 94995 (in 4675 pairings) input.
```

```
## 62969 paired-reads (in 352 unique pairings) successfully merged out of 68743 (in 2220 pairings) input.
```

```
## 65246 paired-reads (in 343 unique pairings) successfully merged out of 71480 (in 2131 pairings) input.
```

```
## 83187 paired-reads (in 457 unique pairings) successfully merged out of 95187 (in 3494 pairings) input.
```

```
## 50839 paired-reads (in 256 unique pairings) successfully merged out of 55903 (in 1585 pairings) input.
```

```
## 68652 paired-reads (in 309 unique pairings) successfully merged out of 76159 (in 2253 pairings) input.
```

``` r
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
```

```
##                                                                                                                                                                                                                                                        sequence
## 1 TACGGAAGGTCCGGGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGCCGGAGATTAAGCGTGTTGTGAAATGTAGACGCTCAACGTCTGCACTGCAGCGCGAACTGGTTTCCTTGAGTACGCACAAAGTGGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCTCACTGGAGCGCAACTGACGCTGAAGCTCGAAAGTGCGGGTATCGAACAGG
## 2 TACGTAGGGTGCAAGCGTTATCCGGAATTATTGGGCGTAAAGGGCTCGTAGGCGGTTCGTCGCGTCCGGTGTGAAAGTCCATCGCTTAACGGTGGATCCGCGCCGGGTACGGGCGGGCTTGAGTGCGGTAGGGGAGACTGGAATTCCCGGTGTAACGGTGGAATGTGTAGATATCGGGAAGAACACCAATGGCGAAGGCAGGTCTCTGGGCCGTCACTGACGCTGAGGAGCGAAAGCGTGGGGAGCGAACAGG
## 3 TACGGAAGGTCCGGGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGCCGGAGATTAAGCGTGTTGTGAAATGTAGATGCTCAACATCTGAACTGCAGCGCGAACTGGTTTCCTTGAGTACGCACAAAGTGGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCTCACTGGAGCGCAACTGACGCTGAAGCTCGAAAGTGCGGGTATCGAACAGG
## 4 TACGGAAGGTCCGGGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGCCGGAGATTAAGCGTGTTGTGAAATGTAGATGCTCAACATCTGCACTGCAGCGCGAACTGGTTTCCTTGAGTACGCACAAAGTGGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCTCACTGGAGCGCAACTGACGCTGAAGCTCGAAAGTGCGGGTATCGAACAGG
## 5 TACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGGTGTGGCAAGTCTGATGTGAAAGGCATGGGCTCAACCTGTGGACTGCATTGGAAACTGTCATACTTGAGTGCCGGAGGGGTAAGCGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAACACCAGTGGCGAAGGCGGCTTACTGGACGGTAACTGACGTTGAGGCTCGAAAGCGTGGGGAGCAAACAGG
## 6 AACGTAGGTCACAAGCGTTGTCCGGAATTACTGGGTGTAAAGGGAGCGCAGGCGGGAAGACAAGTTGGAAGTGAAATCTATGGGCTCAACCCATAAACTGCTTTCAAAACTGTTTTTCTTGAGTAGTGCAGAGGTAGGCGGAATTCCCGGTGTAGCGGTGGAATGCGTAGATATCGGGAGGAACACCAGTGGCGAAGGCGGCCTACTGGGCACCAACTGACGCTGAGGCTCGAAAGTGTGGGTAGCAAACAGG
##   abundance forward reverse nmatch nmismatch nindel prefer accept
## 1      9062       1       1     77         0      0      2   TRUE
## 2      2820       2       3     77         0      0      1   TRUE
## 3      2214       3       1     77         0      0      2   TRUE
## 4      2132       4       1     77         0      0      2   TRUE
## 5      1975       5       4     77         0      0      1   TRUE
## 6      1379       6       2     77         0      0      2   TRUE
```

### Construct sequence table

We can now constuct the ASV table. This table is a matrix of samples (rows) and each ASV identified (columns). 


``` r
# View the dimensions of our merged samples 
seqtab  <- makeSequenceTable(mergers)
dim(seqtab)
```

```
## [1]   35 5902
```

Our table contains 35 samples and 5902 ASVs. Now let's view the lengths of the variants. 


``` r
# Inspect the distribution of sequence lengths
table(nchar(getSequences(seqtab))) # the length should be around 250
```

```
## 
##  200  252  253  254  255 
##    1  389 5466   44    2
```


### Remove Chimeras

DADA2's error model corrects for substitutions and indels, but not Chimeras (two different DNA strands megred together). Chimeras are easier to identify in ASVs, as opposed to OTUs, because of their exact sequences.


``` r
# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```

```
## Identified 4456 bimeras out of 5902 input sequences.
```

*Note: DADA2 uses the term "bimera" specifically because it refers to sequences formed from exactly *two* parent sequences, "bi-" meaning two. While "chimera" is the broader term that can include sequences formed from multiple parent sequences (three or more), most artificial chimeric sequences in amplicon sequencing are actually bimeras. DADA2's algorithm specifically looks for these two-parent artifacts as they are the most common type of chimeric sequence formed during PCR amplification.*


``` r
# New sequence table
dim(seqtab.nochim)
```

```
## [1]   35 1446
```

``` r
# Divide the new number of ASVs by the first table.
# This shows us what percentage was removed
sum(seqtab.nochim)/sum(seqtab)
```

```
## [1] 0.9041567
```

Here we can see that Chimeras accounted for ~10% of the merged sequence reads.


``` r
# Save our data
saveRDS(seqtab.nochim, "seqtabnochim.rds")
```

### Track reads through the pipeline

Let's construct a summary table of our steps to see how many reads were removed in the cleaning process. 


``` r
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim),
               final_perc_reads_retained=round(rowSums(seqtab.nochim)/out[,1]*100, 1))
# formatting table before writing to file
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim","percent")
rownames(track) <- sample.names
track
```

```
##           input filtered denoisedF denoisedR merged nonchim percent
## WP-11O29 114247    39363     39130     38899  35785   35325    30.9
## WP-1641C 180306    63123     62576     61262  53717   50893    28.2
## WP-1J0A0 118458    41544     41085     41021  36690   35654    30.1
## WP-3GH4E 105493    37698     37503     37303  35362   34799    33.0
## WP-4SP0U 245187    83598     82881     82013  71673   64365    26.3
## WP-6RCV3 114269    40530     40341     40292  37855   35649    31.2
## WP-7FAEC 288105    99030     98289     98026  85388   75749    26.3
## WP-8SXNH 115689    39699     39364     39339  36520   35137    30.4
## WP-8V6E3 286174    99900     99028     98606  87580   78614    27.5
## WP-95020 314274   110089    109446    109042  96754   81235    25.8
## WP-9ZJJI 214646    75378     74726     73625  64124   57152    26.6
## WP-ATJLX 101691    36181     35970     35790  33886   32929    32.4
## WP-DFG11 118093    41526     41037     41058  37238   36165    30.6
## WP-EC5AP  90638    31948     31752     31586  29704   29463    32.5
## WP-G375R 105928    36964     36671     36630  33822   32454    30.6
## WP-G9EBL  57168    19481     19278     19033  17383   17321    30.3
## WP-J06MK 223327    77711     77295     76740  67794   63284    28.3
## WP-J0B3X 335549   114355    113660    113368 103320   91594    27.3
## WP-KM8BK 146122    49661     49256     48787  43015   39750    27.2
## WP-KXGT6 132046    45830     45486     45502  41990   40801    30.9
## WP-MRCZG 352149   122888    121987    121730 106122   90633    25.7
## WP-QIGEC 302735   103522    102921    102225  90427   77523    25.6
## WP-RAJES 307883   106278    105429    104844  93042   80190    26.0
## WP-RPV50 131966    45389     45066     44899  41204   40710    30.8
## WP-SUFNU 383563   131341    130268    129981 110386   88913    23.2
## WP-TDY7G 146345    49243     48841     48686  44761   43235    29.5
## WP-U4XI6  85807    29334     29158     28906  27224   26929    31.4
## WP-UIXW3 304129   104203    103303    103001  88260   75767    24.9
## WP-X576L 122710    43671     43368     43158  39647   37792    30.8
## WP-X803O 281887    97132     96371     95702  82639   74452    26.4
## WP-XDI0I 199408    69682     69310     69088  62969   56107    28.1
## WP-YF6XT 210492    72768     72192     72016  65246   59099    28.1
## WP-YRYYF 273738    97360     96780     95716  83187   71752    26.2
## WP-YY8WG 165937    56868     56492     56250  50839   47555    28.7
## WP-Z01LM 229511    77521     77018     76632  68652   63541    27.7
```

``` r
write.table(track, "read-count-tracking.txt", quote=FALSE, sep="\t")
```

We retained ~25-30% of our original reads! This is likely due to the poor quality we saw at the beginning. Although we now have much less data to work with, it is high quality. 

### Read in taxonomy

The last step in this model is to read in a taxonomy database. We must first select a database. Two of the most commonly used are Greengenes and SILVA, each with its strengths and weaknesses.
 + **Greengenes<sup>2</sup>** is a chimera-checked 16S rRNA gene database. It provides chimera screening, standard alignment, and taxonomic classification using multiple published taxonomies. Its most recent update was on July 03, 2017.
 + **[SILVA](https://www.arb-silva.de/documentation/release-138.1/)<sup>3-4</sup>** is a high quality ribosomal RNA database. It is a comprehensive and quality-controlled database for up-to-date ribosomal RNA sequences. Additionally, SILVA also provides many other tools like alignment, phylogenetic tree classification, and probe/primer matching. It was last updated on August 27, 2020.

Since our data has chimeras removed and we want the most up-to-date analysis, we will use SILVA.

*This step will take several minutes to run.*


``` r
# Read in Taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "tax/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "tax/silva_species_assignment_v138.1.fa.gz")
```


``` r
# Let's take a look
taxa.print <- taxa # Removing sequence row names for display only
rownames(taxa.print) <- NULL
head(taxa.print)
```

```
##      Kingdom    Phylum         Class         Order            
## [1,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales"  
## [2,] "Bacteria" "Firmicutes"   "Clostridia"  "Oscillospirales"
## [3,] "Bacteria" "Firmicutes"   "Clostridia"  "Lachnospirales" 
## [4,] "Bacteria" "Firmicutes"   "Clostridia"  "Lachnospirales" 
## [5,] "Bacteria" "Firmicutes"   "Clostridia"  "Oscillospirales"
## [6,] "Bacteria" "Firmicutes"   "Clostridia"  "Lachnospirales" 
##      Family            Genus              Species         
## [1,] "Bacteroidaceae"  "Bacteroides"      "vulgatus"      
## [2,] "Ruminococcaceae" "Faecalibacterium" NA              
## [3,] "Lachnospiraceae" "Blautia"          NA              
## [4,] "Lachnospiraceae" "Agathobacter"     NA              
## [5,] "Ruminococcaceae" "Faecalibacterium" "prausnitzii"   
## [6,] "Lachnospiraceae" "Fusicatenibacter" "saccharivorans"
```

#### Adding Controls

Microbime studies often include mock and/or blank controls to help identify and mitigate possible contamination, ensuring the accuracy and reliability of the data.

 + **Mock Controls:** These are samples with a known microbial community composition, often with a defined mix of microbial strains in known proportions. By including mock controls, researchers can verify the performance of the sequencing pipeline and taxonomic assignment tools. If the analysis accurately identifies the species and proportions in the mock community, it provides confidence that the methods are performing correctly. Any discrepancies may indicate biases in the sequencing or analysis process that should be addressed.
 + **Blank Controls:** These are negative controls that contain no added DNA. Blank controls help detect potential contamination from reagents, the lab environment, or equipment. Since microbial DNA is often present at very low concentrations, contamination can skew results, especially for low-biomass samples. If microbial sequences appear in blank controls, it could indicate contamination, and researchers can account for this by removing these sequences from experimental samples.

Although our study does not have these controls, incorporating both mock and blank controls in the study design can boost quality control. It helps distinguish true biological signals from artifacts introduced during sample processing or sequencing. 


``` r
# Save our data
saveRDS(taxa, "taxa.rds")
```


## Conclusion
In this module, we took a comprehensive journey from the initial stages of data acquisition through taxonomic assignment. We began by downloading and exploring the survey data to understand the metadata associated with each sample. Following this, we downloaded the raw sequencing FASTQ files and preprocessed them using the DADA2 pipeline, which allowed us to perform quality control, error correction, and sequence denoising, resulting in high-quality amplicon sequence variants (ASVs). Finally, we assigned taxonomy to the ASVs using a reference database, enabling us to identify and categorize the microbial communities present in our samples. This workflow provides a foundational dataset ready for downstream analysis, such as community composition comparisons and diversity assessments in the next module.


## References

1. Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP. 2016. DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13:581–583. doi:10.1038/nmeth.3869
2. DeSantis TZ, Hugenholtz P, Larsen N, Rojas M, Brodie EL, Keller K, Huber T, Dalevi D, Hu P, Andersen GL. Greengenes, a chimera-checked 16S rRNA gene database and workbench compatible with ARB. Appl Environ Microbiol. 2006 Jul;72(7):5069-72. doi: 10.1128/AEM.03006-05. PMID: 16820507; PMCID: PMC1489311.
3. Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO (2013) The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucl. Acids Res. 41 (D1): D590-D596.
4. Yilmaz P, Parfrey LW, Yarza P, Gerken J, Pruesse E, Quast C, Schweer T, Peplies J, Ludwig W, Glöckner FO (2014) The SILVA and "All-species Living Tree Project (LTP)" taxonomic frameworks. Nucl. Acids Res. 42:D643-D648
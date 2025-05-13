# Submodule 3: Taxonomic Classification and Diversity Analysis

## Overview

The module focuses on analyzing 16S rRNA sequencing data to explore the composition and diversity of the gut microbiome. The module begins with importing and integrating survey-based metadata, followed by assigning taxonomy to sequencing reads using the phyloseq package. Key concepts such as rarefaction, alpha diversity (e.g., Shannon and Inverse Simpson indices), and beta diversity (e.g., Bray-Curtis NMDS) are introduced to assess microbial diversity within and between samples. Visualizations and statistical analyses are used to examine the relationships between microbiome composition and variables like protein intake, sex, and BMI. Through hands-on exercises, learners gain practical experience in microbiome data processing, visualization, and interpretation, building a strong foundation for future research.

<center>
    <img src="../images/3.jpg" alt="2" width="500"/>
</center>



## Learning Objectives
+ Perform taxonomic classification of 16S rRNA sequences
+ Evaluate and explore rarefaction curves to assess sampling depth and species richness
+ Analyze alpha and beta diversity within microbial communities



# 1. Import Data

### Files
 
 <b>1. Food Frequency Questionnaire (FFQ):</b>
    - This CSV file contains data from a food frequency questionnaire completed by study participants. It includes columns for participant identifiers (e.g., SampleID), demographic information (e.g., SEX, AGE), and food consumption frequencies and quantities (e.g., BREAKFASTSANDWICHFREQ, EGGSFREQ). This data will allow us to analyze dietary habits and link them with microbiome profiles.
    
    
 <b>2. Lifestyle Questionnaire:</b>
    - This text file contains responses to lifestyle questions, such as physical activity levels, smoking status, and other lifestyle factors. It has a row for each question, an answer key, and responses for individual participants in columns.
    
   
 <b>3.  FASTQ Files:</b>
    - These files contain raw sequencing data from 16S rRNA sequencing. Each read in a FASTQ file has a sequence identifier, nucleotide sequence, and quality scores. 
    

Let's start by loading the survey files into DataFrames and performing some initial quality checks to ensure data integrity.


## Import Survey Data

### Step 1: Load the Survey Data
Let's load both survey files:


``` r
# Install and load required packages
if (!requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
}

if (!requireNamespace("vegan", quietly = TRUE)) {
    install.packages("vegan")
}

if (!requireNamespace("ggpubr", quietly = TRUE)) {
    install.packages("ggpubr")
}
```


``` r
library("dplyr")
library("ggplot2")
library("vegan")
library("ggpubr")
```


``` r
#Load in the data
ffq <- read.csv("metadata/FFQ_Data.csv")
lsq <- read.csv("metadata/MetaData.txt", sep = '\t')
# Display the first few rows of the dataframe to understand its structure
head(ffq[, 1:5])
```

```
##   SampleID RESPONDENTID BOOKNUM SEX PREGNANT
## 1 WP_J06MK            1   18980   1        M
## 2 WP_SUFNU            2   18979   2        1
## 3 WP_1VMLX            4   18978   1        M
## 4 WP_J0B3X            5   18976   1        3
## 5 WP_RAJES            6   18977   2        1
## 6 WP_G9EBL            7   18971   1        3
```

Our Food Frequency Questionaire has columns representing the Sample ID, Respondent ID number, and all the questions asked. The heading of the columns represents abbreviations for the questions. Now let's look at the dimensions of our data.


``` r
# Full dimensions
dim(ffq)
```

```
## [1]  100 1074
```

``` r
# Number of row
nrow(ffq)
```

```
## [1] 100
```

``` r
# Number of columns
ncol(ffq)
```

```
## [1] 1074
```

We can see that we have 100 participants (every row is a sample) and 1074 columns which means over 1000 food questions! Now let's look at the lifestyle survey. 


``` r
# Display the first few rows of the dataframe to understand its structure
head(lsq[, 1:5])
```

```
##   SampleID    Sex Age Weight Height
## 1 WP_J06MK   Male  34    241     69
## 2 WP_SUFNU Female  36    120     52
## 3 WP_1VMLX   Male  21    130     65
## 4 WP_J0B3X   Male  27    175     72
## 5 WP_RAJES Female  57    140     61
## 6 WP_G9EBL   Male  56    150     68
```

Our Lifestyle Questionaire has columns representing the Sample ID and all the lifestyle questions asked. The heading of the columns represents abbreviations for the questions. Now let's look at the dimensions of our data.



``` r
# Full dimensions
dim(lsq)
```

```
## [1] 96 33
```

``` r
# Number of row
nrow(lsq)
```

```
## [1] 96
```

``` r
# Number of columns
ncol(lsq)
```

```
## [1] 33
```

We have less questions in this survey, only 33, and 96 participants.


### Step 2: Quick Quality Control Checks
Let's perform several initial checks to assess data quality:


  <b>1. Collect Intersection of Samples:</b> We need to find the Sample IDs in the FFQ data rows and in the LSQ data columns. 
 


``` r
# Find sample intersection
samples <- intersect(ffq$SampleID, lsq$SampleID)
```


``` r
length(samples)
```

```
## [1] 96
```

  <b>2. Basic Summary Statistics:</b> We can get an overview of the data, which can help identify outliers or unexpected values. We will look at the demographic questions in the FFQ and Lifestyle Questionnaire separately.
 
 

``` r
# Summary statistics for the first 10 columns
summary(ffq[, 1:9])
```

```
##    SampleID          RESPONDENTID       BOOKNUM           SEX      
##  Length:100         Min.   :  1.00   Min.   :18942   Min.   :1.00  
##  Class :character   1st Qu.: 26.75   1st Qu.:18977   1st Qu.:1.00  
##  Mode  :character   Median : 52.50   Median :19026   Median :2.00  
##                     Mean   : 52.06   Mean   :19025   Mean   :1.63  
##                     3rd Qu.: 77.25   3rd Qu.:19066   3rd Qu.:2.00  
##                     Max.   :102.00   Max.   :19139   Max.   :2.00  
##                                                                    
##    PREGNANT              AGE            WEIGHT        HEIGHTFEET  
##  Length:100         Min.   :18.00   Min.   : 87.0   Min.   :5.00  
##  Class :character   1st Qu.:25.25   1st Qu.:137.5   1st Qu.:5.00  
##  Mode  :character   Median :40.50   Median :150.0   Median :5.00  
##                     Mean   :41.31   Mean   :161.1   Mean   :5.18  
##                     3rd Qu.:54.75   3rd Qu.:180.0   3rd Qu.:5.00  
##                     Max.   :85.00   Max.   :268.0   Max.   :6.00  
##                     NA's   :2       NA's   :2                     
##   HEIGHTINCHES  
##  Min.   : 0.00  
##  1st Qu.: 2.00  
##  Median : 6.00  
##  Mean   : 5.42  
##  3rd Qu.: 8.00  
##  Max.   :11.00  
## 
```
 

``` r
# Summary statistics for 3 columns
summary(lsq[, 3:5])
```

```
##       Age            Weight          Height     
##  Min.   :19.00   Min.   : 97.0   Min.   :52.00  
##  1st Qu.:26.00   1st Qu.:136.8   1st Qu.:64.00  
##  Median :40.50   Median :150.0   Median :67.00  
##  Mean   :41.62   Mean   :161.4   Mean   :66.90  
##  3rd Qu.:55.00   3rd Qu.:181.2   3rd Qu.:69.62  
##  Max.   :85.00   Max.   :268.0   Max.   :76.00
```
 

  <b>3. Visualize Data:</b> Ensure the data looks accurate with exploratory visualizations.


``` r
# Load necessary library
library(ggplot2)

# Set up the plotting theme
theme_set(theme_minimal())

# 1. Bar plot for Sex distribution
# Bar plot for Sex distribution without modifying the original data
ggplot(ffq, aes(x = factor(SEX, levels = c(1, 2), labels = c("Male", "Female")))) +
  geom_bar(fill = "#0073C2FF") +  # Blue color for bars
  labs(title = "Distribution of Sex", x = "Sex", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 20),
    axis.text = element_text(size = 18))
```

![](../Submodule03_files/figure-html/unnamed-chunk-11-1.png)<!-- -->


``` r
ggplot(lsq, aes(x = factor(Sex))) +
  geom_bar(fill = "#0073C2FF") +  # Blue color for bars
  labs(title = "Distribution of Sex", x = "Sex", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 20),
        axis.text = element_text(size = 18))
```

![](../Submodule03_files/figure-html/unnamed-chunk-12-1.png)<!-- -->


``` r
# 2. Histogram for Age distribution
# Dropping NA values for AGE
ggplot(ffq, aes(x = AGE)) +
  geom_histogram(bins = 10, fill = "lightgreen", color = "black") +
  labs(title = "Age Distribution", x = "Age", y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 20),
    axis.text = element_text(size = 18))
```

![](../Submodule03_files/figure-html/unnamed-chunk-13-1.png)<!-- -->


``` r
ggplot(lsq, aes(x = Age)) +
  geom_histogram(bins = 10, fill = "lightgreen", color = "black") +
  labs(title = "Age Distribution", x = "Age", y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 20),
        axis.text = element_text(size = 18))
```

![](../Submodule03_files/figure-html/unnamed-chunk-14-1.png)<!-- -->


``` r
# 3. Histogram for Weight distribution
# Dropping NA values for WEIGHT
ggplot(ffq, aes(x = WEIGHT)) +
  geom_histogram(bins = 15, fill = "lightcoral", color = "black") +
  labs(title = "Weight Distribution", x = "Weight", y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 20),
    axis.text = element_text(size = 18))
```

![](../Submodule03_files/figure-html/unnamed-chunk-15-1.png)<!-- -->


``` r
ggplot(lsq, aes(x = Weight)) +
  geom_histogram(bins = 15, fill = "lightcoral", color = "black") +
  labs(title = "Weight Distribution", x = "Weight", y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 20),
        axis.text = element_text(size = 18))
```

![](../Submodule03_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

### Step 3: Combine Metadata

We have two different surveys, let's combine them into one so we only have one metadata file to work with. 


``` r
# Load in data from Submodule 2
sample.names <- readRDS("sampleNames.rds")
taxa         <- readRDS("taxa.rds")
reads        <- readRDS("seqtabnochim.rds")
```


``` r
# Merge data
metadata <- merge(ffq, lsq, by = "SampleID", all = TRUE)
# Rename samples to match FastQ files
metadata$SampleID <- gsub("WP_", "WP-", metadata$SampleID)
```

We have a lot more samples in our metadata than we do sequenced, therefore we need to pull out only the samples who have sequenced reads from our metadata.

``` r
# Add in only if they have sequences
dfseqIndo <- data.frame(sampleNames = sample.names,
                       inSeq       = rep('yes', length(sample.names)))
# Only keep those who do
metadata <- merge(dfseqIndo, metadata, by.x = 'sampleNames', by.y = 'SampleID', all.x = TRUE)
colnames(metadata)[1] <- 'SampleID'
rownames(metadata)    <- metadata$SampleID
```



``` r
dim(metadata)
```

```
## [1]   35 1107
```

We now have 35 samples with metadata and sequenced reads.

When merging metadata from two different surveys, it is common to encounter discrepancies due to human error, such as participants providing slightly different responses across surveys for variables like age or weight. While these inconsistencies highlight the importance of careful data handling, for this analysis, we will rely on the FFQ data as it is more standardized and consistent. 


``` r
# Flag discrepancies and missing data
discrepancies <- metadata %>%
  mutate(
    age_discrepancy = ifelse(Age != AGE, TRUE, FALSE),
    weight_discrepancy = ifelse(Weight != WEIGHT, TRUE, FALSE)
  )


# Summarize discrepancies and missing data
summary <- data.frame(
  age_discrepancy = sum(discrepancies$age_discrepancy, na.rm = TRUE),
  weight_discrepancy = sum(discrepancies$weight_discrepancy, na.rm = TRUE)
)

summary
```

```
##   age_discrepancy weight_discrepancy
## 1               9                  4
```

Our data shows some inconsistencies. We will use the FFQ data, as it is a standardized survey. Now, let's examine our data more by creating a table to get an overview of our data by sex and BMI.


``` r
# Order BMI
metadata$BMI_ordinal = factor(metadata$BMI_ordinal, levels = c('Healthy weight', 'Overweight', 'Obese'))
```



``` r
# Group data by Sex
table1 <- metadata %>%
  group_by(Sex) %>%
  summarize(
    Count = n(),
    Mean_Age = mean(AGE, na.rm = TRUE),
    Mean_Weight = mean(WEIGHT, na.rm = TRUE),
    Mean_Height = mean(Height, na.rm = TRUE),
    Mean_BMI = mean(BMI, na.rm = TRUE),
    Mean_Num_Animals = mean(num_animals, na.rm = TRUE),
    Mean_Kcal = mean(DT_KCAL, na.rm = TRUE),
    Mean_Protein = mean(DT_PROT, na.rm = TRUE),
    Mean_Carbohydrate = mean(DT_CARB, na.rm = TRUE),
    Mean_Fat = mean(DT_TFAT, na.rm = TRUE),
    Mean_Alcohol = mean(DT_ALCO, na.rm = TRUE),
    Mean_Sugars = mean(DT_SUG_T, na.rm = TRUE),
    Mean_Fiber = mean(DT_FIBE, na.rm = TRUE)
  )

# View the table
table1
```

```
## # A tibble: 2 × 14
##   Sex    Count Mean_Age Mean_Weight Mean_Height Mean_BMI Mean_Num_Animals
##   <chr>  <int>    <dbl>       <dbl>       <dbl>    <dbl>            <dbl>
## 1 Female    19     35.7        142.        64.7     24.0            0.947
## 2 Male      16     35.2        179.        69.5     26.2            1.25 
## # ℹ 7 more variables: Mean_Kcal <dbl>, Mean_Protein <dbl>,
## #   Mean_Carbohydrate <dbl>, Mean_Fat <dbl>, Mean_Alcohol <dbl>,
## #   Mean_Sugars <dbl>, Mean_Fiber <dbl>
```



``` r
# Group data by BMI
table2 <- metadata %>%
  group_by(BMI_ordinal) %>%
  summarize(
    Count = n(),
    Mean_Age = mean(AGE, na.rm = TRUE),
    Mean_Weight = mean(WEIGHT, na.rm = TRUE),
    Mean_Height = mean(Height, na.rm = TRUE),
    Male_Percent = mean(SEX == 1, na.rm = TRUE) * 100,
    Female_Percent = mean(SEX == 2, na.rm = TRUE) * 100,
    Mean_Num_Animals = mean(num_animals, na.rm = TRUE),
    Mean_Kcal = mean(DT_KCAL, na.rm = TRUE),
    Mean_Protein = mean(DT_PROT, na.rm = TRUE),
    Mean_Carbohydrate = mean(DT_CARB, na.rm = TRUE),
    Mean_Fat = mean(DT_TFAT, na.rm = TRUE),
    Mean_Alcohol = mean(DT_ALCO, na.rm = TRUE),
    Mean_Sugars = mean(DT_SUG_T, na.rm = TRUE),
    Mean_Fiber = mean(DT_FIBE, na.rm = TRUE)
  )

# View the table
table2
```

```
## # A tibble: 3 × 15
##   BMI_ordinal Count Mean_Age Mean_Weight Mean_Height Male_Percent Female_Percent
##   <fct>       <int>    <dbl>       <dbl>       <dbl>        <dbl>          <dbl>
## 1 Healthy we…    21     31.0        146.        67.7         38.1           61.9
## 2 Overweight      9     42          174.        67.1         55.6           44.4
## 3 Obese           5     42.4        189.        63.2         60             40  
## # ℹ 8 more variables: Mean_Num_Animals <dbl>, Mean_Kcal <dbl>,
## #   Mean_Protein <dbl>, Mean_Carbohydrate <dbl>, Mean_Fat <dbl>,
## #   Mean_Alcohol <dbl>, Mean_Sugars <dbl>, Mean_Fiber <dbl>
```

Let's now categorize dietary intake so we can eaisly view differences.


``` r
# Create metadata variables with levels based on dietary standards
metadata <- metadata %>%
  mutate(
    # Calories (kcal/day, general range: 1500-3000 for adults)
    Kcal_Level = case_when(
      DT_KCAL < 1500 ~ "Low",
      DT_KCAL >= 1500 & DT_KCAL <= 3000 ~ "Moderate",
      DT_KCAL > 3000 ~ "High",
      TRUE ~ NA_character_
    ),
    
    # Protein (g/day, general range: 46-56 g for adults)
    Protein_Level = case_when(
      DT_PROT < 46 ~ "Low",
      DT_PROT >= 46 & DT_PROT <= 56 ~ "Moderate",
      DT_PROT > 56 ~ "High",
      TRUE ~ NA_character_
    ),
    
    # Carbohydrates (g/day, general range: 225-325 g based on 2000 kcal diet)
    Carb_Level = case_when(
      DT_CARB < 225 ~ "Low",
      DT_CARB >= 225 & DT_CARB <= 325 ~ "Moderate",
      DT_CARB > 325 ~ "High",
      TRUE ~ NA_character_
    ),
    
    # Total Fat (g/day, general range: 44-77 g based on 2000 kcal diet)
    Fat_Level = case_when(
      DT_TFAT < 44 ~ "Low",
      DT_TFAT >= 44 & DT_TFAT <= 77 ~ "Moderate",
      DT_TFAT > 77 ~ "High",
      TRUE ~ NA_character_
    ),
    
    # Dietary Fiber (g/day, general recommendation: 25-38 g)
    Fiber_Level = case_when(
      DT_FIBE < 25 ~ "Low",
      DT_FIBE >= 25 & DT_FIBE <= 38 ~ "Moderate",
      DT_FIBE > 38 ~ "High",
      TRUE ~ NA_character_
    ),
    
    
    # Categorize alcohol intake
    Alcohol_Level = case_when(
      DT_ALCO == 0 ~ "No Alcohol",
      DT_ALCO > 0 & DT_ALCO <= 14 ~ "Low",
      DT_ALCO > 14 & DT_ALCO <= 30 ~ "Moderate",
      DT_ALCO > 30 ~ "High",
      TRUE ~ NA_character_
    ),
    
    # Categorize sugar intake
    Sugar_Level = case_when(
      DT_SUG_T <= 25 ~ "Low",
      DT_SUG_T > 25 & DT_SUG_T <= 50 ~ "Moderate",
      DT_SUG_T > 50 ~ "High",
      TRUE ~ NA_character_
    )
  )


# Reorder Fiber_Level
metadata$Fiber_Level <- factor(metadata$Fiber_Level, levels = c("Low", "Moderate", "High"))

# Reorder Alcohol_Level
metadata$Alcohol_Level <- factor(metadata$Alcohol_Level, levels = c("No Alcohol", "Low", "Moderate", "High"))

# Reorder Sugar_Level
metadata$Sugar_Level <- factor(metadata$Sugar_Level, levels = c("Low", "Moderate", "High"))


# Reorder Kcal_Level
metadata$Kcal_Level <- factor(metadata$Kcal_Level, levels = c("Low", "Moderate", "High"))

# Reorder Protein_Level
metadata$Protein_Level <- factor(metadata$Protein_Level, levels = c("Low", "Moderate", "High"))

# Reorder Carbohydrate_Level
metadata$Carb_Level <- factor(metadata$Carb_Level, levels = c("Low", "Moderate", "High"))

# Reorder Fat_Level
metadata$Fat_Level <- factor(metadata$Fat_Level, levels = c("Low", "Moderate", "High"))

# Treat the number of animals as discrete not continuous
metadata$num_animals = factor(metadata$num_animals, levels = c('0', '1', '2', '3', '4'))
```


Let's view the distribution of our `Protein_Level`.


``` r
ggplot(metadata, aes(x = Protein_Level, fill = Protein_Level)) +
  geom_bar() +
  theme_classic() +
  labs(
    title = "Distribution of Protein Levels",
    x = "Protein Level",
    y = "Count"
  ) +
  scale_fill_manual(values = c('Low' = '#f7b801',  
                                'Moderate' = '#f35b04', 
                                'High' = '#52006A')) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
```

![](../Submodule03_files/figure-html/unnamed-chunk-26-1.png)<!-- -->

Let's view protein intake by sex.


``` r
ggplot(metadata, aes(x = Protein_Level, fill = Protein_Level)) +
  geom_bar() +
  theme_classic() +
  labs(
    title = "Distribution of Protein Levels",
    x = "Protein Level",
    y = "Count"
  ) +
  scale_fill_manual(values = c('Low' = '#f7b801',  
                                'Moderate' = '#f35b04', 
                                'High' = '#52006A')) +
  facet_wrap(~SEX, labeller = as_labeller(c("1" = "Male", "2" = "Female"))) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
```

![](../Submodule03_files/figure-html/unnamed-chunk-27-1.png)<!-- -->

And now by BMI.


``` r
ggplot(metadata, aes(x = Protein_Level, fill = Protein_Level)) +
  geom_bar() +
  theme_classic() +
  labs(
    title = "Distribution of Protein Levels",
    x = "Protein Level",
    y = "Count"
  ) +
  scale_fill_manual(values = c('Low' = '#f7b801',  
                                'Moderate' = '#f35b04', 
                                'High' = '#52006A')) +
  facet_wrap(~BMI_ordinal) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
```

![](../Submodule03_files/figure-html/unnamed-chunk-28-1.png)<!-- -->


### Step 4: Assign Taxonomy to our Sequencing Data

We will be using the package `phyloseq` to analyze our microbiomes.


``` r
# Install and load libraries
if (!requireNamespace("phyloseq", quietly = TRUE)) {
    install.packages("phyloseq")
}
if (!requireNamespace("ggpubr", quietly = TRUE)) {
    install.packages("ggpubr")
}



library(phyloseq)
library(ggpubr)
```


``` r
# Assigning taxonomy
ps <- phyloseq(otu_table(reads, taxa_are_rows = F),
               tax_table(taxa),
               sample_data(metadata))

ps
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1446 taxa and 35 samples ]
## sample_data() Sample Data:       [ 35 samples by 1114 sample variables ]
## tax_table()   Taxonomy Table:    [ 1446 taxa by 7 taxonomic ranks ]
```



``` r
write.table(otu_table(reads, taxa_are_rows = F), 
            "ASV_table.txt", row.names = TRUE, quote = FALSE, sep = '\t')
```


### Step 5: Data Analysis and Visualization 

#### Rarefaction Curves

Rarefaction curves are essential tools in microbiome analysis used to assess the richness of microbial communities and the adequacy of sequencing depth. These curves plot the number of observed species against the number of sequencing reads, providing insight into species diversity within a sample.

The shape of a rarefaction curve indicates whether sequencing captured the community’s full diversity: a plateau suggests that most species have been identified, while a steep incline indicates that additional sequencing may uncover more species. Rarefaction also enables comparisons across samples with varying sequencing depths, ensuring that observed differences in diversity are not artifacts of unequal sampling effort.

Let's start by looking at the differences between sex and protein intake.



``` r
rare <- rarecurve(reads,step = 1000, tidy = TRUE)
rare <- inner_join(rare, metadata, by=c("Site"="SampleID"))

ggplot(rare, aes(x=Sample, y=Species, col=Protein_Level, linetype=Sex, group=Site))+
  geom_line(linewidth = 2, alpha = 0.8)+
  scale_color_manual(values = c('Low' = '#f7b801',  
                                'Moderate' = '#f35b04', 
                                'High' = '#52006A')) +
  theme_pubr(legend="right") +
  labs(x="Sample Size")
```

![](../Submodule03_files/figure-html/unnamed-chunk-32-1.png)<!-- -->

Now let's plot BMI and protein intake.



``` r
ggplot(rare, aes(x=Sample, y=Species, col=Protein_Level, linetype=BMI_ordinal, group=Site))+
  geom_line(linewidth = 2, alpha = 0.8) +
  scale_color_manual(values = c('Low' = '#f7b801',  
                                'Moderate' = '#f35b04', 
                                'High' = '#52006A')) +
  theme_pubr(legend="right") +
  labs(x="Sample Size")
```

![](../Submodule03_files/figure-html/unnamed-chunk-33-1.png)<!-- -->


#### Diversity in Microbiome

Alpha and beta diversity are key concepts in microbiome studies used to understand the composition of microbial communities.

  - Alpha diversity refers to the diversity within a single sample, capturing the richness (number of species) and evenness (distribution of species). It provides insights into how diverse a community is at a specific site or condition. Common metrics include species richness, Simpson’s index, and the Shannon index.

 
  - Beta diversity, on the other hand, measures the variation in microbial communities between samples. It highlights differences in composition and can reveal patterns of microbial distribution across environments or conditions. Common beta diversity metrics include Bray-Curtis dissimilarity and UniFrac distances.


##### Alpha Diversity

Let's look at the alpha diversity in our samples with the Shannon index. The Shannon index measures both the richness and evenness of species within a sample. A higher Shannon index indicates greater diversity, while a value of 0 means only one species is present in the community.



``` r
shannon <- diversity(reads, index = "shannon") %>%
  as_tibble(rownames = "SampleID")%>%
  inner_join(., metadata, by="SampleID")
```


``` r
ggplot(shannon, aes(x=Sex, y=value, col=Sex, fill=Sex))+
  geom_boxplot(outlier.shape = NA, alpha = 0.5)+
  geom_jitter(width = 0.2) +
  theme_pubr(legend="right") +
  labs(x="",y="Shannon Weaver Index") +
  stat_compare_means(aes(x=Sex, y=value, col=Sex), hide.ns = F, method="t.test",label = "p.signif",
                     label.x = 1.5, label.y= 4.8, show.legend = F) +
  theme(legend.position="none")
```

![](../Submodule03_files/figure-html/unnamed-chunk-35-1.png)<!-- -->

Let's look at BMI within sex.


``` r
# Instal and load ggstance for better visualization
if (!requireNamespace("ggstance", quietly = TRUE)) {
    install.packages("ggstance")
}


library(ggstance)
```



``` r
ggplot(shannon, aes(x=BMI_ordinal, y=value, col=BMI_ordinal, fill=BMI_ordinal))+
  geom_boxplot(alpha=0.5)+  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.5), alpha=0.5) +     
  scale_color_manual(values = c('Healthy weight' = '#f7b801', 'Overweight' = '#f18701', 
                                'Obese' = '#780116')) +
  scale_fill_manual(values = c('Healthy weight' = '#f7b801', 'Overweight' = '#f18701', 
                                'Obese' = '#780116')) +
  theme_pubr(legend="right") +
  labs(x="",y="Shannon Weaver Index") +
  theme(legend.position="none")
```

![](../Submodule03_files/figure-html/unnamed-chunk-37-1.png)<!-- -->


What about alpha diversity by protein intake?


``` r
ggplot(shannon, aes(x=Protein_Level, y=value, col=Protein_Level, fill=Protein_Level))+
  geom_boxplot(alpha=0.5)+  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5), alpha=0.5) +     
  scale_color_manual(values = c('Low' = '#f7b801', 'Moderate' = '#f18701', 
                                'High' = '#780116')) +
  scale_fill_manual(values = c('Low' = '#f7b801', 'Moderate' = '#f18701', 
                                'High' = '#780116')) +
  theme_pubr(legend="right") +
  labs(x="",y="Shannon Weaver Index") +
  theme(legend.position="none")
```

![](../Submodule03_files/figure-html/unnamed-chunk-38-1.png)<!-- -->


Let's now look at the combination of sex or BMI by protein intake. 


``` r
ggplot(shannon, aes(x=Sex, y=value, col=Protein_Level, fill=Protein_Level))+
  geom_boxplot(alpha=0.5)+ 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5), alpha=0.5) +     
  facet_grid(col=vars(Sex), scales = "free_x")+
  scale_color_manual(values = c('Low' = '#f7b801', 'Moderate' = '#f18701', 
                                'High' = '#780116')) +
  scale_fill_manual(values = c('Low' = '#f7b801', 'Moderate' = '#f18701', 
                                'High' = '#780116')) +
  theme_pubr(legend="right") +
  labs(x="",y="Shannon Weaver Index") 
```

![](../Submodule03_files/figure-html/unnamed-chunk-39-1.png)<!-- -->



``` r
ggplot(shannon, aes(x=BMI_ordinal, y=value, col=Protein_Level, fill=Protein_Level))+
  geom_boxplot(alpha=0.5)+   
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5), alpha=0.5) + 
  facet_grid(col=vars(BMI_ordinal), scales = "free_x")+
  scale_color_manual(values = c('Low' = '#f7b801', 'Moderate' = '#f18701', 
                                'High' = '#780116')) +
  scale_fill_manual(values = c('Low' = '#f7b801', 'Moderate' = '#f18701', 
                                'High' = '#780116')) +
  theme_pubr(legend="right") +
  labs(x="",y="Shannon Weaver Index") 
```

![](../Submodule03_files/figure-html/unnamed-chunk-40-1.png)<!-- -->

Next, let's look at the alpha diversity with the Inverse Simpson Index. The Inverse Simpson Index reflects both species richness and evenness within a community. It accounts for the number of species present and their relative abundances, emphasizing the dominance of common species. A higher Inverse Simpson value indicates greater diversity, as it signifies more species evenly distributed in the community. This index is particularly useful for understanding community composition and detecting dominance patterns in microbiome studies. We will focus on our experimental questions, looking at the combination of variables. 



``` r
invsimp<-diversity(reads, index = "invsimpson")%>%
  as_tibble(rownames = "SampleID")%>%
  inner_join(., metadata, by="SampleID")
```


``` r
ggplot(invsimp, aes(x=Sex, y=value, col=Protein_Level, fill=Protein_Level))+
  geom_boxplot(alpha=0.5)+ 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5), alpha=0.5) +     
  facet_grid(col=vars(Sex), scales = "free_x")+
  scale_color_manual(values = c('Low' = '#f7b801', 'Moderate' = '#f18701', 
                                'High' = '#780116')) +
  scale_fill_manual(values = c('Low' = '#f7b801', 'Moderate' = '#f18701', 
                                'High' = '#780116')) +
  theme_pubr(legend="right") +
  labs(x="",y="Inverse Simpson Index") 
```

![](../Submodule03_files/figure-html/unnamed-chunk-42-1.png)<!-- -->


``` r
ggplot(invsimp, aes(x=BMI_ordinal, y=value, col=Protein_Level, fill=Protein_Level))+
  geom_boxplot(alpha=0.5)+ 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5), alpha=0.5) +     
  facet_grid(col=vars(BMI_ordinal), scales = "free_x")+
  scale_color_manual(values = c('Low' = '#f7b801', 'Moderate' = '#f18701', 
                                'High' = '#780116')) +
  scale_fill_manual(values = c('Low' = '#f7b801', 'Moderate' = '#f18701', 
                                'High' = '#780116')) +
  theme_pubr(legend="right") +
  labs(x="",y="Inverse Simpson Index") 
```

![](../Submodule03_files/figure-html/unnamed-chunk-43-1.png)<!-- -->


##### Beta Diversity

Non-metric Multidimensional Scaling (NMDS) is an ordination method used to visualize differences in community composition between samples, making it a key tool for analyzing beta diversity. It represents pairwise dissimilarities (e.g., Bray-Curtis distances) in a low-dimensional space, preserving the rank order of distances rather than their exact values.

The result is a visual representation where samples with similar microbial communities appear closer together, while dissimilar communities are farther apart. NMDS is particularly useful for identifying patterns and grouping structures in complex ecological or microbiome datasets. The quality of the ordination is assessed using a stress value, with lower values indicating a better representation of the original data.



``` r
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
```

```
## Run 0 stress 0.19587 
## Run 1 stress 0.2035157 
## Run 2 stress 0.2451144 
## Run 3 stress 0.2019726 
## Run 4 stress 0.1958701 
## ... Procrustes: rmse 6.881275e-05  max resid 0.000269416 
## ... Similar to previous best
## Run 5 stress 0.1957315 
## ... New best solution
## ... Procrustes: rmse 0.005442248  max resid 0.02678304 
## Run 6 stress 0.2547582 
## Run 7 stress 0.2237335 
## Run 8 stress 0.1957316 
## ... Procrustes: rmse 0.0001455843  max resid 0.0005954312 
## ... Similar to previous best
## Run 9 stress 0.2369315 
## Run 10 stress 0.19587 
## ... Procrustes: rmse 0.005434551  max resid 0.02671452 
## Run 11 stress 0.1957316 
## ... Procrustes: rmse 0.0001236024  max resid 0.0005064555 
## ... Similar to previous best
## Run 12 stress 0.2463371 
## Run 13 stress 0.19587 
## ... Procrustes: rmse 0.005433799  max resid 0.02671089 
## Run 14 stress 0.2609654 
## Run 15 stress 0.1957315 
## ... Procrustes: rmse 4.041635e-05  max resid 0.0001226979 
## ... Similar to previous best
## Run 16 stress 0.2397441 
## Run 17 stress 0.2237335 
## Run 18 stress 0.2060787 
## Run 19 stress 0.2074005 
## Run 20 stress 0.2638578 
## *** Best solution repeated 3 times
```

``` r
# Plot beta diversity                                   
plot_ordination(ps.prop, ord.nmds.bray, color="Sex", title="Bray NMDS") + theme_pubr(legend="right")
```

![](../Submodule03_files/figure-html/unnamed-chunk-44-1.png)<!-- -->

Let's also view by BMI and protein intake.


``` r
plot_ordination(ps.prop, ord.nmds.bray, color="BMI_ordinal", title="Bray NMDS") + scale_color_manual(values = c('Healthy weight' = '#f7b801', 'Overweight' = '#f18701', 
                                'Obese' = '#780116')) + theme_pubr(legend="right")
```

![](../Submodule03_files/figure-html/unnamed-chunk-45-1.png)<!-- -->



``` r
plot_ordination(ps.prop, ord.nmds.bray, color="Protein_Level", title="Bray NMDS") +
scale_color_manual(values = c('Low' = '#f7b801', 'Moderate' = '#f18701', 
                                'High' = '#780116')) +
  theme_pubr(legend="right")
```

![](../Submodule03_files/figure-html/unnamed-chunk-46-1.png)<!-- -->

Next, we will explore the distribution of taxa at the phylum level across our samples using three different datasets: the full dataset, the proportion-transformed dataset, and the subset of the top 20 most abundant taxa. Visualizing these datasets will allow us to examine the overall microbial composition, highlight differences in relative abundance, and focus on the most dominant taxa within the community.


``` r
# Creating top 20 most abundant
top20    <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- prune_taxa(top20, ps.prop)

# Plotting by Sex
plot_bar(ps, x="SampleID", fill="Phylum") + facet_wrap(~Sex, scales="free_x") + theme_bw() +  theme(axis.text.x =element_text(angle=90))
```

![](../Submodule03_files/figure-html/unnamed-chunk-47-1.png)<!-- -->

``` r
plot_bar(ps.prop, x="SampleID", fill="Phylum") + facet_wrap(~Sex, scales="free_x") + theme_bw() +  theme(axis.text.x =element_text(angle=90))
```

![](../Submodule03_files/figure-html/unnamed-chunk-47-2.png)<!-- -->

``` r
plot_bar(ps.top20, x="SampleID", fill="Phylum") + facet_wrap(~Sex, scales="free_x") + theme_bw() +  theme(axis.text.x =element_text(angle=90))
```

![](../Submodule03_files/figure-html/unnamed-chunk-47-3.png)<!-- -->



``` r
# Plotting by BMI
plot_bar(ps, x="SampleID", fill="Phylum") + facet_wrap(~BMI_ordinal, scales="free_x") + theme_bw() +  theme(axis.text.x =element_text(angle=90))
```

![](../Submodule03_files/figure-html/unnamed-chunk-48-1.png)<!-- -->

``` r
plot_bar(ps.prop, x="SampleID", fill="Phylum") + facet_wrap(~BMI_ordinal, scales="free_x") + theme_bw() +  theme(axis.text.x =element_text(angle=90))
```

![](../Submodule03_files/figure-html/unnamed-chunk-48-2.png)<!-- -->

``` r
plot_bar(ps.top20, x="SampleID", fill="Phylum") + facet_wrap(~BMI_ordinal, scales="free_x") + theme_bw() +  theme(axis.text.x =element_text(angle=90))
```

![](../Submodule03_files/figure-html/unnamed-chunk-48-3.png)<!-- -->



``` r
# Plotting by protein intake
plot_bar(ps, x="SampleID", fill="Phylum") + facet_wrap(~Protein_Level, scales="free_x") + theme_bw() +  theme(axis.text.x =element_text(angle=90))
```

![](../Submodule03_files/figure-html/unnamed-chunk-49-1.png)<!-- -->

``` r
plot_bar(ps.prop, x="SampleID", fill="Phylum") + facet_wrap(~Protein_Level, scales="free_x") + theme_bw() +  theme(axis.text.x =element_text(angle=90))
```

![](../Submodule03_files/figure-html/unnamed-chunk-49-2.png)<!-- -->

``` r
plot_bar(ps.top20, x="SampleID", fill="Phylum") + facet_wrap(~Protein_Level, scales="free_x") + theme_bw() +  theme(axis.text.x =element_text(angle=90))
```

![](../Submodule03_files/figure-html/unnamed-chunk-49-3.png)<!-- -->


Now, we will focus on a single phylum of interest, *Firmicutes*, to examine its distribution across our study participants. By extracting this specific phylum, we can explore how its relative abundance varies across key categories such as sex, BMI levels, and protein intake. This detailed investigation will provide insights into potential associations between this microbial group and individual characteristics, helping us better understand the role this phylum may play in the overall microbiome composition.



``` r
# Group data at Phylum level - combines ASVs belonging to same phylum
ps.phylum <- tax_glom(ps.prop, taxrank = "Phylum")

# Convert phyloseq object to long format dataframe for easier manipulation
dat <- psmelt(ps.phylum)

# Find rows where Phylum is "Firmicutes"
idx.keep <- which(dat$Phylum == "Firmicutes")

# Create new dataframe with only Firmicutes data
Firmicutes <- dat[idx.keep,]

# Count number of Firmicutes observations to verify subsetting worked
table(Firmicutes$Phylum)
```

```
## 
## Firmicutes 
##         35
```



``` r
# Plot
ggplot(Firmicutes)+
  geom_col(aes( x=Sample,y=Abundance, fill=Phylum), position="dodge")+
  facet_grid(col=vars(Protein_Level), space="free_x",scales="free")+
  scale_fill_manual(values =  c('#780116')) +
  labs(title = "Abundance by Protein Level") +
  theme_bw() +
  theme(axis.text.x =element_text(angle=90)
  ) 
```

![](../Submodule03_files/figure-html/unnamed-chunk-51-1.png)<!-- -->


We can now use a Kruskal-Wallis test to determine if there are different amounts of *Firmicutes* by protein intake.




``` r
# Create comparison groups 
my_comparisons <- list( c("Low", "Moderate"), c("High", "Low"), c("Moderate", "High") )

# Plot
ggplot(Firmicutes, aes(x=Protein_Level,y=Abundance, fill = Protein_Level, color = Protein_Level))+
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  scale_color_manual(values = c('Low' = '#f7b801', 'Moderate' = '#f18701', 
                                'High' = '#780116')) +
  scale_fill_manual(values = c('Low' = '#f7b801', 'Moderate' = '#f18701', 
                                'High' = '#780116')) +
  stat_compare_means(comparisons = my_comparisons, label.y = c(1, 1.2, 1.1))+
  stat_compare_means(label.y = 1.3) +
  labs(title = "Firmicutes Abundance by Protein Intake", x = "Protein Intake", y = "Firmicutes Abundance") +
  theme_classic() +
  theme(legend.position = "none")
```

![](../Submodule03_files/figure-html/unnamed-chunk-52-1.png)<!-- -->


It seems there is no significant difference in *Firmicutes* abundance between levels of protein intake. 


#### Discussion of Results
This study sought to investigate how dietary patterns, particularly protein intake, and external factors like sex and BMI influence gut microbiome diversity and composition. By examining alpha and beta diversity across these variables, we aimed to uncover meaningful relationships that could shed light on the dynamic interplay between diet, physiology, and microbial communities.


Despite our hypotheses, the findings did not reveal significant differences in microbiome diversity or composition between levels of protein intake, sex, or BMI. The lack of observable trends suggests that these factors may not independently drive major shifts in the microbiome, or that their influence is subtle and requires larger sample sizes to overcome individual variation. It is also possible that interactions among dietary components or unmeasured lifestyle factors play a more significant role, emphasizing the complexity of the diet-microbiome relationship.



Although no clear patterns emerged, these findings highlight the need to further investigate the multifactorial influences on the microbiome to better understand its variability. Future research should consider broader dietary patterns, such as fiber and fat intake, as well as longitudinal approaches to capture temporal shifts in microbial communities. Investigating additional environmental and genetic factors could also provide deeper insights into the interactions shaping the gut microbiome.


Ultimately, while our exploration did not confirm our hypotheses, it underscores the complexity of microbiome research and the intricate interactions between diet, physiology, and lifestyle. This module, designed as a training exercise using a subset of data from a larger study, highlights the importance of refining research questions and methodologies in this evolving field. We hope this serves as a foundation for learners to better understand how to analyze and interpret microbiome data, providing context for future studies and reserch.


## Conclusion
In this module, we successfully explored taxonomic classification and diversity analysis using microbiome sequencing data. By integrating dietary and lifestyle metadata with sequencing results, we evaluated the diversity within and between microbial communities. We analyzed alpha diversity using indices such as Shannon and Inverse Simpson to understand richness and evenness, and beta diversity through Non-metric Multidimensional Scaling to assess community composition differences. Additionally, we visualized microbial composition at various taxonomic levels and identified significant associations, such as the relationship between protein intake and *Firmicutes* abundance. These analyses provide a foundation for understanding how diet and other factors influence microbiome composition and diversity.


$~$








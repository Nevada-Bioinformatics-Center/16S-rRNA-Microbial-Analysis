![course card](images/anchor.jpg)

# Microbial analysis of 16S rRNA sequencing data

## **Contents**

- [Microbial analysis of 16S rRNA sequencing data](#microbial-analysis-of-16s-rrna-sequencing-data)
  - [**Contents**](#contents)
  - [**Overview**](#overview)
  - [**Background**](#background)
  - [**Before Starting**](#before-starting)
  - [**Getting Started**](#getting-started)
  - [**Software Requirements**](#software-requirements)
  - [**Architecture Design**](#architecture-design)
  - [**Data**](#data)
  - [**Funding**](#funding)
  - [**License for Data**](#license-for-data)

## **Overview**

This cloud-based learning module introduces the principles of 16S rRNA sequencing and its applications in microbial community analysis. 16S rRNA gene sequencing is a commonly used method to study the diversity and composition of microbial communities by focusing on a specific region of the ribosomal RNA (rRNA) gene that is present in all bacteria. This sequencing technique generates a vast amount of data. Understanding how to process and analyze this data through a series of computational steps is critical in studies related to the human gut microbiome, among others.

The overall structure of the module is explained below:

- [Submodule 01](Submodule01.ipynb) aims to provide an understanding of the principles of 16S rRNA sequencing and its applications in microbial community profiling.
- [Submodule 02](Submodule02.ipynb) teaches how to preprocess raw 16S rRNA sequencing data and perform quality control checks to ensure reliable results.
- [Submodule 03](Submodule03.ipynb) covers taxonomic classification of 16S rRNA sequences, evaluating rarefaction curves, and analyzing alpha and beta diversity in microbial communities.

## **Background**

The study of microbial communities and their dynamics is essential for understanding ecosystems, human health, and environmental processes. Traditional methods of identifying and characterizing microorganisms, such as culturing, often fail to capture the full diversity of microbes present in a given sample. This is where molecular techniques, particularly 16S rRNA sequencing, have revolutionized microbial analysis.

16S rRNA sequencing focuses on a component of the small subunit of the ribosome found in all bacteria and archaea. This gene is highly conserved among different species, but also contains hypervariable regions that differ between taxa, making it an ideal target for identifying and classifying microbes at various taxonomic levels - from phylum to species.

The process involves amplifying the gene from environmental DNA samples, followed by sequencing and bioinformatic analysis. The resulting data provides a detailed profile of the microbial community present in a sample, including the relative abundance of various species. This has profound applications in fields such as environmental microbiology, clinical diagnostics, gut microbiome research, and biotechnology, where understanding microbial composition can reveal important insights into health, disease, and ecological function.

The objective of this module is to enhance understanding of 16S rRNA sequencing principles by analyzing microbial community composition and diversity. This involves bioinformatics, including quality control, sequence processing, taxonomic assignment, and diversity analysis.

## **Before Starting**

This learning module was designed to be executed on Amazon Web Services (AWS) cloud platform, with the aim of requiring nothing but the files within this GitHub repository. Please, follow the steps highlighed **here** to create a new user account in AWS, open AWS Sagemaker Studio and launch JupyterLab.

## **Getting Started**

 Each learning sub-modules is organized in Jupyter notebooks with step-by-step hands-on practice. This module uses the coding languge R to install necessary tools, obtain data, perform analyses, and visualize the results.

- Step 1: Open the Terminal to clone the GitHub repository with the modules and the data. Use the Git command:

```
git clone https://github.com/Nevada-Bioinformatics-Center/NIGMS-UNR-2024_CloudModule
```

<p align="center">
 <img width= "1280" height= "720" src="images/step11.png">
</p>

- Step 2: After running the command, a new folder will be visible in the left pane labeled "NIGMS-UNR-202..". Select it.

<p align="center">
 <img width= "1280" height= "720" src="images/step12.png">
</p>

- Step 3: Select `Submodule01.ipynb` file in the left pane to open it in the center pane.

<p align="center">
 <img width= "1280" height= "720" src="images/step13.png">
</p>

- Step 4: Once done with `Submodule01.ipynb`, proceed by closing it in the center pane, then opening `Submodule02.ipynb`. Likewise, perform the same operation when proceeding to `Submodule03.ipynb`

- Step 5: When reading through the tutorial, users will encounter code blocks that can be executed in JupyterLab by either clicking on the cell and pressing `Shift + Enter` or by pressing the `Play` button in the top toolbar.

<p align="center">
 <img width= "1280" height= "720" src="images/step14.png">
</p>

- Step 6: When finished with the tutorial, please shutdown the cloud resources. Select `File` in the top bar then `Shutdown`.

<p align="center">
 <img width= "1280" height= "720" src="images/step15.png">
</p>

## **Software Requirements**

The following software are required to run this cloud module.

R Version:
```
R version 4.4.1 (2024-06-14)
```

Required R Packages:

- dada2 (version 1.32.0)
- ggplot2 (version 3.5.1)
- devtools (version 2.4.5)
- patchwork (version 1.3.0.9000) | `devtools::install_github("thomasp85/patchwork")`
- phyloseq (version 1.48.0)
- vegan (version 2.6-8)
- lubridate (version 1.9.3)
- dplyr (version 1.1.4)
- tidyverse (version 2.0.0)
- remotes (version 2.5.0)
- ggvegan (version 0.1.999) | `remotes::install_github("gavinsimpson/ggvegan")`
- ggpubr (version 0.6.0)
- rstatix (version 0.7.2)
- stringr (version 1.5.1)

**Note:** The above listed R packages and the R environment has already been created via a Docker container image `054037098114.dkr.ecr.us-east-1.amazonaws.com/sagemaker/r-notebook` and is deployed via AWS SageMaker when running this cloud module.

## **Architecture Design**


## **Data**

 Gut microbiome data from the [**WOLFPACK Study**](https://freselab.org/wolfpack/) (Wide Open Local Fecal sample collection comparing Pharmaceutical intake, ACtivity, and dietary intaKe) analyzed in this cloud module is made possible thanks to the permission of [_Dr. Steven Frese_](https://www.unr.edu/nutrition/directory/frese-steven), Assistant Professor and the Principal Investigator of the study. Conducted by the University of Nevada, Reno , the [**WOLFPACK**](https://freselab.org/wolfpack/) study investigates the effects of diet, health, and lifestyle impact on the gut microbiome of adults living in Northern Nevada.

 Sequencing for this data was conducted by the [_Idaho State University Molecular Research Core Facility_](https://www.isu.edu/mrcf/), RRID:SCR_012598.

## **Funding**

The project described was fully supported by a grant from the _National Institute of General Medical Sciences_ of the _National Institutes of Health_ under Award Number _GM103440_

## **License for Data**

**Training-Only Data Use Disclaimer:**

 The data, materials, and content provided in this training module are intended solely for educational and instructional purposes. They are not to be used for novel research, publication, or commercial applications.

By accessing and using this module, you agree to the following conditions:

- The data and materials are provided as-is with no guarantees regarding their accuracy, completeness, or suitability for research purposes.
- You may use the materials only within the context of this training and may not repurpose them for independent research or external dissemination.
- Any unauthorized use, including but not limited to publication, presentation, or redistribution of these materials outside of the training environment, is strictly prohibited.

For any questions regarding permitted use or to request an exception, please contact _nbc@unr.edu_ and _sfrese@unr.edu_.
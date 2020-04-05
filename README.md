# Table of Contents
   * [MultiGWAS](#multigwas)
   * [Requirements](#requirements)
   * [Installation](#installation)
      * [General steps to install multiGWAS on a linux system](#general-steps-to-install-multigwas-on-a-linux-system)
      * [Installation of multiGWAS on Ubuntu linux 19.XX](#installation-of-multigwas-on-ubuntu-linux-19xx)
      * [Installation of multiGWAS on Ubuntu linux 18.04](#installation-of-multigwas-on-ubuntu-linux-1804)
   * [Running the examples](#running-the-examples)
   * [General usage](#general-usage)
   * [Genomic data formats](#genomic-data-formats)
      * [Genotype](#genotype)
      * [Phenotype](#phenotype)
   * [Configuration file](#configuration-file)
   * [Considerations](#consideratios)
      * [Implementation](#implementation)
      * [Number of SNPs in Manhattan and QQ plots](#number-of-snps-in-manhattan-and-qq-plots)
      * [Correction for multiple testing](#correction-for-multiple-testing)
   

# MultiGWAS
MultiGWAS is a tool for GWAS analysis that integrates the results of multiple GWAS software tools. Currently, MultiGWAS uses four tools: the GWASpoly R library, SHEsis,  PLINK, and TASSEL. The first two for analysis in polyploids organisms and the other two for diploids. 

This repository includes:
 - The sources of MultiGWAS (R code and bash scripts)
 - The four tools: GWASpoly R library, SHEsis and PLINK (1.9 and 2.0 version) binary programs, and TASSEL Java files 
 - The preinstalled R libraries specific for MultiGWAS.
 - The instructions about how to install MultiGWAS in a Linux system (i.e., general and specific installation).
 - The instructions to run the examples included in MultiGWAS.
 - The description of the input file formats and configuration files.

As an alternative to this minimal installation, there are two other different options to install and test MultiGWAS. The first one is a full installation ([multiGWAS-full](https://github.com/agrosavia-bioinfo/multiGWAS-full)), including a preinstalled Java runtime and the R 3.61 sources to be compiled. And second, there is a ready-to-use Linux virtual machine created with VirtualBox ([MultiGWAS-vm](https://github.com/agrosavia-bioinfo/multiGWAS-vm)).


# Requirements:
  - x86_64 GNU/Linux (Tested with Ubuntu Linux 18.04 LTS, Kernel 5.3.0-40-generic)
  - Java runtime JRE (Tested with OpenJDK-11-JRE)
  - R 6.1  (Tested with R 6.1.3)
  - Pandoc general markup converter (Tested with Pandoc 1.19)
  - Git version control system
  
# Installation
## General steps to install multiGWAS on a linux system
```
# 1. Install pandoc markup converter: 
# 2. Install git tool to clone github repository:
# 3. Install Java runtime:
# 4. Install R 3.61 (included in default Ubuntu 19.XX repositories):
# 5. Clone the multiGWAS repository
git clone https://github.com/agrosavia/multiGWAS-min.git
# 6. Execute the multiGWAS installer:
. INSTALL-MIN.SH
```

## Installation of multiGWAS on Ubuntu linux 19.XX: 
```
# 1. Install pandoc markup converter: 
sudo apt install pandoc

# 2. Install git tool to clone github repository:
sudo apt install git

# 3. Install the default Java runtime:
sudo apt install default-jre

# 4. Install R 3.61 (included in default Ubuntu 19.XX repositories):
sudo apt install r-base

# 5. Clone the multiGWAS repository
git clone https://github.com/agrosavia/multiGWAS-min.git

# 6. Execute the multiGWAS installer:
. INSTALL-MIN.SH
```

## Installation of multiGWAS on Ubuntu linux 18.04:
```
# 1. Install Pandoc markup converter: 
sudo apt install pandoc

# 2. Install git tool to clone github repository:
sudo apt install git

# 3. Install the default Java runtime:
sudo apt install default-jre

# 4. Steps to install R 3.61 (not included in the Ubuntu 18.XX default repositories):
# 4.1. Install the packages necessary to add a new repository over HTTPS
sudo apt install apt-transport-https software-properties-common

# 4.2. Enable the CRAN repository and add the CRAN GPG key to your system:
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'

# 4.3. Update the packages list and install the R package
sudo apt update
sudo apt install r-base

# 5. Clone the multiGWAS repository
git clone https://github.com/agrosavia/multiGWAS-min.git

# 6. Execute the multiGWAS installer and :
. INSTALL-MIN.SH
```
# Running the examples
The multiGWAS folder contains an "example" subfolder with genomic data ("example-genotype.tbl" and "example-phenotype.tbl" files) and two ready-to-use configuration files ("full.config" and "naive.config" files for **Naive** and **Full** GWAS analysis, respectively). Both genotype and phenotype come from the Solanaceae Coordinated Agricultural Project (SolCAP) potato diversity panel described in the paper. 

To run a Full GWAS analysis with the multiGWAS tool, follow the next steps:
 - Open a Linux terminal
 - Change to the "multiGWAS" folder
 - Change to the "example" subfolder
 - Execute the **multiGWAS tool** using as argument the configuration file "full.config":
  ```
        multiGWAS full.config 
  ```
 - An output folder will be created named as "out-XXXX" where XXXX is the prefix of the configuration filename. For the above example, an "out-full" subfolder will be created with the following files and subfolders:
    - A "multiGWAS-report.html" file in HTML format with the full report from the multiGWAS results.
    - A "report" subfolder with the resulting tables and graphics included in the previous report file.
    - An "out" subfolder that contains temporary files created by multiGWAS and the other GWAS tools.
    - A "logs" subfolder that contains the log outputs from the different tools.
    
# General usage
  - Create a new folder (e.g. "test" folder).
  - Copy the phenotype and genotype files to new folder (see data formats below)
  - Create a configuration file (e.g. test.config) and copy it to the new folder (see configuration file below)
  - Open a Linux terminal
  - Change to the new folder
  - Execute the multiGWAS tool using as argument the configuration file
```
      multiGWAS test.config
```
  - Results will be saved in the "out-test" folder
 
# Genomic data formats
## Genotype:
The genotype file is formatted as a table separated by commas with the names of the variables (columns) in the first line. The first three columns corresponding to the marker name, chromosome number, and position in the chromosome. The next columns contain the marker data for each individual in the population codified in the ”ACGT” format (e.g., AATT, CCGG, AAAT, GGCG). 

For example, for a genotype with five individuals and four markers (SNPs), the file looks like this:
```
 Marker,Chrom,Position,ACBrador,AdirondackBlue,AllBlue,AlpineRusset,Alturas
 c2_41437,0,805179,AAAG,AAGG,AAGG,AAAA,AAAA
 c2_24258,0,1252430,AAGG,AGGG,GGGG,GGGG,GGGG
 c2_21332,0,3499519,TTCC,TTCC,TTCC,CCCC,CCCC
 c2_21320,0,3810687,TTTT,TTTT,TTTT,TTTT,TTTT  
```
## Phenotype 
The phenotype file is formatted as a table separated by commas with the names of the variables (columns) in the first line and with two columns. The first one containing the name (or ID) of the individual, and the second the trait value. 

For example, for a quantitative phenotype with five individuals the file looks like this:
```
Name,tuber_shape
ACBrador,3.59
AdirondackBlue,4.07
AllBlue,4.73
AlpineRusset,4.85
Alturas,4.46
```
# Configuration file
The configuration file is a text file with a list of parameter names and their values separated by a colon (":"). The first line contains the text "default:" and the next lines contain the parameters included in the configuration file, which are: 
|Parameter name| Description|
|--------------|------------|
|genotypeFile      | The genotype filename (full path)|
|phenotypeFile     | The phenotype filename (full path)|
|significanceLevel | The genome-wide significance threshold α (commonly 0.01 or 0.05)|
|correctionMethod  | The method for multiple testing correction (”Bonferroni” or ”FDR”)|
|gwasModel         | The type of GWAS analysis (”Naive” or ”Full”)|
|filtering         | TRUE or FALSE whether to use quality control (QC) filters or not (see below) |
|MAF               | Minor allele frequency QC filter |
|MIND              | Individual missing rate QC filter |
|GENO              | SNP missing rate QC filter |
|HWE               | Hardy-Weinberg threshold QC filter|

For example, this is the contents for the above "full.config" configuration file:
```
default:
 genotypeFile         : "example-genotype.tbl"
 phenotypeFile        : "example-phenotype.tbl"
 significanceLevel    : 0.05
 correctionMethod     : "Bonferroni"     # FDR, Bonferroni
 gwasModel            : "Naive"
 filtering            : TRUE
 MAF                  : 0.01   #0.001    # Default 0.01
 MIND                 : 0.1    #0.5      # Default 0.1
 GENO                 : 0.1    #0.5      # Default 0.1
 HWE                  : 1e-10  #1e-50    # Default 1e-10
```

# Considerations
## Implementation
Most of the code uses the R language. However, some scripts that calling the GWAS tools are writing in bash. The version of the four tools are GWASpoly 1.3 (R library), SHEsis 1.0 (binary program), PLINK 1.9 and 2.0 (binary programs), and TASSEL 5.0 (Java packages). PLINK 1.9 is used for GWAS analysis (association between SNPs and quantitative traits), and PLINK 2.0 is used to account for cryptic relatedness (estimating kinship coefficientes).

## Number of SNPs in Manhattan and QQ plots
The Manhattan and QQ plots for the different GWAS tools show a different number of markers (SNPs). Two reasons explain this pattern. First, the GWASpoly software uses four models for the marker effect (i.e., additive, general, simplex dominance, and duplex dominance). Therefore,  the plots show the SNPs four times, one for each model. Second, MultiGWAS is using scores instead of raw p-values, and scores are the -log10(p) results. So, when p-values are high, the scores have a negative value, and because the y-axes in the plot start in zero, they are not shown.

## Correction for multiple testing
MultiGWAS is using two methods for correction for multiple testing, the Bonferroni correction and adjusting the False Discovery Rate (FDR). MultiGWAS calculates the Bonferroni correction using the number of non-missing genotypes (NMISS) included in the analysis instead of the whole genotypes. Only SHEsis, PLINK, and TASSEL give the NMISS number. In contrast, GWASpoly does no show the NMISS number, but it uses it internally to calculate the corrections.



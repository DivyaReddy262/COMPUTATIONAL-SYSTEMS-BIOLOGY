# Assignment 2: SNP Association and eQTL Analysis

## Overview
This assignment focuses on statistical and bioinformatics analysis of Single Nucleotide Polymorphisms (SNPs) in the context of Alzheimer’s disease. It includes association testing using chi-square analysis and functional interpretation through eQTL analysis.

## Objectives
- Perform SNP association analysis using Chi-square test
- Identify significant SNPs from GWAS Catalog
- Analyze SNP positions relative to genes (upstream, downstream, within)
- Study chromosome distribution of SNPs
- Perform eQTL analysis using BRAINEAC database

## Dataset
- GWAS Catalog SNP data (Alzheimer’s disease)
- Significant SNPs (p-value < 5E-8)
- eQTL dataset (BRAINEAC)

## Tools & Technologies
Python
Jupyter Notebook
Libraries: Pandas, NumPyMatplotlib
Excel (data preprocessing)

## Methodology
# Part 1: SNP Association Analysis
- Constructed case-control contingency table
- Calculated expected counts under independence
- Computed Chi-square statistic
- Evaluated statistical significance (p < 0.05)

# Part 2: eQTL Analysis
1. Extracted significant SNPs from GWAS Catalog
2. Mapped SNPs to nearest genes
3. Categorized SNP locations:
- Upstream
- Downstream
- Within genes
4. Generated visualizations:
- Bar plot of SNP locations
- Chromosome distribution plot
- Histogram of SNP-gene distances
5. Identified eQTLs:
- Used BRAINEAC database
- Filtered significant eQTLs (p < 0.05)
- Recorded associated genes and tissues

## Output
- SNP location distribution plots
- Chromosome-wise SNP distribution
- Distance histogram
- eQTL table (SNP, gene, p-value, tissue)

## Files
- assignment-2-analysis.ipynb --> Python implementation
- .xlsx --> SNP and eQTL datasets
- assignment-2.docx --> report

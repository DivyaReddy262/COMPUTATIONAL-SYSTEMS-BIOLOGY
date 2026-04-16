# Assignment 3: Bayesian Networks and WGCNA Analysis

## Overview
This assignment focuses on probabilistic modeling using Bayesian Networks and gene co-expression analysis using Weighted Gene Co-expression Network Analysis (WGCNA). The goal is to understand gene dependencies and identify biologically meaningful gene modules associated with traits.

## Objectives
## Part 1: Bayesian Networks
- Compute joint probabilities using chain rule
- Analyze dependency structures between variables
- Calculate conditional probabilities
- Determine minimum parameters for network representation

## Part 2: WGCNA Analysis
- Perform gene co-expression network analysis
- Identify gene modules using clustering
- Merge similar modules
- Associate gene modules with phenotypic traits

## Dataset
Gene expression data (HW03_expression.csv)
Trait data (HW03_Traits.csv)

## Tools & Technologies
R Programming
WGCNA package
Hierarchical clustering
Statistical correlation analysis

## Methodology
## WGCNA Workflow
1. Data preprocessing and filtering
2. Selection of soft-thresholding power
3. Construction of adjacency matrix
4. TOM-based dissimilarity calculation
5. Hierarchical clustering of genes
6. Dynamic tree cutting to detect modules
7. Merging of similar modules
8. Calculation of module eigengenes
9. Correlation of modules with traits

## Results
- Soft Threshold Selection
Chosen power ensures scale-free topology
Balance between R² fit and mean connectivity
- Gene Clustering
Hierarchical clustering using TOM dissimilarity
Modules identified using dynamic tree cutting
- Module Detection & Merging
Multiple gene modules detected
Similar modules merged based on correlation threshold
- Module–Trait Relationships
Strong associations observed between specific modules and traits
Heatmap visualization highlights significant correlations

## Analysis 
- Gene modules represent groups of co-expressed genes
- Modules significantly correlated with traits may indicate biological relevance
- WGCNA effectively reduces high-dimensional gene data into interpretable structures

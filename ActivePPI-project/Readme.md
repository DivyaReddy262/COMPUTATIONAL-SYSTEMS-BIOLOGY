# ActivePPI-Based Pathway Activity Analysis
## Overview:
This project implements a pathway activity analysis framework inspired by the ActivePPI method, which integrates protein-protein interaction (PPI) networks with expression data to identify biologically significant pathways.
The objective is to evaluate pathway activity across different biological conditions and compare the performance of multiple distance metrics.

## Objectives:
- Apply ActivePPI-based pathway activity analysis
- Compare different distance metrics:
Canberra
Cosine
Euclidean
Manhattan
Minkowski
- Evaluate pathway significance using p-values
- Analyze results across two datasets:
- Breast Cancer (BRCA)
- SARS-CoV-2

## Methodology
The approach is based on the ActivePPI framework:
- Protein interactions are modeled using a network-based approach
- Pathways are evaluated based on activity scores
- Multiple distance metrics are used to compare pathway behavior
- Statistical significance is assessed using p-values
- Results are visualized using plots and heatmaps

## Results & Observations
Distance Metric Comparison
- Euclidean & Manhattan:
Lower variance
More stable results
- Canberra:
Higher variability
Captures extreme differences
- Cosine & Minkowski:
Moderate sensitivity

Biological Insights
- Pathway activity is consistent across datasets
- SARS-CoV-2 dataset shows:
Lower variability
- More uniform biological response
BRCA dataset shows:
Higher variation across pathways

## Visualizations
- Distribution of p-values
Comparison across distance metrics
- Pathway Activity
log10(p-value) representation
Significance threshold highlighted
- Heatmaps
Top pathways across metrics
Comparative visualization

## Note on Implementation
- This project is based on the ActivePPI framework described in the original research paper.
- The analysis reproduces and extends the methodology using custom scripts for visualization and comparison of distance metrics.
- Original work:
https://github.com/zpliulab/ActivePPI

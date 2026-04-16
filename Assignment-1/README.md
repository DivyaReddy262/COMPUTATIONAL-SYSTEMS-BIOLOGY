# Assignment 1: Protein-Protein Interaction Network Analysis

## Overview
This assignment focuses on analyzing the Human Protein-Protein Interaction (PPI) network using data from the BioGRID database. The objective is to explore network properties and identify important proteins using graph theory concepts.

## Objectives
- Compute node degree and identify top hub proteins
- Analyze degree distribution to determine if the network is scale-free
- Compute shortest path lengths between node pairs
- Calculate centrality measures (betweenness and closeness)
- Identify key proteins based on network importance

## Dataset
Source: BioGRID database
Data: Human Protein-Protein Interaction network

## Tools & Technologies
Python
Jupyter Notebook
Libraries: NetworkX, Pandas, Matplotlib

## Methodology
1. Imported PPI network data from BioGRID
2. Constructed graph using NetworkX
3. Calculated degree for each node
4. Plotted degree distribution
5. Computed shortest paths between nodes
6. Calculated:
      - Betweenness centrality
      - Closeness centrality
## Analysis
- The network shows characteristics of a scale-free network, where few nodes (hub proteins) have very high connectivity
- High centrality nodes play a crucial role in maintaining network structure
- Hub proteins such as EGFR and PARK2 are biologically significant

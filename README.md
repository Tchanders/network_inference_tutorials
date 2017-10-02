# README

## Introduction
This directory contains scripts and test data for experimenting with PIDC network inference, introduced in [Gene Regulatory Network Inference from Single-Cell Data Using Multivariate Information Measures](http://www.cell.com/cell-systems/fulltext/S2405-4712(17)30386-1) (Cell Systems, 2017).

Data were simulated using GeneNetWeaver, as described in _Methods_. Our analyses on experimental datasets can be reproduced by downloading the datasets from the relevant references, putting them into the same format as our simulated datasets, and pointing the __Infer network__ script to them.

### Directories

goldstandards: True network structures, from GeneNetWeaver, to compare with inferred networks.

networks: Networks inferred from the various datasets, using the various algorithms. (NB The R package MINET can be used to reproduce our results for the ARACNE, CLR and MINET algorithms.)

notebooks: Scripts for inferring networks and speed benchmarking.

simulated_datasets: Single-cell datasets, simulated using GeneNetWeaver, from the goldstandards.


## Requirements

### Julia version
Julia 0.6

[Download](https://julialang.org/downloads/)

[Documentation](https://docs.julialang.org/en/stable/)

### Packages
InformationMeasures (version 0.1.1 or higher)

`Pkg.add("InformationMeasures")`

PyPlot

`Pkg.add("PyPlot")`

LightGraphs

`Pkg.add("LightGraphs")`

GraphPlot

`Pkg.add("GraphPlot")`

NetworkInference

`Pkg.add("NetworkInference")`

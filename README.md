# Covid 19 Drug Discovery Repo
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3986979.svg)](https://doi.org/10.5281/zenodo.3986979)
[![GitHub issues](https://img.shields.io/github/issues/AliSajid/Covid19)](https://github.com/AliSajid/Covid19/issues)
[![GitHub forks](https://img.shields.io/github/forks/AliSajid/Covid19)](https://github.com/AliSajid/Covid19/network)
[![GitHub stars](https://img.shields.io/github/stars/AliSajid/Covid19)](https://github.com/AliSajid/Covid19Covid19/stargazers)

![GitHub release (latest SemVer including pre-releases)](https://img.shields.io/github/v/release/AliSajid/Covid19?include_prereleases&label=latest-release)
![GitHub release (latest SemVer)](https://img.shields.io/github/v/release/AliSajid/Covid19?label=latest-stable)
[![GitHub license](https://img.shields.io/github/license/AliSajid/Covid19)](https://github.com/AliSajid/Covid19/blob/main/LICENSE)

![GitHub language count](https://img.shields.io/github/languages/count/AliSajid/Covid19)
![GitHub top language](https://img.shields.io/github/languages/top/AliSajid/Covid19)
![Lines of code](https://img.shields.io/tokei/lines/github/AliSajid/Covid19)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/AliSajid/Covid19)
![GitHub repo size](https://img.shields.io/github/repo-size/AliSajid/Covid19)

## Introduction

This repository contains the code for the COVID-19 Drug Discovery Paper by O'Donovan et al published in *journal*. This repository contains all the code andd supplementary files needed to complete the analysis from start to finish. It also contains the data files that we used to generate particular results from each version for quicker verification.

## Quickstart

This analysis runs best on *nix or macOS operating systems. Run the following commands to run the analysis including the data download.

```bash
$ git clone https://github.com/AliSajid/Covid19.git
$ cd Covid19
$ make
```

## System Requirements

This analysis requires R (> 4.0.0)

## Repository Structure

```bash
.
├── ...                  	  # Boilerplate files and primary scripts
├── data                 	  # Contains all the data downloaded and generated from iLINCS
│ ├── A549-10uM-24h      	  #
│ ├── A549-10uM-6h       	  #
│ ├── HA1E-10uM-24h      	  #
│ ├── HT29-10uM-24h      	  # Filtered drug and group signatures and connected perturbagens
│ ├── MCF7-10uM-24h      	  # for all 8 cell line combinations
│ ├── PC3-10uM-24h       	  #
│ ├── VCAP-10uM-24h      	  #
│ ├── VCAP-10uM-6h       	  #
│ ├── disease            	  # Filtered disease signatures and connected perturbagens
│ └── signatures         	  # Downloaded signatures for diseases and drugs and drug groups
├── figures              	  # Contains all the generated figures
| ├── ...                	  # All primary figures in PNG, JPG and PDF formats
│ ├── densityplots       	  # Density plots for all group signatures signifying normal distribution
│ ├── histograms         	  # Histograms for all group signatures signifying normal distribution
│ └── qqplots            	  # Quartile-Quartile plots for all group signatures signifying normal distribution
├── maps                 	  # Drug signature to drug name maps for all 8 cell line combinations
├── raw                  	  # Files that are imported from outside and shouldn't be deleted
| ├── ...                	  # General files
│ ├── annotation         	  # Annotation and gene counts for the GSE56192 dataset
│ └── drug_signature_list	  # List of signatures for each drug in our list
├── renv                 	  # Renv Dependency Management directory
└── results              	  # Final Results folder
```

## Script Descriptions

### 0-setup.R

This file is used to set up the directory structure and other initialization settings.

### 1-explore_drug_sig_lists.R

This script takes in the drug-specific signatures and creates a file which combines all of that and outputs the groups with maximum representation in the selected drugs

### 2-*.R

The two scripts download the drug and the group data and connected perturbagens

### 3-*

These scripts process the influenza, mers and sars datasets and generaate their signatures

### 4-generate_disease_data.R

This script downloads the perturbagnes connected to the disease signatures, including the SARS-CoV-2 signatures

### 5-*.R

These scripts take the lists of generated perturbagens and create intersections between them at given thresholds

### 6-*.R

These scripts perform the final steps of the analysis including outputting summarized datasets and threshold tables

### f-*.R

These scripts generate all the individual figures

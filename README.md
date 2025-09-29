# Mapping-Structural-Diversity-Using-GEDI
This repository contains the Matlab code used for mapping the structural diversity of Central African and Western US forests using GEDI. 

Please read and reference (cite) the following scientific paper when using this repository:

Fabian D. Schneider*, Morgan Dean, Elsa M. Ordway, Moses B. Libalah, & Antonio A. Ferraz. Mapping the structural diversity of Central African and Western US forests using GEDI. In Review at Remote Sensing of Environment.  
*fabian.schneider@bio.au.dk; Section for Ecoinformatics & Biodiversity, Department of Biology, Aarhus University, Ny Munkegade 114, DK-8000 Aarhus, Denmark

## Usage Information
This code has been run and tested with Matlab version 2024a.

We include example data to test the code in the 'data' folder. Please place this folder together with the 'geotiff-out' folder in the main working directory in Matlab, from which you are running the code. This is necessary, so that relative paths to 'data/GEDI_traits_ranges.mat' will work.

**mapGEDI_Diversity.m**: this is the main script to run the code. This script will call other functions and load the data.
**getFRicFEveFDiv_PDFadapt.m**: this is the main function to calculate structural diversity in multivariate feature space. This function will call akde.m to calculate multivariate probability densities.
**akde.m**: this is a Matlab function by Zdravko Botev to calculate kernel density in high dimensional feature space. See [Mathworks](https://se.mathworks.com/matlabcentral/fileexchange/58312-kernel-density-estimator-for-high-dimensions) for more details.
**scaleZeroOne_cols_absolute.m** and **scaleNoOutliers_cols.m**: these are helper functions to scale data from 0 to 1 in different ways.

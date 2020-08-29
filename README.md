# Phylogenomic analysis of all available COVID-19 genomes and its influence on mortality

<div align='center'>
  <img src='https://github.com/MountainMan12/GISAID_phylo/blob/master/Images/Tree.png'>
</div>

## Package information
### Python 3.7.4
- pandas 1.0.5
- biopython 1.77

### CODE IMPLEMENTATION

1) Data_filt.py : The python script performs the data cleaning and writes the cleaned data to a new file. This code gets rid of NAs and unknowns from the clinical dataset. The dataset originally contained 4592 samples which after running the python script were 3620. 
  - Input: The script takes input the GISAID clinical data containing all samples, along with the complete GISAID genomic data for 4592 samples.
  - Output: Two files containing filtered samples.

2) For OR calculations see [ORCaD](https://github.com/MountainMan12/GISAID_phylo/tree/master/ORCaD)


### CONFLICTS

Please contact pawan12394@gmail.com in case of any issues


### LICENSE 

The project was licensed under the Apache 2.0 license

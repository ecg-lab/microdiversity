# Examination of microbial diversity patterns

This method has now been published in PNAS!

A null model for microbial diversification  
Timothy J. Straub and Olga Zhaxybayeva  
PNAS July 3,2017 vol. 114 no. 27 E5414-E5423  
10.1073/pnas.1619993114  
http://www.pnas.org/content/114/27/E5414.abstract


## Pre-requisities
* Python 2.7, including Biopython and Numpy
* R (https://www.r-project.org/)
* MCL (https://micans.org/mcl/)
* Mothur (https://mothur.org/)
* MUSCLE (http://drive5.com/muscle/)
* Seq-Gen (http://tree.bio.ed.ac.uk/software/seqgen/)

## General steps in the analysis pipeline
1. Run the pipeline on empirical data. See `pipeline/pipeline.sh` for details.  
a. Fast alignment-free similarity comparisons of all vs all protein sequences using afree  
b. Clustering on similarity using MCL to generate gene families  
c. Filtering on gene families to generate near single-copy core  
d. Furthest-neighbor clustering within gene families using Mothur  
2. Simulate neutral gene families based on empirical data. See `sim/simulations.sh` for details.  
a. Simulate phylogenies using a Moran model  
b. Generate DNA sequences based on these phylogenies using Seq-Gen  
c. Clustering with Mothur  
3. Comparisons between simulated and empirical data using R. See `statistics.R` for details.  
a. Load in both empirical and simulated clustering results  
b. Compare empirical gene family clustering patterns to simulated null model  
c. Investigate gene families that are significantly different from null model  

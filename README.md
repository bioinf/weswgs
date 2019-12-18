# Determinants of WES and WGS sequencing coverage 

This repository contains all code used to analyze the data and plot figures in the paper:

*Barbitoff Y.A., Polev D.E., Shcherbakova E.A., Kiselev A.M., Glotov A.S., Serebryakova E.A., Kostareva A.A., Glotov A.S., Glotov O.S., and Predeus A.V. (2019)* Systematic dissection of biases in whole-exome and whole-genome sequencing reveals major determinants of coding sequence coverage. *Scientific Reports*

<img align="center" width="800" height="800" src="https://github.com/bioinf/weswgs/blob/master/img/circa1.png">



## Subfolder contents:

`./coverage_analysis` - all scripts used to make alignment and coverage data manipulations

`./coverage_analysis/multimap/` - scripts to analyze coverage difference upon MQ > 10 filtering
    
`./coverage_analysis/norm_curves/` - scripts to calculate normalized coverage profiles from BEDGRAPH and histogram files generated by `collect_coverage_data.sh`
    
`./coverage_analysis/wie_profiles/` - scripts to make mean WIE profiles for a selection of samples, per-platform

`./Fig_1 - Fig_5` - R scripts and data files used to create figures

`./variant_analysis` - scripts used to analyze variant calling results

`./linear_predictions` - scripts and dataset for running GLM and random forest predictions of normalized coverage

## Additional 

For Fig_3, some larger data files are available through Google Drive:

https://drive.google.com/drive/folders/179iwL44LDdSeII6D_iy_uXYOOBNyyP58?usp=sharing

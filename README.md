# Determinants of WES and WGS sequencing coverage 

This repository contains all code used to analyze the data and plot figures in the paper:

Barbitoff YA, Polev DE, Shcherbakova EA, Kiselev AM, Glotov AS, Serebryakova EA, Kostareva AA, Glotov AS, Glotov OS, Predeus AV [Systematic dissection of biases in whole-exome and whole-genome sequencing reveals major determinants of coding sequence coverage](https://www.nature.com/articles/s41598-020-59026-y), *Scientific Reports* **2020**, 10(1), 1-13.

Preprint is available at [bioRxiv](https://www.biorxiv.org/node/117933.abstract).

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

For Fig\_3, some larger data files are available for download via [Google Drive](https://drive.google.com/drive/folders/179iwL44LDdSeII6D_iy_uXYOOBNyyP58?usp=sharing).

## Contacts

If you have any questions, please contact Yury A Barbitoff (barbitoff at bioinf me) or Alexander V Predeus (predeus at bioinf me).

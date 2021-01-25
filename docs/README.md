# BGCViz
BGCViz is a shiny application (can be run locally, or via [website](https://biopavlohrab.shinyapps.io/BGCViz/)), which uses genomic coordinates of BGC annotations from different sources to visualize interception between them in one sequence (genome). This integrative approach could point to both:
 - "regions of interest", annotated with more that one tool.
 - novel regions, annotated only by one of the methods.

This app is written as a part of **Cambridge Bioinformatics Hackaton 2020** ([link](https://cambiohack.uk)). 

![Biocircos](images/biocircos.png)

Currently app supports annotations from such sources:
1. Antismash (both json and csv)
2. PRISM (both json and csv)
3. SEMPI (only csv)
4. DeepBGC (raw tsv file)
5. RRE-Finder (modified txt file)
6. ARTS (raw csv files)

For more information about which files need to be preparated please see [input files](#Input_files_options) section [TEST](#TEST TEST)
# Where to go from here

# Contributing

# TO-DO list

## Prerequisites
App is written entirely in R. <br>
 **Note:**  the [web version](https://biopavlohrab.shinyapps.io/BGCViz/) is also available, with no installation needed. <br>
Needed for a local run:
- R ( v. 4.0.2 )  [site](https://www.r-project.org)
- Rstudio ( v. 1.3.959 )  [site](https://rstudio.com)
- shiny ( v. 1.5.0 )
- tidyverse ( v. 1.3.0 )
- plyr ( v. 1.8.6 )
- IntervalSurgeon ( v. 1.0 )
- plotly ( v. 4.9.2.1 )
- BioCircos ( v. 0.3.4 )
- ggplot2 ( v. 3.3.2 )



## Installation
All packages can be installed via R console in Rstudio via:
```R
install.packages(c("BioCircos", "ggplot2", "plotly", "IntervalSurgeon", "plyr", "tidyverse", "shiny" ))
``` 
If conda package manager is installed (see [Anaconda](https://www.anaconda.com) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)) you can use the following command to install R and Rstudio:

`conda install r-base rstudio`

## Usage
Local run:
1. Open app.R in Rstudio
2. Press "Run app" button on the upper right corner of the script window.

The app should open in a viewer panel in Rstudio, or a default web browser. The plots will be visible as soon as input files are uploaded. 
### Input files
- Antismash and PRISM are comma-separated files, with 4 columns -> "Cluster", "Start", "Stop", "Type", where "Cluster" column holds clusted IDs as numbers. <br> Example: (for data see example_data folder) <br>
![antismash](images/antismash.png) <br>
- deepbgc.tsv file is an standart tsv file output from DeepBGC. <br>
Example: (data in in example_data folder) <br>
![deepbgc](images/deepbgc.png) <br>
- RREFinder.txt file is standart txt output file from RREFinder tool. But as a delimiter double underscore is used "__". So for proper parsing of this file, expected sequence name is: "Sequence-name__Coordinates__Gene-ID",  possible "__Coordinates__product-description__Gene-ID" is reformated with product-decription deleted. <br>
Example: (data example is in example_data folder) <br>
![rrefinder](images/rrefinder.png) <br>
### Controls
The app consist of three input categories:
1. Files upload and chromosome length input (mandatory inputs) <br>
![files](images/upload.png) <br>
2. DeepBGC and antismash data comparison plots controls. Those used for the first two result plots <br>
![comp](images/deepbgc_expl.png) <br>
3. DeepBGC filtering controls. Applied globally + download datasets button is available <br>
![deep_filt](images/deepbgc_filt.png) <br>
### Result plots
1. Barplot with a comparison of DeepBGC and antismash data. Have three categories -> "Annotated only by antismash", "Annotated by antismash and deepbgc", "Annotated only by deepbgc"
![deepbgc](images/barplot1.png)
2. Connected scatterplot with Novelty, Annotation and Skip rates, where:
   - Novelty rate = "# of BGC annotated only by deepbgc"/("# clusters annotated with only by antismash" + "# clusters annotated with antismash and deepbgc"). This rate points to how much clusters are annotated only by DeepBGC.
   - Annotation rate = "# of BGC annotated by antismash and deepbgc"/"total number of antismash annotated BGC". This rate points to how much DeepBGC annotated clusters alongside with antismash. 
   - Skip rate = "# of BGC annotated only by antismash"/"total number of antismash clusters". This rate points of how much clusters DeepBGC missed, assuming, that antismash is a reference annotation <br>

![rates](images/rates.png) <br>
3. "BGC on chromosome" plots. On hover additional data is available 
   - First row is antismash BGC data plotted
   - Rows under have the following structure: Reference_data_with_only_intercepted_results, data_with_which_the_interception_was_performed. The example can be "A_vs_D" and "D" rows. So "A_vs_D" holds antismash clusters (we use this annotation as a reference), which intercepts with the DeepBGC data (the row below). "D" row is all DeepBGC data plotted (not only intercepted. 

![ref](images/reference.png)
4. BioCircos plot. Here the all-vs-all comparison is performed. Hover on links gives a little bit of additional info, as IDs of the cluster of connected annotation tools, names of those tools and types of clusters.
<br>

![biocirc](images/biocircos.png)
5. Barplot, which shows how many times the given cluster is annotated by other tools.
![barplot2](images/barplot2.png)
## Citation
This project is still a work in progress, so there is no official publication. If it was useful you can cite it as: <br>
P. Hrab & B. Ostash 2020: BGCViz, GitHub repository: https://github.com/pavlohrab/BGCViz, doi: 10.13140/RG.2.2.23431.01444

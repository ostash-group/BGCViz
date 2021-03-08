# Installation

The app itself is written in R and is available as a [shiny app server](https://ostash-group.shinyapps.io/BGCViz/) with no istallation required. However local run have some advantages(See Why local run section) over remote one. Considering ease of installation and UI indentity, we are suggetting local installation

## Base packages
All packages can be installed via R console in Rstudio via:
```R
install.packages(c("BioCircos", "ggplot2", "plotly",  "plyr", "tidyverse", "shiny", "DT" ,"rjson", "stringr", "shinyjs" ))
``` 
After, to install GenomicRanges package please run
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GenomicRanges")
```

If conda package manager is installed (see [Anaconda](https://www.anaconda.com) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)) you can use the following command to install R and Rstudio:

`conda install r-base rstudio`

## Packages for input files conversion

To covert input files from the program's native output format to BGCViz-compatible csv one, a set of R scripts in the /scripts folder is available. They depends on `rjson` library to convert antismash and PRISM data,  and `RSQLite` library to covert the SEMPI data. Also common libraries in all three scripts are: `dplyr`, `tidyr`, `stringr`. Therefore one line installation from R console is:
```R
install.packages(c("dplyr", "tidyr", "stringr", "RSQLite", "rjson" ))
```

## Why local run
The local run of an app is identical to the website one, with the same UI in the browser. Howewer, it have several advantages:
- Fast upload of files
- Better response time
- Control over execution. Quick modifications to the code.

Besides speed and flexibility of a local run, ability of upload raw  json output files from PRISM and AntiSMASH is a big plus (more details [here](Input_files_options.md)). For server usage we are suggesting to convert the json files to the csv ones, using provided scripts (more details [here](Input_files_options.md)). 

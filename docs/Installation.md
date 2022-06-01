# Installation

The app itself is written in R and is available as a [shiny app server](https://ostash-group.shinyapps.io/BGCViz/) with no installation required. However local run have some advantages(See Why local run section) over remote one. Considering ease of installation and UI indentity, we are suggetting local installation

## Installation
Install `remotes` package first with R console. (for example via Rstudio)
```R
install.packages("remotes")
``` 
After run the following command in R console:
```R
remotes::install_github("ostash-group/BGCViz")
```
This will install BGCViz as an R package. Running an app is as easy as:
```R
BGCViz::run_app()
```

You computer should have [R](https://cloud.r-project.org) and [Rstudio](https://www.rstudio.com/products/rstudio/download/) installed

If you have troubles with installation, do not hesitate to open an [issue](https://github.com/ostash-group/BGCViz/issues)

## Why local run
The local run of an app is identical to the website one, with the same UI in the browser. Howewer, it have several advantages:
- Much faster upload of files. Especially bigger json files or archives.
- Better response time
- Control over execution. Quick modifications to the code.
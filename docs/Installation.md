# Installation

The app itself is written in R and indented to be run locally, or on own server. Feel free to play with our demo version [shiny app server](https://ostash-group.shinyapps.io/BGCViz/) with no installation required. 

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
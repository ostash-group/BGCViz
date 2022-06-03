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

## Hosting free version on shinyapps.io

The BGCViz can also be hosted as usual shiny app, using any solution you like. This example uses free version of shinyapps.io (suitable for one user). 

First you should clone the github repository:

```bash
git clone https://github.com/ostash-group/BGCViz
```

Then the hosting options are yours to choose. Use `app.R` as application file. Following this simple example you can open the file in Rstudio and run in R console:

```R
>install.packages('rsconnect')

>options(repos = BiocManager::repositories())
```
Then, change the working directory to the one, app.R file is located and run:

```R
rsconnect::deployApp()
```

**In case you are jumped to Markers tab with warnings, return to the R console and respond to the "Filepaths are case-sensitive on deployment server. Do you want to proceed with deployment? [Y/n]:" -> Y**

You can follow the same procedure with GUI. Just click the Rconnect button on right top corner of the script:


![rconnect](/images/rconnect.png)


Then follow the instructions [here](https://shiny.rstudio.com/articles/shinyapps.html)


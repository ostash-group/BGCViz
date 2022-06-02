# Building a Prod-Ready, Robust Shiny Application.
# 
# README: each step of the dev files is optional, and you don't have to 
# fill every dev scripts before getting started. 
# 01_start.R should be filled at start. 
# 02_dev.R should be used to keep track of your development during the project.
# 03_deploy.R should be used once you need to deploy your app.
# 
# 
###################################
#### CURRENT FILE: DEV SCRIPT #####
###################################

# Engineering

## Dependencies ----
## Add one line by package you want to add as dependency
usethis::use_package( "ggplot2" )
usethis::use_package( "shiny" )
usethis::use_package( "shinydashboard" )
usethis::use_package( "shinydashboardPlus" )
usethis::use_package( "shinyjqui" )
usethis::use_package( "shinydisconnect" )
usethis::use_package( "BioCircos" )
usethis::use_package( "plotly" )
usethis::use_package( "plyr" )
usethis::use_package( "dplyr" )
usethis::use_package( "DT" )
usethis::use_package( "rjson" )
usethis::use_package( "stringr" )
usethis::use_package( "shinyjs" )
usethis::use_package( "GenomicRanges" )
usethis::use_pipe()
usethis::use_package( "shinycssloaders" )
usethis::use_package( "sortable" )
usethis::use_package( "tidyr" )
usethis::use_package( "stringi" )
usethis::use_package( "RSQLite" )


## Add modules ----
## Create a module infrastructure in R/
golem::add_module( name = "deep_reference" ) # Name of the module
golem::add_module( name = "deep_reference_2" ) # Name of the module
golem::add_module( name = "barplot_rank" ) # Name of the module
golem::add_module( name = "group_table" ) # Name of the module
golem::add_module( name = "download" ) # Name of the module
golem::add_module( name = "deepbgc_plots" ) # Name of the module
golem::add_module( name = "gecco_plots" ) # Name of the module
golem::add_module( name = "biocircos" ) # Name of the module


## Add helper functions ----
## Creates fct_* and utils_*
golem::add_fct( "biocircos" )
golem::add_fct( "deep_reference" )
golem::add_fct( "group_table" )
golem::add_fct( "validation" )
golem::add_fct( "filtering" )
golem::add_fct( "reading_processing" )
golem::add_fct( "helpers" )
golem::add_fct( "format_transformation" )
golem::add_utils( "deep_reference" )
golem::add_utils( "helpers" )

## External resources
## Creates .js and .css files at inst/app/www
#golem::add_js_file( "script" )
#golem::add_js_handler( "handlers" )
#golem::add_css_file( "custom" )

## Add internal datasets ----
## If you have data in your package
usethis::use_data_raw( name = "anti_data", open = FALSE ) 
usethis::use_data_raw( name = "prism_data", open = FALSE ) 
usethis::use_data_raw( name = "prism_supp_data", open = FALSE )
usethis::use_data_raw( name = "gecco_data", open = FALSE )
usethis::use_data_raw( name = "sempi_data", open = FALSE )


## Tests ----
## Add one line by test you want to create
usethis::use_test("fct_biocircos")

# Documentation

## Vignette ----
usethis::use_vignette("BGCViz")
devtools::build_vignettes()

## Code Coverage----
## Set the code coverage service ("codecov" or "coveralls")
usethis::use_coverage()

# Create a summary readme for the testthat subdirectory
covrpage::covrpage()

## CI ----
## Use this part of the script if you need to set up a CI
## service for your application
## 
## (You'll need GitHub there)
usethis::use_github()

# GitHub Actions
usethis::use_github_action() 
# Chose one of the three
# See https://usethis.r-lib.org/reference/use_github_action.html
usethis::use_github_action_check_release() 
usethis::use_github_action_check_standard() 
usethis::use_github_action_check_full() 
# Add action for PR
usethis::use_github_action_pr_commands()

# Travis CI
usethis::use_travis() 
usethis::use_travis_badge() 

# AppVeyor 
usethis::use_appveyor() 
usethis::use_appveyor_badge()

# Circle CI
usethis::use_circleci()
usethis::use_circleci_badge()

# Jenkins
usethis::use_jenkins()

# GitLab CI
usethis::use_gitlab_ci()

# You're now set! ----
# go to dev/03_deploy.R
rstudioapi::navigateToFile("dev/03_deploy.R")


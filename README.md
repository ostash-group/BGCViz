
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BGCViz

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
![Visitor
count](https://shields-io-visitor-counter.herokuapp.com/badge?page=ostash-group.BGCViz&style=for-the-badge)
![GitHub](https://img.shields.io/github/license/ostash-group/BGCViz?style=for-the-badge)
![GitHub release (latest by date including
pre-releases)](https://img.shields.io/github/v/release/ostash-group/BGCViz?include_prereleases&style=for-the-badge)
![GitHub
issues](https://img.shields.io/github/issues/ostash-group/BGCViz?style=for-the-badge)
![GitHub Repo
stars](https://img.shields.io/github/stars/ostash-group/BGCViz?style=for-the-badge)
![GitHub
contributors](https://img.shields.io/github/contributors/ostash-group/BGCViz?style=for-the-badge)

![DOI](http://img.shields.io/badge/DOI-10.13140%20%2F%20RG.2.2.23431.01444-blue.svg?style=for-the-badge)
<!-- badges: end -->

BGCViz is a shiny application (can be run locally, or via
[website](https://ostash-group.shinyapps.io/BGCViz/)), which uses
genomic coordinates of BGC annotations from different sources to
visualize interception between them in one sequence (genome). This
integrative approach could point to both:

 - “regions of interest”, annotated with more than one tool.

 - novel regions, annotated only by one of the methods.

This app is written as a part of **Cambridge Bioinformatics Hackathon
2020** ([link](https://cambiohack.uk)). 

**All the documentation is available on our
[site](https://ostash-group.github.io/BGCViz)**

![Biocircos](images/biocircos.png)

Currently, the app supports annotations from such sources:

1.  Antismash (both json and csv)

2.  PRISM (both json and csv)

3.  SEMPI (only csv)

4.  DeepBGC (raw tsv file)

5.  RRE-Finder (modified txt file)

6.  ARTS (raw tsv files)

7. RiPPMiner-genome (cured txt file)

# Contributing

There are no contributing guidelines yet. But feel free to resolve any
posted issue on
[repository](https://github.com/ostash-group/BGCViz/issues) or implement
anything from our TO-DO list. 

If you have any questions, suggestions, or bugs with the BGCViz please
let us know via Issues section of the repo on
[GitHub](https://github.com/ostash-group/BGCViz/issues)

# TO-DO list

1.  Write all group information in one GenBank file (for know separate
    files are generated. More info [here](Additional_analysis.md))

2.  Add BigFAM information parsing

3.  Add generation of json files for an input to `--sideload` flag on
    antismash v. 6.0.  Then the data can be analyzed in details using
    antismash interface.

4.  Provide parsing of SQLite file from SEMPI (raw output). This is
    possible locally, but the web version crashes.

# Citation

This project is still a work in progress, so there is no official
publication. If it was useful you can cite it as:

P. Hrab & B. Ostash 2020: BGCViz, GitHub repository:
<https://github.com/ostash-group/BGCViz>, doi:
10.13140/RG.2.2.23431.01444

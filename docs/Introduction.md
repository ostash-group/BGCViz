# BGCViz

BGCViz is a shiny application (can be run locally, or via [website](https://ostash-group.shinyapps.io/BGCViz/)), which uses genomic coordinates of BGC annotations from different sources to visualize interception between them in one sequence (genome). This integrative approach could point to both:

 - "regions of interest", annotated with more than one tool.

 - novel regions, annotated only by one of the methods.

  

This app is written as a part of **Cambridge Bioinformatics Hackathon 2020** ([link](https://cambiohack.uk)). 

  

![Biocircos](images/biocircos.png)

  

Currently, the app supports annotations from such sources:

1. Antismash (both json and csv)

2. PRISM (both json and csv)

3. SEMPI ( all results in zip archive and csv)

4. DeepBGC (raw tsv file)

5. RRE-Finder (modified txt file)

6. ARTS (all results in zip archive and csv)

  

For more information about which files need to be preparated please see [input files](Input_files_options.md) section

# What's next?

1. See our [installation](Installation.md) guide if you decided to run the app locally.

2. You can use our [guide](Quick_start.md) for a quick start.
  

# Contributing

There are no contributing guidelines yet. But feel free to resolve any posted issue on [repository](https://github.com/ostash-group/BGCViz/issues) or implement anything from our TO-DO list. 

If you have any questions, suggestions, or bugs with the BGCViz please let us know via Issues section of the repo on [GitHub](https://github.com/ostash-group/BGCViz/issues)

# TO-DO list

1. Write all group information within one GenBank file (for now separate files are generated. More info [here](Additional_analysis.md))

2. Add BigFAM information parsing

3. Add generation of json files for an input to `--sideload` flag on antismash v. 6.0.  Then the data can be analyzed in details using antismash interface.

4. Provide parsing of SQLite file from SEMPI (raw output). This is possible locally, but the web version crashes.

  

# Citation

This project is still a work in progress, so there is no official publication. If it was useful you can cite it as: 

P. Hrab & B. Ostash 2020: BGCViz, GitHub repository: https://github.com/ostash-group/BGCViz, doi: 10.13140/RG.2.2.23431.01444

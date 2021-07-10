# Input files
We tried to make integration of various services as seemless as possible. However, due to computation limitation of shiny server we sugest to transform at least antismash output json file to the csv dataframe. More information is below.
All the example data is in the `/example_data` directory.

# Antismash
There are two options for antismash data upload:
1. Raw output json file.  You should download all your antismash results from the server, if the run was remote, and in the root results directory find the json file. Then be sure that box, which states "My AntiSMASH data is a dataframe, not json results file from antismash" is UNTICKED before actual upload. 
****We highly suggest to use this option only with local installation****

2. The csv format file. This file can be easily generated from the json file, using antismash.R script in the /scripts directory. Before be sure that all necessary packages are installed (more [here](Installation#packages-for-input-files-conversion)). Then just run the following command from the command line:
 ```R
 Rscript antismash.R <path-to-the-json>
 ```
Script should generate the the csv file in the /scripts directory.
The structure of the csv file is the following:

![anti_data](/images/antismash_data.png)

File itself  can be created manually from the antismash results.

**Note: If the region contain several components we  are expecting the double underscore delimiter between them "__". Alternatively they will not be splited for renaming and will maintened as a one type cluster. See more on renaming [here](BGCViz_renaming_and_coloring_options.md)**

**This field supports the csv file upload. It can be any csv file, as long as the format is the satisfied. Therefore you can upload any results in this field, in place of antismash, but the result plot will label these as "Antismash"**

# PRISM
Similarly to the AntiSMASH, there are two options for the PRISM input:
1. Use the json report file. It can be downloaded after the PRISM run. Because of the small  size of json output file, it can be directly uploaded to the BGCViz. Be sure the checkbox "My PRISM data is a dataframe, not json results file"  is UNTICKED.  **We are expecting that the run was made on version >= 4.4.5 of the PRISM.  If earlier version was used, please provide the csv instead of json one.**
2. The csv format file. Identical to one, described in the [Antismash](#antismash) section above. Can also be created manually. The command for generation from json to csv format is the following:
```R
 Rscript prism.R <path-to-the-json>
```


**If json file was provided, the additional option of using Regulatory and Resistance genes, identified by PRISM can be used from the "Data manipulation option"  menu option in BGCViz. Those genes will be treated as an additional input, and chromosome will be named "PRISM_SUPPORT"**

**This field supports the csv file upload. It can be any csv file, as long as the format is the satisfied. Therefore you can upload any results in this field, in place of PRISM, but the result plot will label these as "PRISM"**

# SEMPI
SEMPI data can be used only as a csv file in the format, described in the [Antismash](#antismash) section above.  Can be created manually.
After SEMPI run please download all the results. 
The db file, script can transform to the .csv file is the   `/genome_browser/main/Tracks.db` one, from the SEMPI results folder.

```R
 Rscript sempi.R <path-to-the-db-file>
```

**This field supports the csv file upload. It can be any csv file, as long as the format is the satisfied. Therefore you can upload any results in this field, in place of SEMPI, but the result plot will label these as "SEMPI"**
# DeepBGC
We are expecting the default DeepBGC tsv output file. Please see the example in /example_data folder

# RRE-Finder
We the expecting the modified default txt results file from RRE-Finder, which is in the format as below (in case of exploratory mode):

![rre_data](/images/rre_data.png)

Or like this (in case of precise mode):

![rre_data_precise.png](/images/rre_data_precise.png)

(These files are created from genbank input to RRE-Finder (which is an expected input by default) )

The expected list of modifications:
1.  The first column ("Gene_name")  is divided in 3 parts by double underscore as a delimiter ("__".)  
2.  First part is a contig name. (Here "S136_genome"). 
3.  Second part is the most valuable one - coordinates.
4.  Third part is a locus name, or other sequence identifier.

**Please be sure that any other data from Gene_name column is deleted. The data itself is divided by double underscores (should be made manually) into three clean parts. Please delete all other data in Gene_name column (product names or gene names), if any**

# ARTS
There are two uploads for the ARTS data. And only two files are provided the data will be included in the analysis.  Firstly, download the ARTS results (all results). 
1. The first expected tsv file is knownhits.tsv (holding the hits to existing resisntance models) in the `<ARTS-results-dir>/tables/knownhits.tsv`
2. The second input is duptable.tsv (holding the duplicated housekeeping genes data) in the  `<ARTS-results-dir>/tables/duptable.tsv`

# GECCO
We are expecting the default GECCO tsv output file. Please see the example in /example_data folder
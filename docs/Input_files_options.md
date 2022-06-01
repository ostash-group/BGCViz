# Input files
We tried to make integration of various services as seemless as possible. However, this is comes with a cost of rather long uplod times of raw data. Therefore we suggest install an application locally (more [here](Installation.md))
# Antismash
There are two options for antismash data upload:
1. Raw output json file.  You should download all your antismash results from the server, if the run was remote, and in the root results directory find the json file. 

2. The csv format file. This file can be easily generated from the json file, using `BGCViz::antismash_to_csv()` function of a package after it's [installation](Installation.md) as following:

```R
BGCViz::antismash_to_csv(file="antismash.json")
```

Script should generate the the csv file in the current directory
The structure of the csv file is the following:

![anti_data](/images/antismash_data.png)

File itself also can be created manually from the antismash results.

**Note: If the region contain several components we  are expecting the double underscore delimiter between them "__". Alternatively they will not be splited for renaming and will maintened as a one type cluster. See more on renaming [here](BGCViz_renaming_and_coloring_options.md)**

**This field supports the csv file upload. It can be any csv file, as long as the format is the satisfied. Therefore you can upload any results in this field, in place of antismash, but the result plot will label these as "Antismash"**

# PRISM
Similarly to the AntiSMASH, there are two options for the PRISM input:
1. Use the json report file. It can be downloaded after the PRISM run.  **We are expecting that the run was made on version >= 4.4.5 of the PRISM.  If earlier version was used, please provide the csv instead of json one.**
2. The csv format file. Identical to one, described in the [Antismash](#antismash) section above. Can also be created manually. The function to generate csv file from json one is:
```R
BGCViz::prism_to_csv(file="prism.json")
```


**If json file was provided, the additional option of using Regulatory and Resistance genes, identified by PRISM is used. You can hide them from analysis in Global options, under "Prism supplement+ARTS options".Chromosome is named "PRISM-Supp"**

**This field supports the csv file upload. It can be any csv file, as long as the format is the satisfied. Therefore you can upload any results in this field, in place of PRISM, but the result plot will label these as "PRISM"**

# SEMPI
SEMPI data can be uploaded as:
1. Zip archive (whole project), downloaded directly from SEMPI.
2. Csv file in the format, described in the [Antismash](#antismash) section above.  Can be created manually. After SEMPI run please download all the results. 
Then you can use the following function:

```R
BGCViz::sempi_to_csv(project_archive="project.zip")
```
SEMPI project archive can be downloaded from the site with "Project" button:

![sempi_res_export](/images/sempi_res_export.png)

**This field supports the csv file upload. It can be any csv file, as long as the format is the satisfied. Therefore you can upload any results in this field, in place of SEMPI, but the result plot will label these as "SEMPI"**
# DeepBGC
We are expecting the default DeepBGC tsv output file. Please see the example in inst/extdata folder

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
ARTS data can be uploaded as:
1. Zip arhive of ARTS results
2. Csv file, which was extracted from zip zrchive, using the following function:
```R
BGCViz::arts_to_csv(project_archive="arts.zip")
```
Zip archive should be downloaded from ARTS Export tab (Zip all files):

![arts_res_server.png](/images/arts_res_export.png)
# GECCO
We are expecting the default GECCO tsv output file. Please see the example in inst/extdata folder
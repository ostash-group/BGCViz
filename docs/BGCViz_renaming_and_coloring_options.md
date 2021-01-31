# BGCViz options
The user can explicitly change the renaming scheme and coloring options. 

# Renaming 

The renaming schema is stored is a csv file, called rename.csv. This file looks like this:

| Code      | Group |  Group_color | Color | Vacant_colors | Hierarchy |
| ----------- | ----------- |  ----------- |  ----------- |  ----------- |  ----------- |
polyketide|pks|pks|#a6cee3||Antismash
nonribosomal_peptide|nrps|nrps|#1f78b4||PRISM
melanin|melanin|ripp|#b2df8a||SEMPI
ectoine|other|saccharide|#33a02c||DeepBGC
pentangular_polyphenol|other|melanin|#fb9a99||RRE
nrps-independent_siderophore_synthase|other|other|#e31a1c||PRISM-supp
angucycline-type|pks|terpene|#fdbf6f||ARTS
angucycline|pks|alkaloid|#ff7f00||
butyrolactone|ripp|hybrid|#cab2d6||
class_i_lantipeptide|ripp|core|#6a3d9a||
lasso_peptide|ripp|regulatory|#ffff99||
nis_synthase|other|resistance|#b15928||
acyl_amino_acids|other|base|#d4ced6||
aminocoumarin|other|||| 
... | ... ||||

The data used for renaming is stored under "Code" and "Group" columns. The default file is available is the BGCViz directory or in the [Glossary](Glossary.md). The one is free to use any renaming scheme, but the whole csv file should be uploaded in the file input in BGCViz:

![renaming_opt.png](/images/renaming_opt.png)

# Coloring renamed data

The renamed categories in the "Group" can be then colored. The categories itself should be specified in the "Group_color", and "Color" column holds the corresponding color in hexadecimal format. The default colors are from [ColorBrewer palette](https://colorbrewer2.org/#type=qualitative&scheme=Paired&n=12) and can be easily changed. 

The "Vacant_color" column is not used in any analysis.  The logic behind it is to store any of the colors, that is not used in the coloring now(due to lack of groups). The selected colors play together nicely, (example: [ColorBrewer palletes](colorbrewer2.org)), so the one wants to store them if the number of groups increases.  As there are 12 (base color excluded) groups in "Group_color" column, so no vacant colors are available from ColorBrewer palette.  

# Hierarchy

The "Hierarchy" column defines the order of the link coloring in the 'Hierarchy-based' mode. More on Biocircos link coloring is available [here](Logic_of_the_output.md#biocircos-plot)
# Additional analysis
After the analysis is done the one can download uploaded data. This option is particularly valuable for DeepBGC and GECCO data cleaning, whereas for other inputs the data remains unchanged.  Alongside the inputs, "group_by.csv" file is generated, which is the table from "Summarize interception" tab and group.py python script. 

Here we depend on the usage of a clinker software, which is an awesome app. Check out their GitHub [page](https://github.com/gamcil/clinker), or published [paper](
https://doi.org/10.1093/bioinformatics/btab007) for more details.

# Logic
The structure of the group_by.csv looks like this:

 |Antismash|DeepBGC|RRE-Finder|PRISM|SEMPI|PRISM-supp|Group |
| ----------- | ----------- |  ----------- |  ----------- |  ----------- |  ----------- | ----------- |
1|NA|NA|1|1|NA|group_1
2|6|NA|NA|NA|NA|group_2
3|11|1|2,19|NA|NA|group_3
4|NA|NA|3|3|NA|group_4
6|45|NA|NA|NA|11,12|group_5
7|NA|NA|5|NA|NA|group_6
8|66|NA|7|NA|NA|group_7
9|NA|NA|8|NA|NA|group_8
10|70|NA|9|12|36,92,93,94|group_9
11|97|NA|10|18|56,57,110,111|group_10
13|100|NA|11|19|NA|group_11
15|NA|NA|12|20|119,120,121|group_12
19|117|NA|13|23|65,122,123,124,125|group_13
20|119|NA|14|25|126,127|group_14
21|NA|NA|15|NA|NA|group_15
23|NA|NA|16|28|70|group_16
24|131|2|17|NA|NA|group_17
26|157|NA|NA|32|NA|group_18
27|161|NA|18|33|75|group_19
NA|13,15,51,64,30,64||6,4,6|4,7,7,10|3,35,78,80,128|group_20

The logic of this post-analysis is to extract the GenBank files, which corresponds to the clusters,  from the master GenBank file. Then put them in the separate Group folders.

**Tip: It's better to use the feature-rich GenBank file, as from antismash annotation for example. This analysis will not add any information to the actual GenBank records**

# Dependencies
For running the grouping "group.py" script, several should be installed in your environment - `pandas`, `biopython`
Also, clinker software is required. They all can be installed via pip from the command line:

```bash
pip install biopython pandas clinker
```

Pip should be already installed if you have [Python](https://www.python.org) installed in your system
# Step 1. Group the GenBank records
## Inputs

There is one input -> genome sequence, which was used for BCG annotations in GenBank format. 

## Usage
The usage is pretty straightforward - you need to specify only one input - master GenBank file :

```bash
python group.py <location-of-your-gb-file>
```

## Results
The script is working rather slow. The grouping can take up to 20-30 min. The result of the grouping is several folders, which are named as "group_1", "group_2", etc. These folders hold extracted records in GenBank format.

# Step 2. Run clinker

## Input
Inputs are the folder with the .gb files. Therefore you can choose the group you want to visualize

## Usage
For a comprehensive guide on how to use clinker software please refer to their GitHub [page](https://github.com/gamcil/clinker). To interactively visualize the certain group, please run:

```
clinker <Group-folder/*.gb> --plot
```

## Results

Then, after the analysis is done, the default browser will open with the interactive visualization. For example, for the default S.coelicolor data, results for group 3 (grouped by antismash data) can be viewed as follows (while using antismash results GenBank file) :

![clinker](/images/clinker_example.png)

This is a result of running clinker with the following command:

```bash
clinker group_3/*.gb --plot
```
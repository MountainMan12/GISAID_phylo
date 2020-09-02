# ORCaG - **O**dds **R**atio **Ca**lculations for **G**isaid data
## Introduction
This is small shine application for Odds Ration calculation for GISAID dataset. Made as a part of HackBio internship. 

<div align='center'>
  <img src='https://github.com/pavlohrab/GISAID_phylo/blob/master/ORCaG/Images/intro.png'>
</div>

## Dependencies 

This application is ready to use. However if you want to run it locally (not on the [shinyapps](https://biopavlohrab.shinyapps.io/ORCaG/)) you need to manually install
- R v 4.02
- Rstudio v. 1.3.959
- Bioconductor
- R libraries

#### R libraries 
- shiny
- tidyverse
- plyr
- ggplot2
- hrbrthemes
- forcats
- ggtree
- ape
- epiR
- tidytree
- gridExtra
- formattable
- grid
- ggplotify
- knitr
- ggpubr
- ggsignif
- phylocanvas
- getopt
- optparse
- plotrix
- readODS

You can install majority using the following code from inside Rstudio :
```R
install.packages("shiny", "tidyverse", "plyr", "ggplot2", "hrbrthemes", "forcats", "ape", "epiR", "tidytree", "gridExtra", "formattable", "grid", "ggplotify", "knitr", "ggpubr", "ggsignif", "phylocanvas"  )
```
Then with Bioconductor:
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("ggtree", "optparse", "plotrix", "getopt", "readODS"))
```

After all libraries are installed you are ready to go!

## Usage
If you are starting app locally, you need to press run app button in the Rstudio (upper border of the script on the right). Or you can [use web version](https://biopavlohrab.shinyapps.io/ORCaG/) from shinyapps.io
App main screen :
<div align='center'>
  <img src='https://github.com/pavlohrab/GISAID_phylo/blob/master/ORCaG/Images/main_screen.png'>
</div>
Then the usage of the GUI is pretty the same. Inputs should be:
<div align='center'>
  <img src='https://github.com/pavlohrab/GISAID_phylo/blob/master/ORCaG/Images/inputs.png'>
</div>
1.GISAID dataset is a must. We assume you are downloading 14-column dataset with `\t` as separator (tsv file). Names of the columns doesn't matter, but we are using the first one as Virus.name, second - as Accession.ID, and so on. The example is below:
<br>

>"Virus.name, Accession.ID, Collection.date, Location, Host, Additional.location.information,Gender, Patient.age, Patient.status,Passage, Specimen, Additional.host.information, Lineage, Clade"
 <br> 
If your dataset is different please note **App will use columns as declared. So if the Location column is empty you with not be able to get Country data, if Clade is empty - no Clade information, and so on..** . The only required columns for OR calculations are Patient.age and Patient.status - names doen't matter but they should be 8th and 9th column of a dataset. <br>
2. Status_assignment.csv file. This is comma-delimited file, where the first column is `Patient.status` categories from GISAID dataset, and the second - categories you want to rename into. If you want to delete category completely. leave cell in the second column blank , or type "-" (without ""). <br> The example of this file is in the folder. Note: **The app will use only the first two columns, and second column should have only two categories for OR calculations**  <br>
3. Tree is a newick format. This input is optional, only required if you want to generate iTOL annotation files or map Patient.status categories to the nodes <br>

When files are uploaded you can do plenty of things:
<div align='center'>
  <img src='https://github.com/pavlohrab/GISAID_phylo/blob/master/ORCaG/Images/controls.png'>
</div>
1. Change Hospitalized patient status (High risk/Low risk/Default(as uploded)/Remove from dataset) (select Menu) <br>
2. Use or not to use Alive/Live/Symptomatic categories of GISAID dataset (check box) <br>
3. Download the currently used dataset (Button) <br>
4. Generate iTOL annotation files and download them (two distinct buttons) <br>
5. Select which category to use for node colouring in iTOL annotation (select menu) <br>
6. Choose with labes are used for newick tree (for tree representation) (select menu) <br>
7. Age cutoff slider (slider) <br>

The side panel have two links:
<div align='center'>
  <img src='https://github.com/pavlohrab/GISAID_phylo/blob/master/ORCaG/Images/links.png'>
</div>
1. Tree upload windows in iTOL <br>
2. Our annotated tree in the iTOL <br>

**Calculations will begin only when gisaid dataset is uploded and status_assignment.csv is uploded. No loading animation is provided, so please wait up to 20s**

## Ouput
The app is generating several output plots:
1. OR calculations. This plot should be interpreted as every point represents OR values for this age cutoff.
<div align='center'>
  <img src='https://github.com/pavlohrab/GISAID_phylo/blob/master/ORCaG/Images/OR_plot.png'>
</div>
2. Text, which corresponds to the OR values with given age cutoff slider input. The OR value is represented with CI (95%)
<div align='center'>
  <img src='https://github.com/pavlohrab/GISAID_phylo/blob/master/ORCaG/Images/OR_text.png'>
</div>
3. (optional) Tree with mapped patient.status categories to the nodes (please select which node labels to use - virus.name, accession.id)
<div align='center'>
  <img src='https://github.com/pavlohrab/GISAID_phylo/blob/master/ORCaG/Images/ggtree.png'>
</div>
4. Age density plot with two categories. Vertical red line dynamically shows current age cutoff
<div align='center'>
  <img src='https://github.com/pavlohrab/GISAID_phylo/blob/master/ORCaG/Images/age_dens.png'>
</div>
5. Country count barplots, which shows countries with patient counts that are younger than selected age cutoff or older.
<div align='center'>
  <img src='https://github.com/pavlohrab/GISAID_phylo/blob/master/ORCaG/Images/barplots.png'>
</div>
6. Age over categories distribution barplot. On the right computed statistics for those categories
<div align='center'>
  <img src='https://github.com/pavlohrab/GISAID_phylo/blob/master/ORCaG/Images/age_box.png'>
</div>
7. Distribution of age over Clades (with p_value computation among them) and Distribution of age over Clades with patient.status categories division. The significance is calculated between means in each clade for two categories
<div align='center'>
  <img src='https://github.com/pavlohrab/GISAID_phylo/blob/master/ORCaG/Images/clade_box.png'>
</div>
8. Phylocanvas interactive tree. 
<div align='center'>
  <img src='https://github.com/pavlohrab/GISAID_phylo/blob/master/ORCaG/Images/phylocanvas.png'>
</div>

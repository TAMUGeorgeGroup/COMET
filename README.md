# COMET
COMET is a computational framework for inferring EMT trajectories from timecourse single-cell RNA sequencing data. Here, we are providing COMET as an R package that can be installed and run for scoring scRNA sequencing data. 
<img src="https://github.com/TAMUGeorgeGroup/COMET/assets/62211977/16e1d9a1-bbe7-4138-a9d5-4ba2aab8a35d" alt="model-diagram" width="500" height="600">
<h2>Guidelines</h2> 

Download <a href="https://drive.google.com/file/d/1JHspm15geKSlaAvaU9fwvZRnnNVOMqPu/view?usp=sharing">this zip file</a> and place it in a directory of your choice, then set your working directory to that path:
<br>
setwd("~/Desktop/COMET_STAR_protocol/")
<br>
Step I. Load relevant libraries
```
library(dplyr)
library(ggplot2)
library(tidyverse)
library(data.table)
library(tidyr)
library(reshape2)
library(Rmagic)
library(umap)
library(readxl)
library(Seurat)
library(archetypes)
library(phateR)
library(dtw)
library(pracma)
```
Step II. Make sure that you set these three variables to the correct directory names
<ul>
  <li>Tables directory, should store a csv table with three columns, "Sample", "DataPath", and "MetaData" each of these should store the name of the input files (with no space), path to count matrix, and path to metadata respectively.</li>
  <li>Input Data directory, should store the metadata and the count matrix both stored in csv.gz formats, the meta.data should store a column named "Time" and have timepoints stored in numeric values.</li>
  <li>Input Files directory contains the necessary files for the pipeline to run such as list of important EMT genes for scoring, note that we are providing them separately so you can modify these lists directly.</li>
</ul>

```
tables.dir <- "Tables/"
input.data.dir <- "Input_Data/"
input.files.dir <- "Input_Files/"
data.inputs <- read.csv(paste0(tables.dir, "DataTableTest.csv"), sep=",")
```

Step III. Run the pipeline

```
COMET::start.comet(tables.dir, input.data.dir, input.files.dir)
COMET::generate_pipeline_files(data.inputs, tables.dir, input.data.dir, input.files.dir)
COMET::calculate_conf_intervals(data.inputs)
COMET::DTW_calculate(data.inputs,  c(7.33, 8, 10))
fit.all.data(data.inputs, c(7.33, 8, 10)) ->final.result
```

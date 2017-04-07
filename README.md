# ShinyCiona
Shiny app for our *Ciona* cardiopharyngeal specification single cell RNA-seq data.

This shiny app is for local runs. You can also check the web interface at [ShinyCiona](https://ciona.shinyapps.io/shinyciona).

## Installation
First, you will need to install the following packages:

`install.packages('shiny')`

`install.packages('gplots')`

`install.packages('RColorBrewer')`

Follow the instructions to install [Seurat](http://satijalab.org/seurat/install.html).

## Usage
Download the ShinyCiona project from [github](https://github.com/stevexniu/ShinyCiona) and the precomputed [R objects](https://drive.google.com/drive/folders/0B7ZaAVsPgGYfZW5BUHFCRy1JZms?usp=sharing).

Unzip the files and put the ***spatial.Robj*** and ***temporal.Robj*** under the **ShinyCiona-master** folder.

Open either ***ui.R*** or ***server.R*** in Rstudio and click **Run App** button. Now you get the ShinyCiona running on your laptop :tada:!

## Welcome to the repository for the **Interspecific hybridization in insects in times of climate change** book chapter
### by
### Rosa Ana Sánchez-Guillén, Luis Rodrigo Arce-Valdés, Andrea Viviana Ballén-Guapacha, Jesús Ordaz-Morales & Miguel Stand-Pérez

In this repository you will find raw literature search results, raw data and scripts for bioinformatic analyses for this book chapter. To have a quick view of supporting tables including all hybrid and barriers species pairs take a look at `RawData_Summaries.xlsx`.

#### Repository structure

Here is the structure of the directory of the book chapter, several folders of intermediary files have been gitignored. Within each folder `bin` includes the scripts that conform our general pipeline. Scripts in the bins folders are written either in bash or in R.

```
+-RawData_Summaries.xlsx					#This excel file contains summarized tables including all the species pairs and references we used in our chapter.
+---Linux									#Linux directory contains all scripts and results for figures and statistical testing.
|   +---01_Consensus						#Contains scripts to download and make a COI consensus sequence
|   |   +---bin
|   |   +---data
|   |   \---results
|   +---02_Genetic_Distances				#Contains scripts for the first set of figures included in the book chapter.
|   |   +---bin
|   |   \---figures
|   +---03_Additional_Species_Distances		#We needed to estimate genetic distances in some additional species. We did this directory for that.
|   |   +---bin
|   |   +---data
|   |   \---results
|   +---04_Barriers_vs_Distance				#In this directory we make analysis between reproductive isolation and genetic distance.
|   |   +---bin
|   |   +---data
|   |   \---results
|   +---05_Barriers_Calculation				#Now in this directory we plotted results from the previous directory.
|   |   +---bin								
|   |   \---figures
|   +---06_COIs_Meta						#We used this directory to create a list of all COI sequences access numbers we downoladed and used from BOLD and Genbank.
|   |   +---bin
|   |   +---data
|   |   \---results
|   +---07_Consequences						#Directory for consequences of hybridization plots.
|   |   +---bin
|   |   \---figures
|   \---08_Publication_BIAS					#Directory for figures for publication bias.
|       +---bin
|       +---data
|       \---figures
\---Windows									#In the windows directory we saved information on literature seraches, our results and filtering and raw data.
    +---Busquedas							#Our directory containing literature searches.
    |   +---Barriers						#This directory contains our results and data colection on reproductive barriers. 
    |   |   +---Diptera
    |   |   +---Hymenoptera
    |   |   +---Lepidoptera
    |   |   +---Odonata
    |   |   |   +---bin
    |   |   |   +---data
    |   |   |   \---results
    |   |   +---Odonata2
    |   |   \---Orthoptera
    |   +---Consequences					#This directory contains our results and data colection on hybridization consquences.
    |   +---Hybrid							#This directory contains our results and data colection on hybridization reports.
    |   |   +---Diptera
    |   |   +---Hymenoptera
    |   |   +---Lepidoptera
    |   |   +---Odonata
    |   |   +---Odonata2
    |   |   \---Orthoptera
    |   \---R								#This directory contains additional r scripts we used for analyse windows data.
```

#### Need help?
If you are seeing this repository and need to do something similar as in our book chapter feel free to email [Luis Rodrigo](mailto:bio.l.rodrigo.arce@gmail.com). 
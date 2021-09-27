## Overview

The data in this app come from performing differential expression testing on treatment
conditions relative to their controls. PCA is performed on the logFC values from
this differential expression testing. The tab-panel on the left allows the user to
select the data to be used in the analysis. 

- **Datasets:** Whether to use protein coding genes or repetitive elements 
- **Data:** Use the logFC of a gene/RE in PCA analysis
- **Conditions:** i.e. contrasts or treatment vs controls. Allows the user to select or remove individual data points.
- **Experiments:** NCBI BioProject IDs for all experiments.
- **Cell Lines:** Cell Lines
- **Drugs:** Drugs
- **Tissues:** Tissues
- **Diseases:** Diseases
- **Epigenetic Class:** One of the epigenetic classes defined by the [HEDDS database](http://hedds.org/search.jsp)

After the data has been selected from the drop downs the user can select the parameters
to be used when computing PCA. **Remove Variance** allows the user to remove a proportion
of low-variance features prior to PCA. **Scale/Center** data specifies whether or not
data should be scaled/centered propr to PCA. **Use complete cases** indicates whether or
not you would like to use genes/REs detected in ALL features. 

Once the data has been selected and the PCA parameters chosen click the **Run** button. PCA
will be performed on the selected data and the plots will update when the tab is selected. 
Meta-Analysis will also be performed on the contrasts selected and used in the PCA.

## PCA Biplot

The PCA Biplot tab has an interactive PCA biplot. You can select the metadata
category to color the points by as well as which axes to plot by using the drop down
menus. The user can also select points on the plot by circling them with the mouse. 
information about the circled points will appear below the plot. 

## PCA loadings

The PCA Loadings tab shows a loadings plot that indicates which features are most influential in driving the PCs. 

## Scree Plot

the Scree Plot tab show the amount of variance explained by each PC.

## Diff Exp.

The Diff Exp. tab allows the user to select a single contrast and view its
differential expression results. The user can select the FDR cutoff and logFC
cutoff and the plots will dynamically update.

## Metadata

The Metadata tabs allows the user to browse all metadata associated with the 
contrasts.

## Features

The Features tab allows the user to search for the differential expression 
data for a given gene/RE across all samples by typing a gene/RE name into the 
search bar.

## Meta-Vote

The Meta-Vote tab implements a vote-counting meta-analysis for all of the 
selected contrasts. The **logFC cutoff** indicates the logFC above which a gene/RE
is considered significant in any study. The **FDR cutoff** indicates the FDR value
below which a gene/RE is considered significant. The **Proportion of experiments
gene is DE** indicates the proportion of experiments where a gene/RE was detected
as differentially expressed in order to be called "meta-differentially expressed". 
**Use only genes common to all experiments** allows the user to select only genes 
detected in all experiments to be used in the meta-analysis. **Sign consistency**
(in the plot) is the number of studies a gene is up-regulated minus the number of
studies in which a genesis down-regulated. 

The user can hover over the plot to display information about each point. 

## Meta-Pcombine

The Meta-Pcombine tab implements a p-value combination meta-analysis technique.
the user can select any one of the p-value combination techniques to combine p-values:

- Fisher: Uses the sum of the log of the p-values. Very sensitive to very small p-values
therefore a very significant single study can lead to a significant combined p-value.
- Pearson: Similar to Fisher but most sensitive to large p-values; therefore more false
negative are obtained.
- Tippet: Uses the minimum of the p-values. Prone to false positives.
- Wilkinson: Uses the max of the p-values: Prone to false negatives.

The user can also select a function to summarize the logFC values. The **FDR cutoff**
and **logFC cutoff** below the plot filter the combined results based on FDR 
and logFC, respectively.
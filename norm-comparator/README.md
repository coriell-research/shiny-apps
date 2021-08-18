# Norm Comparator App

The goal of this app is to allow the user to visualize different RNA-seq 
filtering and normalization methods in order to determine which method
is most appropriate for their data. 

**I had originally designed this app for personal use so I'm expecting it to be a 
bit buggy**, however, I think it is still generally usable. If you run into 
problems using the app please either raise an issue on GitHub (this repo) or 
send me a message on Slack. 

## How do I use this?

The app requires a `SummarizedExperiment` object as input (see 'How do I create a `SummarizedExperiment` object below).
**The `SummarizedExperiment` object that you upload MUST have an assay named 
'counts' containing the raw count data as well as a column in the colData named
'group'.** Other assays and metadata columns are allowed but these two things
must be present. 

To use the app:

1. Sign onto the VPN and navigate your web browser to the app url.
2. Upload the `SummarizedExperiment` object by clicking the "Browse" button under 
"Import SummarizedExperiment"
3. Load the object by clicking the "Load Dataset" button. You should see the 
"select sample" drop-down on the filtering tab update with the name of the first 
sample in your `SummarizedExperiment` object. 
4. Fill in the filtering and normalization parameters and upload a list of 
control genes to use in `RUVg` normalization (not yet implemented) . I'll add 
the control gene functionality later. Right now the control genes used in `RUVg` 
normalization are hard-coded so take the `RUVg` results very lightly. 
5. Click the "Run filtering and normalization" button and *wait*. The app will now
attempt to filter and normalize your data with each method. Upon completion you
should see the filtering tab update with two histograms. Now that the data has 
been processed you can click through the other tabs and compare how the 
different normalization methods affect your data. 

**Please note** that if you change any of the filtering or normalization parameters
on the left-hand panel of the app then you will have re-run the filtering and 
normalization steps again by clicking the "Run filtering and normalization" button. 

If you want to upload a different dataset then it is probably best to refresh the 
page and run through the steps again. 

## How do I create a `SummarizedExperiment` object?

The `SummarizedExperiment` object is just a container for your raw count matrix 
and the metadata attached to your samples. Creating a `SummarizedExperiment` object
is simple if we have these two parts. 

Below, I will create a fake matrix of count data and fake data.frame of metadata 
associated with each sample and then use these two objects to create my 
`SummarizedExperiment`. You will typically already have a count 
matrix and metadata in which case you would simply read these objects into R 
and and create the `SummarizedExperiment` as shown below.

Copy and paste the code below into a fresh R session if 
you are unsure about any of the code.

```r
library(SummarizedExperiment)


# creating fake data -----------------------------------------------------------
# create a fake count matrix
count_matrix <- round(matrix(runif(1e3, 0, 5e3), ncol = 10, nrow = 100))

# name the fake rows gene.1, gene.2, gene.3...etc and
# name the fake columns sample.1, sample,2, sample.3... etc.
dimnames(count_matrix) <- list(paste("gene", 1:100, sep = "."),
                               paste("sample", 1:10, sep = "."))

# create fake sample metadata
# sample.1-5 are 'controls' and sample.6-10 are 'treatment'
metadata <- data.frame(group = rep(c("control", "treatment"), each = 5))
metadata$group <- factor(metadata$group, levels = c("control", "treatment"))
rownames(metadata) <- colnames(count_matrix)

# creating the SummarizedExperiment object -------------------------------------
se <- SummarizedExperiment(assays = list(counts = as.matrix(count_matrix)), 
                           colData = metadata)

# save the se object (this .rds file what gets uploaded into the app)
saveRDS(se, file = "my-experiment.rds")
```

## What does each graph show?

### Filtering

This tab shows a histogram of the count-per-million (CPM) of the count data 
before and after filtering for the specified sample. 

### RLE

This tab shows Relative-log-expression (RLE) boxplots for each sample in the 
dataset. The RLE is calculated by taking the expression value for a single gene
in a sample over the median expression of that gene across samples. 

The RLE boxplots are useful for examining whether or not your data has 
any unwanted variation. Boxplots centered around zero with minimal variation 
indicate appropriate normalization.

### MA (Smear plots)

These plots can be used to compare sample vs sample variation. On the x-axis is
the mean expression and on the y-axis is the log expression of sample2 / sample1. 
The yellow points represented by the smear are representative of genes with very 
low counts, especially those that are zero for one of the samples. Good normalization
is indicated when samples are relatively zero-centered with the mass of expression
being greater than zero.

A smooth scatter can be created to show the density of expression values and a 
lowess line can be plotted to show the trend in expression. For samples in the 
same treatment group, the lowess line should have a slope roughly equal to zero.

### Scatter

These plots show a sample-vs-sample scatter plot. Samples in the same treatment
group should align with a slope roughly equal to 1. Larger deviations may be seen
when comparing samples from different groups. With good normalization there should
not be a global shift away from the 45 degree line.

### PCA

PCA plot on the logCPM normalized data. With adequate normalization your samples
should cluster by treatment group (or whatever factor you're interested in). 

### Corr

`corrplot` of the correlation in expression between the selected samples. Samples
in the same treatment group should correlate highly whereas sames in different
treatment groups should be less correlated.

By default, only the first 10 samples are selected. Use the selectInput widget to
add/remove samples. 

### Heatmap

A heatmap of the row-scaled logCPM expression data. The heatmap function uses the
unsupervised variance to get the top N most variable features and then creates a
heatmap from those features. 

By default, only the first 10 samples are selected. Use the selectInput widget to
add/remove samples. 

**The heatmap may display the 'NA/Inf values in matrix' ERROR** If this happens,
play around with adding more samples or increasing the number of features in the
heatmap. 

## How do I reproduce what I see?

Check out the `helpers.R` file in this repository. It contains all of the functions
used for filtering and normalization as well as plotting. 


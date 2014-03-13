# Plan for Schooner

## Inputs

### Desirata

We want this package to be THE de facto standard for analyzing single-cell
data. Thus, we need to be compatible with a variety of expression data
inputs, such as:

* Cufflinks
* DESeq

For now, we will only be compatible with MISO for splicing scores. I'm not
sure what other programs out there give a 0-1 score for splicing.

### Concerns and Considerations

One sample info file, either a google spreadsheet or a tab-delimited file.
With at least the following columns:

* Sample ID
* Celltype/Group
    * How to deal with within-groups like Motor Neurons and Stressed Motor
    Neurons?
* Gene expression file path
    * however, need to also be compatible with a matrix of gene expression
    values
* Gene expression filetype
    * valid values: our RPKM files, Cufflinks, DESeq files
* Miso summary file path
* Miso splice type

Then the following are optional, but will be filled in with running of the
program.

* Color
    * If they specify "red" or "blue" or something,
    map this to the ColorBrewer Set1 red/blue instead of using the
    `matplotlib` red/blue.
* Plotting symbol

#### Remaining questions

* What about mapping statistics like input reads and mapped reads? Won't
people want to know how many splicing events they got with their highly
sequenced samples versus their lowly sequenced samples?
* We need to consider the tradeoffs between having a really flexible API and
very beginner-user friendly.
    * For example, we would like most of the calls to be just `object.pca()`
    for the most part, but we also want to be flexible to power users who
    want to change something.
* What about downsampling?
* What about the minimum number of cells sharing a particular gene/event?
    * Should specify this in the initialization of the data object.

#### Examples for use

Provide both a fully filled-out spreadsheet with color and plotting symbol,
and a more typical input without those values.

* Need to download sequencing data for existing single cell datasets for
examples.
* Also need to provide data somehow, but don't want to have a gigantic repo.
    * Maybe a supplementary `.zip` file

## Outputs

Overall, want to separate the calculation from the visualization. This way,
can keep the values from the calculation around if they're needed for
something else.

### Data model

Both RPKM and splicing data will be transformed to be in the "tidy" format,
such as:

event_name  sample_id   miso_posterior_mean celltype    gene_name   ensembl

These can be transformed to `event_name` x `sample_id` matrices easily using
`pandas`'s `pivot()` function. But this format makes

It is possible that storing data in memory in this format becomes infeasible
for older computers and we may need to move to a database instead.

### Calculations

User-facing

* Outlier detection
* PCA
* JSD
* Bimodal gene expression
* Splicing modality
* Correlations? Could take a long time..
* Canonical correlation analysis?
* Clustering
* ANOVA?
* RPKM barplots saturation analysis (ask Patrick)

Internal

* Switchy score (needs a new name..)

#### Plotting objects

* PCA
    * Scatterplot with vectors as usual
    * Percent explained by PCs (up to 10-20?)
    * Loadings of each gene in the top 2 PCs
    * Violin plots of the top explaining genes
* JSD
    * Heatmap
    * Violin plots of high JSD events
    * Loadings of high JSD events to compare JSD values
* Lavalamp (you're welcome, Mike)
    * General use for stuff like looking at all events in a celltype,
    or events of different modalities
* Bar plots for binning RPKM events to show high vs low expression
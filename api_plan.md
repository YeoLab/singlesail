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

We must have at a minimum, all the functionality of Fluidigm's "SingulaR" and
 much more.

### Input files

One sample info file, either a google spreadsheet or a tab-delimited file.
With at least the following columns:

* Sample ID
* Celltype/Group
    * How to deal with within-groups like Motor Neurons and Stressed Motor
    Neurons?
* Gene expression file path
    * however, need to also be compatible with a matrix of gene expression
    values
* Splicing filetype
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

#### Tunable parameters

* Number of cells sharing an event/gene for it to be valid
* number of vectors plotted for PCA
* Number of modalities detected?

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
* What about mapping gene IDs to MISO and whatever cufflinks/DESeq uses?
* What about the minimum number of cells sharing a particular gene/event?
    * Should specify this in the initialization of the data object.
* Filtering miso events? Using the ci_halves_max or just ci_diff? How to
expose this to the user?
* how to integrate other data such as conservation of exons or genes? Exon
info files?
* What about the Biomark data? Are we just going to not support it?
* How to output the data of each plot such that users can change the
visualization themselves?
* Statistical analysis to be confident that certain events truly ARE high or
low variance.

#### Examples for use

Provide both a fully filled-out spreadsheet with color and plotting symbol,
and a more typical input without those values.

* Need to download sequencing data for existing single cell datasets for
examples.
* Also need to provide example data somehow,
but don't want to have a gigantic repo.
    * Maybe a supplementary `.zip` file
    * Host the files on Sauron? Maybe the file paths in the google
    spreadsheet could be a URL...


## Outputs

Overall, want to separate the calculation from the visualization. This way,
can keep the values from the calculation around if they're needed for
something else.

### Data model

Both RPKM and splicing data will be transformed to be in the "tidy" format,
such as:

    event_name  sample_id   miso_posterior_mean celltype    gene_name   ensembl

These can be transformed to `event_name` x `sample_id` matrices easily using
`pandas`'s `pivot()` function. But this format makes it very easy for
creating violin plots and counting numbers of events via `groupby`s.

It is possible that storing data in memory in this format becomes infeasible
for older computers and we may need to move to a database instead.

### Calculations

User-facing

* PCA
* JSD
* Bimodal gene expression
* Splicing modality 
* *de novo* modality clustering
* Correlations? Could take a long time..
* Canonical correlation analysis?
* Clustering
* ANOVA?
* RPKM barplots saturation analysis (ask Patrick)
* histograms of psi scores?
* histograms of expression?
* Scatterplots of rpkm vs psi?
    * High variance events vs RPKMs
* Percent mapping

Internal

* Switchy score (needs a new name..)
* Outlier detection (visualized on PCA)

#### Plotting objects

Except when specified, these will be used for both splicing and gene
expression.

* PCA (This should be in seaborn)
    * Scatterplot with vectors as usual
    * Percent explained by PCs (up to 10-20?)
    * Loadings of each gene in the top 2 PCs
    * Violin plots of the top explaining genes
    * Show outliers
* JSD (seaborn?)
    * Heatmap
    * Violin plots of high JSD events
    * Loadings of high JSD events to compare JSD values
* Lavalamp (you're welcome, Mike)
    * General use for stuff like looking at all events in a celltype,
    or events of different modalities
* Bar plots for binning RPKM events to show high vs low expression
* Boxplots of expression for each sample (like the bar plots above)
* Clusterplot
* Splicing modality detection
    * Lavalamp of the different modalities
    * PCA

## Example/Tutorial

Need to have a really complete example ipython notebook. Should have these
actions in the following order:

1. Install library
2. Import library
3. Instantiate data object
    1. Filter gene expression and splicing events based on number of cells
    sharing this event in that celltype?
4. Gene Expression analysis
    1. Histograms of expression for each cell (sharex=True, sharey=True)
    2. Breakdown of RPKMs of different events.
    3. Detect outliers
    4. Show outliers on gene expression PCA
    5. Show clusterplot of genes?
    6. Detect bimodal gene expression, show violin plots
    7. Detect high-JSD events, show violin plots
    8. ANOVA?
5. Splicing analysis
    1. Histogram of splicing scores for each cell
    2. Subset on variance, high medium low?
    1. Lava lamp plots of all events.
    2. Show PCA
    3. Detect modalities
        1. Histograms of splicing scores within modalities
        2. Lava lamps
        3. PCA?
        4. Barplots of modalities
    4. High JSD events, show violin plots
6. Splicing + Gene expression analysis
    1. Scatterplot of RPKMs vs expression (need miso event to ensembl ID
    matching)
    2. Mike's RBPs and splicing analysis?

## Python dependencies

* pandas (numpy, scipy, matplotlib)
* matplotlib
* gspread
* brewer2mpl
* seaborn (statsmodels)
* IPython?

## Naming scheme

* CamelCase for classes only
* separated_by_underscores for functions
* Full names as much as possible. E.g. `gene_symbol` instead of `gsymbol`

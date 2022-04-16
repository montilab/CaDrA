
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CaDrA

Candidate Drivers Analysis: Multi-Omic Search for Candidate Drivers of
Functional Signatures

CaDrA is an R package that supports a heuristic search framework aimed
towards the identification of candidate drivers of oncogenic activity.
Given a binary genomic dataset (where the rows are 1/0 vectors
indicating the presence/absence of genomic features such as somatic
mutations or copy number alteration events), together with an associated
sample ranking (where the samples are ranked by a certain phenotypic
readout of interest such as protein expression, pathway activity etc.),
CaDrA implements a step-wise search algorithm to determine a set of
features that, together (based on their occurence union or ‘logical
OR’), is most-associated with the observed ranking, making it useful for
finding mutually exclusive or largely non-overlapping anomalies that can
lead to the same pathway phenotype.

For more information, please see the associated manuscript [Kartha et
al. (2019)](https://www.frontiersin.org/articles/10.3389/fgene.2019.00121/full)

## (1) Installation

    devtools::install_github("montilab/CaDrA")
    library(CaDrA)

## (2) Quickstart

### Test run code on simulated data

    data(sim.ES)
    data(topn.list)

    # Plot the results from a top-N evaluation by passing the resulting ESet from a specific run
    # To find the combination of features that had the best score
    best.meta <- topn.best(topn.list)

    # Now we can plot this set of features
    meta.plot(best.meta$ESet)

### Running on your own (actual) data

First, you need to have your data (we’ve used copy number variation,
mutation calls or methylation data - anything that can basically be
thresholded / binarized as 0/1) in an `ExpressionSet` object. In the
example below, the expression set object is named `ES.GISTIC.Mut`. The
only other thing you need is a variable used to rank the matrix by
(below, we have a variable called `score`). The ranking of this variable
is what dictates how CaDrA will search for grouped meta-features. By
default, it will look for things that are left-skewed, so you can rank
the variable in either ascending order or descending order to look for
features that are enriched at either end (in the example below, we are
interested in looking at features enriched for higher scores in samples,
hence ranking them in descending order of score first)

The function calls and pre-processing example steps are as follows:

    # Binarize ES to only have 0's and 1's
    exprs(ES.GISTIC.Mut)[exprs(ES.GISTIC.Mut)>1] <- 1

    # Pre-filter ESet based on occurrence frequency
    ES.GISTIC.Mut.filt <- prefilter_data(ES = ES.GISTIC.Mut,max.cutoff = 0.7,min.cutoff = 0.03) # no more than 70%, no less than 3%

    # Order of samples in decreasing order of supplied continuous variable
    sample.order <- order(score,decreasing=TRUE)

    # Number of top starting seed features to test and evaluate over  
    top_n <- 7

    # Metric used for stepwise greedy search
    # Either ks or wilcox is supported
    method <- "ks"

    topn.l <- topn.eval(ESet = ES.GISTIC.Mut.filt, 
                        method=method,
                        N = top_n,
                        do.plot = FALSE, #We will plot it AFTER finding the best hits
                        best.score.only = FALSE,
                        ranking=sample.order,
                        verb=FALSE)

    # Now we can fetch the ESet and feature that corresponded to the best score over the top N search
    topn.best.meta <- topn.best(topn.l)

    # Visualize best result
    meta.plot(ESet = topn.best.meta$ESet,
              var.score = score[sample.order],
              var.name = "Activity score") #Y-axis label for plot


    # You can also evaluate how robust the results are depending on which seed feature you started with

    topn.plot(topn.l) 

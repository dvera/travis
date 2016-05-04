# travis

**R** **U**tilities for **B**ed and **BE**dg**R**aphs

---

## Overview

travis is a suite of tools to process, analyze, and visualize bed and bedGraph data in R. Most of travis's functions can be categorized as either operating on beds or on bedGraphs.

Functions that operate on bed files:
- `bedBreadth` calculates the total number of bases covered in a bed file
- `bedCat` concatenates multiple bed files into one
- `bedCenters` takes a bed file and outputs a bed file of the midpoints of the intervals
- `bedDivideMeta` divides bed intervals (and optionally surrounding region) into an equal number of windows per interval, regardless of size
- `bedDivide` divides equally-sized bed invervals into an equal number of windows per interval 
- `bedEnds` takes a bed file and outputs a bed file of the ends of the intervals
- `bedHist` plots a histogram of interval sizes
- `bedOverhangs` calculates frequencies of pairwise distances between the ends of intervals in a single bed file 
- `bedOverlap` calculates overlap among two sets of bed files
- `bedParseLengths` parses bed intervals into separate files based on interval sizes
- `bedParseStrands` prases bed intervals into separate files based on strand
- `bedPoissonTest` perform a poisson test of interval counts between two bed files with a given set of windows
- `bedRecenter` alter the coordinates of a bed file
- `bedRemoveChrom` remove intervals on a specific chromosome
- `bedSample` randomly sample intevals of a bed file
- `bedSizes` load into R sizes of bed intervals
- `bedSort` sort a bed file by chromosome/coordinate
- `bedStructures` break up a bed12 file into separate components
- `bedWords` calculate nucleotide word frequencies for intervals in a bed file

Functions that operate on bedGraph files:
- `bgAcf` calculates autocorrelations of bedGraph files
- `bgCorMatrix` generates a matrix of pairwise correlations among a set of bedGraph files
- `bgHist` plots a histogram of bedGraph scores
- `bgIqrNorm` normalizes bedGraphs to have a specified interquartile range
- `bgIqr` calculate interquartile ranges of bedGraph files
- `bgKsTest` performs a KS test between two bedGraph files
- `bgLoess` loess smooths a bedGraph files
- `bgOps` perform various operations on single bedGraph files (e.g. transformations) or multiple bedGraph files (averages, differences, variance, ratios)
- `bgParseScores` parses bedGraph intervals into multiple bedGraphs based on score
- `bgPlot` plot bedGraph scores of specified intervals
- `bgPoissonTest` perform a poisson test between two bedGraph files
- `bgQuantileNorm` normalizes bedGraphs with quantile normalization
- `bgQuantiles` calculates specified quantiles of bedGraph files
- `bgRead` reads bedGraph files into R
- `bgScatterGrid` calculate a grid of pairwise density scatter plots using a set of bedGraph files
- `bgScores` read bedGraph scores into R
- `bgSort` sort bedGraph files by score
- `bgThreshold` outputs bedGraphs containing intervals exceeding a specified score
- `bgUnify` make multiple bedGraph files share an identical set of intervals 
- `bgZtest` perform a Z-test between two bedGraph files

There are several functions that use bam or sam files:
- `bamCount` counts the number of reads in bam files
- `bamParseLengths` parses paired-end reads into multiple files based on insert size
- `samStats` calcaultes and plots summary statistics for alignment data in sam files

Finally, there are three plotting functions that operate directly on R objects:
- `intervalBox` draws a boxplot of data that is binned based on its score in another dimension
- `rageHist` plots overlayed histograms of data organized in lists or data frames
- `scatterdens` generates a density plot of points in a scatter plot, useful for generating a scatter plot of many points

## Installation

travis heavily depends on [bedtools](https://github.com/arq5x/bedtools2). Other programs travis uses include several GNU programs that are typically preinstalled in most Linux distributions, including `awk`, `sed`, `sort`, `shuf`, and `cut`. For older linux distributions, you may need to install a more recent version of [GNU core utilities](http://www.gnu.org/software/coreutils/coreutils.html) that supports parallelization in `sort`.

In addition, there are R dependencies. Install devtools if not installed already:
```R
install.packages("devtools")
```

Then install conifur (convenience functions for R), converge (conversion tools for genomic data), gyro (wrappers for genomic tools), and travis:
```R
devtools::install_github("dvera/conifur")
devtools::install_github("dvera/converge")
devtools::install_github("dvera/gyro")
devtools::install_github("dvera/travis")
```

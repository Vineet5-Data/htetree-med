# htetree: Causal Inference with Tree-Based Machine Learning Algorithms

[![CRAN status](https://www.r-pkg.org/badges/version/htetree)](https://CRAN.R-project.org/package=htetree)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/htetree)](https://CRAN.R-project.org/package=htetree)

## Overview

`htetree` is an R package for estimating heterogeneous treatment effects (HTE) with tree-based machine learning algorithms and visualizing estimated results in flexible and presentation-ready ways. The package provides powerful tools for causal inference analysis using decision trees and forest-based methods.

## Key Features

- **Causal Trees**: Estimate heterogeneous treatment effects using tree-based methods
- **Causal Forests**: Forest-based approaches for more robust HTE estimation
- **Matching Methods**: Propensity score matching within leaves
- **IPW Methods**: Inverse probability weighting for treatment effect estimation
- **Visualization Tools**: Flexible plotting functions for presenting results
- **Interactive Interface**: Shiny app for dynamic exploration of results

## Installation

### From CRAN

You can install the released version of `htetree` from CRAN:

```r
install.packages("htetree")
```

### Development Version

To install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("Vineet5-Data/htetree-med")
```

## Quick Start

```r
library(htetree)

# Load example data
data(simulation.1)

# Fit a causal tree
tree <- causalTree(y ~ x1 + x2 + x3 + x4 + x5,
                   data = simulation.1,
                   treatment = simulation.1$treatment,
                   split.Rule = "CT",
                   cv.option = "CT",
                   split.Honest = TRUE,
                   cv.Honest = TRUE)

# Plot the tree
plot(tree)
text(tree)

# Estimate heterogeneous treatment effects
hte_results <- hte_causalTree(tree, data = simulation.1)
```

## Main Functions

- `causalTree()`: Fit causal trees for HTE estimation
- `causalForest()`: Fit causal forests for robust HTE estimation
- `hte_causalTree()`: Estimate HTE using causal trees
- `hte_forest()`: Estimate HTE using forest methods
- `hte_match()`: Estimate HTE using matching methods
- `hte_ipw()`: Estimate HTE using inverse probability weighting
- `hte_plot()`: Visualize HTE estimates
- `runDynamic()`: Launch interactive Shiny app

## Documentation

For detailed documentation and examples, see:
- Package documentation: `help(package = "htetree")`
- Function help: `?causalTree`, `?hte_causalTree`, etc.

## Citation

If you use this package in your research, please cite:

Brand, J.E., Xu, J., Koch, B., and Geraldo, P. (2021). "Uncovering Sociological Effect Heterogeneity using Tree-Based Machine Learning." *Sociological Methodology*, 51(2), 189-223. doi:10.1177/0081175021993503

```r
citation("htetree")
```

## Authors

- **Jiahui Xu** (Maintainer) - jiahuixu@ucla.edu
- **Tanvi Shinkre** - tanvishinkre@ucla.edu
- **Jennie Brand** - brand@soc.ucla.edu

## License

This package is licensed under GPL-2 | GPL-3.

## Acknowledgments

This package started as a fork of the `causalTree` package on GitHub. We greatly appreciate the original authors for their extremely useful and free package.

## Support

For bug reports and feature requests, please open an issue on the [GitHub repository](https://github.com/Vineet5-Data/htetree-med/issues).

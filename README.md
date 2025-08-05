# OmixM6A

**OmixM6A** offers a suite of computational tools for the analysis and visualization of m6A omics data in R. The package specializes in the site calling, classification and normalization of m6A methylation states within quantitative m6A samples, providing researchers with the tools needed to conduct high-level analysis in the field of epitranscriptomics.

## Installation

Install OmixM6A directly from GitHub:

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install OmixM6A from GitHub
devtools::install_github("ZW-xjtlu/OmixM6A")
```

## Overview

OmixM6A provides:

- **Rich set of statistical models**: Fit a range of models, including (zero & one inflated) beta-binomial mixtures, binomial mixtures, binomial-uniform mixtures, and simple binomial, to m6A count data.
- **Flexible data handling**: Efficiently analyze both individual count vectors and `SummarizedExperiment` objects from m6AConquer database, allowing for versatile data analysis workflows.
- **Robust classification and normalization**: Accurately classify and normalize m6A methylation states, calculating posterior probabilities and p-values.
- **High-level visualization**: Generate visualizations to compare model fits, assess goodness of fit, and interpret the results of different statistical models.
- **Scalability**: Parameter initiation of models are optimized for large-scale m6A datasets, ensuring efficient processing and analysis for epitranscriptome research with million number of sites / peaks.

## Key Features

- **P-Value calculation**: Generate p-values to assess the statistical significance of m6A sites in site calling.
- **Posterior probability calculation**: Estimate the posterior probabilities of m6A modification states across biological samples using selected statistical models, which are useful as bayes classifier for methylation states or normalized methylation levels across platforms.
- **Model comparison and visualization**: Compare and evaluate the fit of various statistical models through Q-Q plots, BIC scores, and marginal alignments, enhancing transparency in model selection.
- **Versatile data processing**: Handle both simple data inputs and m6AConquer data-sharing framework within `SummarizedExperiment` objects, making the package adaptable to different research needs.

## Usage

```r
# Load the OmixM6A package
library(OmixM6A)
library(SummarizedExperiment)

# Example usage with count vectors
# Load an example SummarizedExperiment object
m6A_se <- readRDS(system.file("extdata", "example_se.rds", package="OmixM6A"))

# Extract m6A and total counts from the SummarizedExperiment object
m6A_counts <- assays(m6A_se)$m6A[,1]
total_counts <- assays(m6A_se)$Total[,1]

# Apply OmixM6A to count vectors (fitting BBmix/beta-binomial mixture)
result_df <- OmixM6A(m6A_counts, total_counts, method = "bbmix") 

# You can set method = "binomial" to reproduce the binomial test used by default in m6AConquer
# result_df <- OmixM6A(m6A_counts, total_counts, method = "binomial") 

# Display the results
print(result_df)

# Example usage with a SummarizedExperiment object
# Apply OmixM6A directly to the SummarizedExperiment object
result_se <- OmixM6A(se = m6A_se, method = "bbmix")

# Display the results
print(result_se)

# Check fitted model parameters stored in metadata
metadata(result_se)
```
## Documentation

Comprehensive documentation for each function is available within the package. Access it using:

```r
help(package = "OmixM6A")
```

## Contributing

Contributions to OmixM6A are welcome! Report bugs, suggest features, or contribute code by [creating an issue](https://github.com/ZW-xjtlu/OmixM6A/issues) or submitting a pull request.

## License

This package is licensed under the MIT License. See the [LICENSE](https://github.com/ZW-xjtlu/OmixM6A/blob/main/LICENSE) file for details.

## Acknowledgments

OmixM6A is developed to be compatible with the m6AConquer database project. 

We appreciate the contributions and feedback from the epitranscriptomics community.

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

- **Advanced Statistical Models**: Fit a range of models, including zero-inflated beta-binomial mixtures, beta-binomial mixtures, binomial mixtures, and binomial-uniform mixtures, to m6A count data.
- **Flexible Data Handling**: Efficiently analyze both individual count vectors and `SummarizedExperiment` objects from m6AConquer database, allowing for versatile data analysis workflows.
- **Robust Classification and Normalization**: Accurately classify and normalize m6A methylation states, calculating posterior probabilities, p-values, and false discovery rates (FDR).
- **High-Level Visualization**: Generate visualizations to compare model fits, assess goodness of fit, and interpret the results of different statistical models.
- **Scalability**: Parameter initiation of models are optimized for large-scale m6A datasets, ensuring efficient processing and analysis for epitranscriptome research with million number of sites / peaks.

## Key Features

- **P-Value and FDR Calculation**: Generate p-values and false discovery rates to assess the statistical significance of m6A sites in site calling.
- **Posterior Probability Calculation**: Estimate the posterior probabilities of m6A modification states across biological samples using selected statistical models, which are useful as bayes classifier for methylation states or normalized methylation levels across platforms.
- **Model Comparison and Visualization**: Visual tools to compare and evaluate the fit of various statistical models through multiple metrics, enhancing transparancy in model selection.
- **Versatile Data Processing**: Handle both simple data structures and multi-sample data framework within `SummarizedExperiment` objects, making the package adaptable to different research needs.

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

# Apply OmixM6A to count vectors (fitting BBmix/beta-binomial mixture by default)
result_df <- OmixM6A(m6A_counts, total_counts)

# Display the results
print(result_df)

# Example usage with a SummarizedExperiment object
# Apply OmixM6A directly to the SummarizedExperiment object
result_se <- OmixM6A(se = m6A_se)

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

OmixM6A is developed as part of the m6AConquer database project, with the goal of advancing the analysis and understanding of m6A RNA modifications through systematic integration. We appreciate the contributions and feedback from the epitranscriptomics community.

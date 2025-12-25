# Random-Effects Approach to Generalized Linear Mixed Model Analysis of Incomplete Longitudinal Data Package

In this document, we illustrate the main features of the `cdeMissing` R package through examples. Additional information on the statistical methodology and computational details are provided in the accompanying documentation and research articles.

## Cite the package

The package applies methods introduced in the [paper](https://doi.org/10.1002/sim.70343):

T. Nguyen, J. Zhang, and J. Jiang, “ A Random-Effects Approach to Generalized Linear Mixed Model Analysis of Incomplete Longitudinal Data,” Statistics in Medicine 44, no. 28-30 (2025): e70343, https://doi.org/10.1002/sim.70343.

## Install

Open the R console and run the following command to install the package from source:

```r
install.packages("devtools") # When you have not installed devtools package
devtools::install_github("Celaeno1017/cdeMissing")
```

## Tutorial

First, load the R package.

```r
library(cdeMissing)
```

To illustrate the main features of the R package `cdeMissing`, let's first generate some data. We have built in a few functions directly into the R package for this purpose.

```r
######lmm#########
dat<-sim_cde_data(model = 'lmm',beta = c(1,0.5,0.2,0.2,0.2))
#####glmm with binomial response##########
dat_glmm<-sim_cde_data(model = 'glmm',beta = c(-1,0.5,0.2,0.2,0.2))
```

To fit the lmm/glmm model with missingness treated as random effects, we use the cde_fit function which implements the proposed methodology.

```r
####lmm###########
fit <- cde_fit(
    formula=y ~ x1 + Treatment + Time + x2,
    data = dat$data_mis,
    id = "id",
    time = "Time",
    random = ~1|id,
    correlation=NULL,
    missing = list(baseline = c("x1"), time_dependent = c("x2")),
    family = gaussian(),
    engine = "nlme")

####glmm with binomial response###########
fit <- cde_fit(
    formula=y ~ x1 + Treatment + Time + x2,
    data = dat$data_mis,
    id = "id",
    time = "Time",
    random = ~1|id,
    correlation=NULL,
    missing = list(baseline = list("x1"), time_dependent = list("x2")),
    family = binomial(),
    engine = "lme4")
```
Key parameters include:
- `formula`: The fixed effect model formula for implementation.
- 'random`: The random effect for implementation. Should be written in `nlme` package style.
- `data`: The dataframe which is in long format contains id, response, and covariates.
- `id`: The variable in the data indicates the id of subjects. Cannot be NULL.
- `time`: The variable in the data indicates the time for longitudinal data. If is Null, make sure there is no time_dependent variable input in `missing`.
- `missing`: Which variables contain missingness. Here, you need to indicate whether the variable is baseline or time_dependent.
- `correlation`: additional correlation structure passed to `nlme` package for lmm. Not available for glmm.
- `family`: The error distribution and link function to be used in the model. For lmm, it should be gaussian.
- `engine`: The package performs model fitting. for lmm, use `nlme`; for glmm, use `lme4`.

The result will contain the following values:
- `fit`:	the final fitted parameters of lmm/glmm model.


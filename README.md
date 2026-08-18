# README


Code and data for reproducing the analyses presented in:

> Author et al. (YEAR). Bayesian model-stacking improves somatic growth 
estimates in Everglades Cyprinodontid fishes.[Journal, DOI]

This repository contains the R code and data products required to
reproduce the analyses, figures, and tables presented in the manuscript.


For questions about this data set or analysis, please contact:

```bash
Name: William K. Annis

Email: wannis@fsu.edu, williamkannis@gmail.com

OrcID: 0009-0003-3541-8503
```

If you use this code or data, please cite:

>BLANK. Bayesian model-stacking improves somatic growth estimates in 
Everglades Cyprinodontid fishes. in review

Additionally, if you use associated R functions, also cite:

>BLANK


## Data
**Directory:** ```data/```

## Model outputs

## Figures
**Directory:** ```figures/```

Contains raw R plots and .csv tables used to create the figures and tables in 
the main manuscript and appendices. Most figures were edited in Adobe Illustrator
for aesthetic purposes, and we included Illustrator files as well.


## Workflow
**Directory:** ```scripts/```

The following scripts provide all R code necessary to replicate the 
manuscript's analyses, and are named sequentially in workflow order. Scripts
call in custom functions designed to fit growth models, make predictions, and
summarize results. These functions and their documentation can be found 
[here](BLANK).

### Install custom functions
<br>**Script**: ```00_install_growthstack_pkg.R```

<ins>Purpose:</ins> Installs package contain custom functions used for
the manuscript's analyses.

### Data preparation 

<br>**Script**: ```01_age_length_data_cleaning.R```

<ins>Purpose:</ins> Prepares otolith-derived age-at-length data for growth 
modelling. Retains only female specimens for consistency among species, and
removes missing data.

<ins>Output:</ins> X

<br>**Script**: ```02_growth_predictor_data_prep.R```

<ins>Purpose:</ins> Prepares growth parameter predictor variables, and conducts
principal component analysis (PCA) to create composite hydrology variables.

<ins>Output:</ins> X


### Stan model fitting

<br>**Script**: ```03_stan_batch_run_linear.R```

<ins>Purpose:</ins> Fits continuous second-level predictor models of all
three growth forms to each species. Random effect only model is fit to Flagfish
(JORFLO) due to smaller sample size.

<ins>Output:</ins> X

<br>**Script**: ```03_stan_batch_run_cat.R```

<ins>Purpose:</ins> Identical as ```03_stan_batch_run_linear.R```,
but for categorical predictor model.

<ins>Output:</ins> X

### Model stacking and predictions

<br>**Script**: ```04_batch_loo_stacking_linear.R```

<ins>Purpose:</ins> X

<ins>Output:</ins> X

<br>**Script**: ```04_batch_loo_stacking_cat.R```

<ins>Purpose:</ins> Identical as ```04_batch_loo_stacking_linear.R```,
but for categorical predictor model outputs.

<ins>Output:</ins> X

### Model summary statistics

<br>**Script**: ```05_model_stacking_summary_stats.R```

<ins>Purpose:</ins> X

<ins>Output:</ins> X

<br>**Script**: ```05_model_stacking_summary_stats_cat.R```

<ins>Purpose:</ins> Identical as ```05_model_stacking_summary_stats.R```,
but for categorical predictor model outputs.

<ins>Output:</ins> X

### Plotting

<br>**Script**: ```06_growth_curve_plots.R```

<ins>Purpose:</ins> X

<ins>Output:</ins> X

<br>**Script**: ```06_growth_curve_plots_cat.R```

<ins>Purpose:</ins> Identical as ```06_growth_curve_plots.R```,
but for categorical predictor model outputs.

<ins>Output:</ins> X

### Model fit
<br>**Script**: ```07_model_fit_appendix.R```

<ins>Purpose:</ins> X

<ins>Output:</ins> X

<br>**Script**: ```07_model_fit_appendix_cat.R```

<ins>Purpose:</ins> Identical as ```07_model_fit_appendix.R```,
but for categorical predictor model outputs.

<ins>Output:</ins> X

## Growth models
**Directory:** ```stan/```

We estimated growth rates using three model forms: von Bertalanffy, Gompertz, 
and Logistic; and three effect structures: Random effect only, Categorical 
second-level predictors, and continuous second-level predictors. While model
stan files were called in using the [associated R package](BLANK),
we also provide the stan scripts for each model code in this repository.
Additionally, we provided R script used to validate each model with simulated
data.See below for more information of growth forms and model type:

### Growth forms

**von Bertalanffy (vb)**

**Logistic (lg)**

**Gompertz (gz)**

### Random effect only model
**Files:** ```vb_BLANK.stan```, ```gz_BLANK.stan```, ```lg_BLANK.stan```

### Categorical predictors
**Files:** ```vb_BLANK.stan```, ```gz_BLANK.stan```, ```lg_BLANK.stan```

### continuous predictors
**Files:** ```vb_BLANK.stan```, ```gz_BLANK.stan```, ```lg_BLANK.stan```

### Model dvelopment
**Script:** ```BLANK.R```

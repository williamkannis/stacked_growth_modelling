# Source code for *Bayesian model-stacking for improved somatic growth modelling: A case study examining hydrology-dependent growth of six cyprinodontid fish species in the Florida Everglades*

Code and data for reproducing the analyses presented in:

> Author et al. (YEAR). Syndromes of multidimensional beta diversity change in
invaded metacommunities.[Journal, DOI]

This repository contains the R code and data products required to
reproduce the analyses, figures, and tables presented in the manuscript.


For questions about this data set or analysis, please contact:

```bash
Name: William K. Annis

Email: wannis@fsu.edu, williamkannis@gmail.com

OrcID: 0009-0003-3541-8503
```

If you use this code or data, please cite:

>BLANK. Syndromes of multidimensional beta diversity change in 
invaded metacommunities. in review


## Data
**Directory:** ```data/```

## Model outputs

## Figures
**Directory:** ```figures/```

## Growth models
**Directory:** ```stan/```

We estimated growth rates using 

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

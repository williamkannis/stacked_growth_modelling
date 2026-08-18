# Source code for *Bayesian model-stacking for improved somatic growth modelling: A case study examining hydrology-dependent growth of six cyprinodontid fish species in the Florida Everglades*

## Contact information and citation

``` bash
Name: William K. Annis

Email: wannis@fsu.edu, williamkannis@gmail.com

OrcID: 0009-0003-3541-8503
```

Cite as: 
> CITE

## Data

## Figures

## Growth models
We estimated growth rates using 

## Workflow

### Data preparation 

**Script**: ```01_age_length_data_cleaning.R```

<ins>Purpose:</ins> X

<ins>Output:</ins> X

**Script**: ```02_growth_predictor_data_prep.R```

<ins>Purpose:</ins> X

<ins>Output:</ins> X


### Stan model fitting

**Script**: ```03_stan_batch_run_linear.R```

<ins>Purpose:</ins> X

<ins>Output:</ins> X

**Script**: ```03_stan_batch_run_cat.R```

<ins>Purpose:</ins> Identical as ```03_stan_batch_run_linear.R```,
but for categorical predictor model.

<ins>Output:</ins> X


### Model stacking and predictions

**Script**: ```04_batch_loo_stacking_linear.R```

<ins>Purpose:</ins> X

<ins>Output:</ins> X

**Script**: ```04_batch_loo_stacking_cat.R```

<ins>Purpose:</ins> Identical as ```04_batch_loo_stacking_linear.R```,
but for categorical predictor model outputs.

<ins>Output:</ins> X

### Model summary statistics

**Script**: ```05_model_stacking_summary_stats.R```

<ins>Purpose:</ins> X

<ins>Output:</ins> X

**Script**: ```05_model_stacking_summary_stats_cat.R```

<ins>Purpose:</ins> Identical as ```05_model_stacking_summary_stats.R```,
but for categorical predictor model outputs.

<ins>Output:</ins> X

### Plotting

**Script**: ```06_growth_curve_plots.R```

<ins>Purpose:</ins> X

<ins>Output:</ins> X

**Script**: ```06_growth_curve_plots_cat.R```

<ins>Purpose:</ins> Identical as ```06_growth_curve_plots.R```,
but for categorical predictor model outputs.

<ins>Output:</ins> X

### Model fit
**Script**: ```07_model_fit_appendix.R```

<ins>Purpose:</ins> X

<ins>Output:</ins> X
<brk>
**Script**: ```07_model_fit_appendix_cat.R```

<ins>Purpose:</ins> Identical as ```07_model_fit_appendix.R```,
but for categorical predictor model outputs.

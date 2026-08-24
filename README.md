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

If you use these, models, code or data, please cite:

>BLANK. Bayesian model-stacking improves somatic growth estimates in 
Everglades Cyprinodontid fishes. in review

Additionally, if you use associated R functions, also cite:

>BLANK


## Data
**Directory:** ```data/```

The data below can be found at the manuscript's [Zenodo repository:](BLANK)
These can be used to replicate the full manuscript workflow

### Age at length data
**Files:**

### Species key
**File:**

### Predictor data
**Files:**

## Outputs
**Directory:** ```outputs/```

File containing the following outputs can be found at the manuscript's 
[Zenodo repository:](BLANK). These can be used to replicate specific
sections of the workflow, as indicated below:

### Growth model outputs
**Files:**

### Loo outputs
**Files:**

### Model stacking output
**Files:**

### Curve predictions
**Files:**

## File directory
Download the entire repository. Then download required data from Zenodo 
Repository, and data sources listed here, and unzip into the following file 
structure:

```bash

│
├── data
│   │── raw_data
│   │   │── fs_age.rds*
│   │   │── fs_predictors.rds*
│   │   └── fs_species_key.csv*
│   │
│   └── analysis_data
│       │── fs_age_final.rds*
│       └── fs_pred_final.rds*
│
├── outputs
│   │── stan_outputs*
│   │── loo_outputs*
│   └── curve_predictions*
│
├── scripts
├── stan_scripts
└── figures


(*) directories or files downloaded from Zenodo
```


## Workflow
**Directory:** ```scripts/```

The following scripts provide all R code necessary to replicate the 
manuscript's analyses, and are named sequentially in workflow order. Scripts
call in custom functions designed to fit growth models, make predictions, and
summarize results. These functions and their documentation can be found 
[here](BLANK). 
Script whose numerical id is followed by "b" refer to scripts used in
relation to categorical (cat) effect version of the models. These only
need to be ran to replicate the results in Appendix 2.

### 1. Install custom functions
<br>**Script**: ```00_install_growthstack_pkg.R```

<ins>Purpose:</ins> Installs package contain custom functions used for
the manuscript's analyses.

### 2. Data preparation 

<br>**Script**: ```01_age_length_data_cleaning.R```

<ins>Purpose:</ins> Prepares otolith-derived age-at-length data for growth 
modelling. Retains only female specimens for consistency among species, and
removes missing data.

<ins>Output:</ins> X

<br>**Script**: ```02_growth_predictor_data_prep.R```

<ins>Purpose:</ins> Prepares growth parameter predictor variables, and conducts
principal component analysis (PCA) to create composite hydrology variables.

<ins>Output:</ins> X


### 3. Stan model fitting

<br>**Script**: ```03_stan_batch_run.R```

<ins>Purpose:</ins> Fits continuous second-level predictor models of all
three growth forms to each species. Random effect only model is fit to Flagfish
(JORFLO) due to smaller sample size. Creates age and length summary tables.

<ins>Output:</ins> X

<br>**Script**: ```03b_stan_batch_run_cat.R```

<ins>Purpose:</ins> Identical as ```03_stan_batch_run.R```,
but for categorical predictor model.

<ins>Output:</ins> X

### 4. Model stacking and predictions

<br>**Script**: ```04_batch_loo_stacking.R```

<ins>Purpose:</ins> X

<ins>Output:</ins> X

<br>**Script**: ```04b_batch_loo_stacking_cat.R```

<ins>Purpose:</ins> Identical as ```04_batch_loo_stacking.R```,
but for categorical predictor model outputs.

<ins>Output:</ins> X

### 5. Model fit
<br>**Script**: ```05_model_fit_appendix.R```

<ins>Purpose:</ins> X

<ins>Output:</ins> X

<br>**Script**: ```05b_model_fit_appendix_cat.R```

<ins>Purpose:</ins> Identical as ```05_model_fit_appendix.R```,
but for categorical predictor model outputs.

<ins>Output:</ins> X

### 6. Model summary statistics

<br>**Script**: ```06_model_stacking_summary_stats.R```

<ins>Purpose:</ins> X

<ins>Output:</ins> X

<br>**Script**: ```06b_model_stacking_summary_stats_cat.R```

<ins>Purpose:</ins> Identical as ```06_model_stacking_summary_stats.R```,
but for categorical predictor model outputs.

<ins>Output:</ins> X

### 7. Plotting

<br>**Script**: ```07_growth_curve_plots.R```

<ins>Purpose:</ins> X

<ins>Output:</ins> X

<br>**Script**: ```07b_growth_curve_plots_cat.R```

<ins>Purpose:</ins> Identical as ```07_growth_curve_plots.R```,
but for categorical predictor model outputs.

<ins>Output:</ins> X


## Growth models
**Directory:** ```stan/```

We estimated growth rates using three model forms: von Bertalanffy, Gompertz, 
and Logistic; and three effect structures: Random effect only, Categorical 
second-level predictors, and continuous second-level predictors. While model
stan files were called in using the [associated R package](BLANK),
we also provide the stan scripts for each model code in this repository. Stan
files are named by growth form and effect structure as defined below.
Additionally, we provided R script used to validate each model with simulated
data. See below for more information of growth forms and model type:

### Growth forms
For each growth form, the equation for length (L) at age (t), and the differential 
equation for instantaneous growth (G) at length are given below. Asymptote terms 
(L<sub>∞</sub>) describe the maximum length for the average fish, the scaling terms 
(g<sub>1-3</sub>)  describe the slope of the growth curve, and the inflection term 
(t<sub>0</sub>, t<sub>inf</sub>) describes the age at which maximum growth rate occurs

**von Bertalanffy (vb)**
$$L(t) = L_{\infty} (1 - e^{-g_1 (t - t_0)})$$
$$G(L) = g_1 (L_{\infty} - L)$$

**Gompertz (gz)**
$$L(t) = L_{\infty} e^{-e^{-g_2 (t - t_{\text{inf}})}}$$
$$G(L) = g_2 \times L \times \ln\left(\frac{L_{\infty}}{L}\right)$$

**Logistic (lg)**
$$L(t) = \frac{L_{\infty}}{1 + e^{-g_3 (t - t_{\text{inf}})}}$$
$$G(L) = g_3 \times L \left(1 - \frac{L}{L_{\infty}}\right)$$

### Effect structure

**Random effect only model (random)**
$$L_{i,j} = \text{FUN}_x (\text{age}_{i,j} \mid L_{\infty j}, g_j, t_j) + \varepsilon_{i,j}$$

$$\varepsilon_{i,j} \sim \text{StudentT}(\nu, 0, \sigma^2)$$

$$\log \begin{pmatrix} L_{\infty j} \\ g_j \\ t_j + 10 \end{pmatrix} \sim \text{MVN}(\mu, \Sigma)$$

$$\mu = \log(\bar{L}_{\infty}, \bar{g}, \bar{t})$$

Where FUN is the length-at-age function for growth form x (Table 1), 
*L<sub>i,j</sub>* and  *age<sub>i,j</sub>* are length and age of fish *i* at 
sampling event *j*, *L<sub>∞j</sub>*, *g<sub>j</sub>*, and *t<sub>j</sub>* are 
respectively the asymptote, scaling (g<sub>1-3,j</sub>), and inflection 
(t<sub>0,j</sub> or t<sub>inf,j</sub>) parameters, and ε<sub>i,j</sub> are the 
Student’s t distributed errors with ν degrees of freedom. μ contains the 
grand means of the global parameters and Σ is a covariance-variance matrix. 
To aid in model convergence, growth parameters were estimated on the natural 
log scale, with the addition of 10 to the inflection parameter to allow for 
negative values. For the von Bertalanffy model, the inflection parameter 
t<sub>0</sub> was often close to -10 and biased by the 
addition of 10, so we left t<sub>0</sub> on the untransformed scale.  


**Categorical predictors (categorical)**

Same structure as random effect model but:
$$\mu = \log \begin{pmatrix} \bar{L}_{\infty} + \gamma_{1,h} \times \text{hydr}_j \\ \bar{g} + \gamma_{2,h} \times \text{hydr}_j \\ \bar{t} + \gamma_{3,h} \times \text{hydr}_j \end{pmatrix}$$

Where γ<sub>1-3,h</sub> are the fixed effects coefficients for hydroperiod 
classification *h* (1 = short, 2 = intermediate, 3 = long) of sampling event *j*
on the three growth parameters. 

**Continuous predictors (continuous)**

Same structure as random effect model but:
$$\mu = \log \begin{pmatrix} \bar{L}_{\infty} + \gamma_{1,k} \times \text{PC}_{k,j} \\ \bar{g} + \gamma_{2,k} \times \text{PC}_{k,j} \\ \bar{t} + \gamma_{3,k} \times \text{PC}_{k,j} \end{pmatrix}$$

Where γ<sub>1-3,h</sub> are the fixed effects coefficients for environmental 
PC *k* on the three growth parameters.

### Model development
**Script:** ```BLANK.R```

## Figures
**Directory:** ```figures/```

Contains raw R plots and .csv tables used to create the figures and tables in 
the main manuscript and appendices. Most figures were edited in Adobe Illustrator
for aesthetic purposes, and we included Illustrator files as well.

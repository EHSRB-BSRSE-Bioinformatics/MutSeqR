# BMD modeling using PROAST

Calculate the benchmark dose (BMD) of continuous, individual-level data
with optional model averaging. This function is intended to model the
dose-response of mutation frequency. This function is an extension of
the PROAST software (copyright RIVM National Institute for Public Health
and the Environment).

## Usage

``` r
bmd_proast(
  mf_data,
  dose_col = "dose",
  response_col = "mf_min",
  covariate_col = NULL,
  bmr = 0.5,
  adjust_bmr_to_group_sd = FALSE,
  model_averaging = TRUE,
  num_bootstraps = 200,
  plot_results = FALSE,
  output_path = NULL,
  raw_results = FALSE,
  seed = 125
)
```

## Arguments

- mf_data:

  A data frame containing the data to be analyzed. Data should be
  individual for each sample. Required columns are the column containing
  the dose `dose_col` the column(s) containing the mutation frequency
  `response_col`, and the column containing the covariate
  `covariate_col`, if applicable.

- dose_col:

  The column in `mf_data` containing the dose data. Values must be
  numeric. Default is "dose".

- response_col:

  The column(s) in `mf_data` containing the mutation frequency. Multiple
  `response_col`s can be provided. Default is "mf_min".

- covariate_col:

  The column in `mf_data` containing the covariate. If no covariate is
  present, set to `NULL` (default).

- bmr:

  The Benchmark Response value. The BMR is defined as a bmr-percent
  change in mean response relative to the controls. Default is 0.5 (50%
  change).

- adjust_bmr_to_group_sd:

  A logical value indicating whether the group standard deviation should
  be used as the BMR. If TRUE, the BMR will be bet set to one standard
  deviation above the control group mean. Default is FALSE.

- model_averaging:

  A logical value indicating whether confidence intervals should be
  calculated using model averaging. Default is TRUE (recommended).

- num_bootstraps:

  The number of bootstrap resamples to be used in the model averaging.
  Default is 200 (recommended).

- plot_results:

  A logical value indicating whether to plot the BMD models and/or the
  Cleveland plots. Default is FALSE. Plots may be exported directly to
  an output_path, or returned within a list to the user.

- output_path:

  The file path indicating where to save the plots. If NULL, the plots
  will automatically be displayed to the graphics window and then
  returned as a list alongside the bmd results.

- raw_results:

  A logical value indicating whether to return the raw results from the
  PROAST analysis. If FALSE, data is returned as a summary table.

- seed:

  An integer value to set the random seed for reproducibility when using
  model averaging. Default is 125.

## Value

A summmary data frame of final results. If plots or raw results are
selected, all data will be returned within a list.

The summary will include the following for each response variable and
covariate subgroup (if applicable):

- `Model`: The m3 or m5 model selected for each model family
  (Exponential, Hill, Inverse Exponential, LogNormal).

- `Response`: The response variable.

- `Covariate`: The covariate subgroup, if applicable.

- `bmr`: The specified Benchmark Response.

- `BMD`: The Benchmark Dose, in original dose units, estimated for the
  given model.

- `BMDL`: The lower bound of the 90% confidence interval for the BMD,
  calculated by the profile likelihood method.

- `BMDU`: The upper bound of the 90% confidence interval for the BMD,
  calculated by the profile likelihood method.

- `AIC`: The Akaike Information Criterion for the selected model. Lower
  values indicate a better fit. It is advised to choose the BMD value
  from the model with the lowest AIC.

- `weights`: The weight of the model in the model averaging process, if
  applicable.

- `Model averaging`: The BMDL and BMDU calculated by the bootstrap
  method if `model_averaging = TRUE`.

If there is no significant response in the data, the function will
return an empty data frame.

If `plot_results = TRUE` the function will create the following plots
for each response variable. The plots will be saved to the output_path.
If no output_path is provided, then they will be returned within a list
alongside the summary data frame.

- Model Plots. The following plot will be created for each model family
  (Exponential, Hill, Inverse Exponential, LogNormal): The fitted curve
  of the selected (3 or 5) model. Data is log-transformed. Individual
  data points are plotted using small triangles. The geometric mean
  (median) at each dose is plotted as a large triangle. The BMD is
  indicated by the dotted line. If applicable, the covariate subgroup is
  indicated by color. When output_path = NULL, plots are returned
  alongside results as recordedPlot objects. View plots with
  replayPlot().

- bootstrap_curves If `model_averaging = TRUE`, the bootstrap curves
  based on model averaging. The geometric mean (median) at each dose is
  plotted as a large triangle. Data is log-transformed. When output_path
  = NULL, plots are returned alongside results as recordedPlot objects.
  View plots with replayPlot().

- cleveland plot if `model_averaging = TRUE` The BMD estimate for each
  model is plotted as a red point alongside the 90% confidence
  intervals. The size of the BMD point represents the model weight
  assigned during model averaging, based on the AIC. When output_path =
  NULL, plots are returned alongside results as ggplots. View plots by
  calling object name.

If `raw_results = TRUE`, the function will return the raw results of the
PROAST analysis alongside the summary data frame. PROAST raw_results is
a list of variables and data that is continuously modified as it is
passed through the proast functions. It can be given to f.proast() to
resume analysis.

## Details

This function is a modified version of the original interactive PROAST
software (<https://www.rivm.nl/en/proast> that allows for batch
processing of data. The function is designed to be used with the output
of `calculate_mf` for the purpose of calculating the Benchmark Dose of
mutation frequency data. As such, some functionality of the original
PROAST software has been removed.

This function will accept continuous data, with an observation for each
individual subject. It is assumed that data are lognormally distributed.
The response data is log-transformed, then back-transformed after the
statistical analysis. The function will fit model 3 or 5 from various
families of models (Exponential, Hill, Inverse Exponential, LogNormal).
It will then compare the fits of models 3 and 5 for each model family
and select the model with the lowest AIC. The BMD 90% confidence
intervals will be calculated based on the selected model (3 or 5) for
each model family using the profile likelihood method. The BMD 90%
confidence interval may also be calculated using the bootstrap method if
`model_averaging = TRUE`. It is recommended to use 200 bootstraps for
model averaging.

To replicate these results in the PROAST interactive software, select
the following menu options:

1.  f.proast(mf_data)

2.  What type of response data do you want to consider? *1: continuous,
    individual data*

3.  Do you want to fit a single model or fit various nested families of
    models? *3: select model 3 or 5 from various families of models*

4.  Q1: Which variable do you want to consider as the independent
    variable? *\# : dose_col*

5.  Give number(s) of the response(s) you want to analyse. *\# :
    response_col*

6.  Give number of factor serving as potential covariate (e.g.sex) type
    0 if none. *\# : covariate_col*

7.  Do you want to adjust CES to within group SD? *1: no, 2: yes \|
    adjust_bmr_to_group_sd: FALSE/TRUE*

8.  Give value for CES (always positive) type 0 to avoid calculation of
    CIs. *bmr*

9.  Do you want to calculate the BMD confidence interval by model
    averaging? *1: no 2: yes \| model_averaging: FALSE/TRUE*

10. give number of bootstrap runs for calculating BMD confidence
    interval based on MA (e.g. 200) *num_bootstraps*

11. Which models do you want to be fitted? *4 : previous option with
    lognormal DR model added*

## Examples

``` r
# Example data consists of 24 mouse bone marrow DNA samples imported
# using import_mut_data() and filtered with filter_mut().
# Data was summarized per sample using calculate_mf() (see relevant
# examples).
mf_example <- readRDS(system.file("extdata/Example_files/mf_data_global.rds",
  package = "MutSeqR"
))
# We will calculate the BMD for a 50% increase in mutation frequency from
# control with Model averaging.
# For the purpose of this example, num_bootstraps is set to 3 to reduce
# run time. 200 bootstraps is recommended.
bmd <- bmd_proast(
  mf_data = mf_example,
  dose_col = "dose",
  response_col = "mf_min",
  bmr = 0.5,
  model_averaging = TRUE,
  num_bootstraps = 3
)
#> Independent variable: dose
#> Running in batch/command line mode...
#> Independent variable: dose
#> Independent variable: dose
#> Warning: Coercing LHS to a list
#> Setting data type to continuous, individual data
#> Selected model option: 'select model 3 or 5 from various families of models' (3)
#> indep_var_choice: dose
#> dose
#> Independent variable:
#> dose
#> moving on...
#> Using dose as the independent variable.
#> 
#> 1 sample 
#>  2 dose 
#>  3 dose_group 
#>  4 sum_min 
#>  5 sum_max 
#>  6 mf_min 
#>  7 mf_max 
#>  8 group_depth 
#> Give number(s) of the response(s) you want to analyse
#> Model averaging enabled.
#> ans.all$CES: 0.5
#> 
#> 
#> response:  mf_min
#> ANALYSIS WITH EXPONENTIAL MODELS
#> 
#> model     converged   npar    loglik      aic
#>  full model   1   5   14.23   -18.46
#>  null model   1   2   -24.44      52.88 
#>  Expon. m3-   1   4   8.23    -8.46
#> warning: lower constraint on parameter cc has been increased to 1.515 
#> due to value of CES
#> This might result in suboptimal fit of model 15 
#>  
#>  Expon. m5-   1   5   14.23   -18.46
#>  ------------------------------------------------------- 
#> selected model: Expon. m5-  
#> 
#>  estimate for  var- :  0.01789 
#>   estimate for  a- :  1.718e-07 
#>   estimate for  CED- :  8.45 
#>   estimate for  c- :  5.611 
#>   estimate for  d- :  1.629 
#> ------------------------------------------------------ 
#> 
#>  
#> The Grubb outlier for this sample size is: 2.80155 
#>  
#> No outliers detected 
#> 
#> 
#>  calculating confidence intervals ....
#> 
#> 
#> 
#> the CED (in orig. units) and the 90 % confidence interval is: 
#>  8.45 
#>  6.47 
#>  10.3 
#> 
#> 
#> 
#> response:  mf_min 
#> ANALYSIS WITH HILL MODELS
#> 
#> model     converged   npar    loglik      aic 
#>  Hill m3-     1   4   8.28    -8.56
#> warning: lower constraint on parameter cc has been increased to 1.515 
#> due to value of CES
#> This might result in suboptimal fit of model 25 
#>  
#>  Hill m5-     1   5   14.23   -18.46
#>  ------------------------------------------------------- 
#> selected model: Hill m5-  
#> 
#>  estimate for  var- :  0.01789 
#>   estimate for  a- :  1.718e-07 
#>   estimate for  CED- :  9.03 
#>   estimate for  c- :  6.236 
#>   estimate for  d- :  2.292 
#> ------------------------------------------------------ 
#> 
#>  
#> The Grubb outlier for this sample size is: 2.80155 
#>  
#> No outliers detected 
#> 
#> 
#>  calculating confidence intervals ....
#> 
#> 
#> 
#> the CED (in orig. units) and the 90 % confidence interval is: 
#>  9.03 
#>  7.02 
#>  10.7 
#> 
#> 
#> 
#> response:  mf_min 
#> ANALYSIS WITH INVERSE EXPONENTIAL MODELS
#> model     converged   npar    loglik      aic 
#>  Inv.Expon. m3-   1   4   9.09    -10.18
#> warning: lower constraint on parameter cc has been increased to 1.515 
#> due to value of CES
#> This might result in suboptimal fit of model 52 
#>  
#>  Inv.Expon. m5-   0   5   -8.86   27.72
#>  ------------------------------------------------------- 
#> selected model: Inv.Expon. m3-  
#> 
#>  estimate for  var- :  0.02745 
#>   estimate for  a- :  1.701e-07 
#>   estimate for  CED- :  4.714 
#>   estimate for  d- :  0.1632 
#> ------------------------------------------------------ 
#> 
#>  
#> The Grubb outlier for this sample size is: 2.80155 
#>  
#>  1 outliers detected, with values for x and y: 
#> 25 9.25584825420905e-07  
#> 
#> 
#> 
#>  calculating confidence intervals ....
#> 
#> 
#> 
#> the CED (in orig. units) and the 90 % confidence interval is: 
#>  4.714 
#>  2.9 
#>  6.96 
#> 
#> 
#> 
#> response:  mf_min 
#> ANALYSIS WITH LOGNORMAL DR MODELS
#> model     converged   npar    loglik      aic 
#>  LN m3-   1   4   8.72    -9.44
#> warning: lower constraint on parameter cc has been increased to 1.515 
#> due to value of CES
#> This might result in suboptimal fit of model 54 
#>  
#>  LN m5-   1   5   14.23   -18.46
#>  ------------------------------------------------------- 
#> selected model: LN m5-  
#> 
#>  estimate for  var- :  0.01789 
#>   estimate for  a- :  1.718e-07 
#>   estimate for  CED- :  9.171 
#>   estimate for  c- :  5.956 
#>   estimate for  d- :  1.467 
#> ------------------------------------------------------ 
#> 
#>  
#> The Grubb outlier for this sample size is: 2.80155 
#>  
#> No outliers detected 
#> 
#> 
#>  calculating confidence intervals ....
#> 
#> 
#> 
#> the CED (in orig. units) and the 90 % confidence interval is: 
#>  9.171 
#>  7.27 
#>  10.7 
#> 
#> 
#> 
#>  -----------  CES =  0.5  --------------------------------- 
#> 
#> Calculating confidence intervals by model averaging, this may make some time ....
#> 
#> 
#>  The weights used in model averaging are:
#>    model weight
#> 1    EXP 0.3316
#> 2   HILL 0.3316
#> 3 INVEXP 0.0053
#> 4   LOGN 0.3316
#> 
#> Start of MA bootstrap runs ...
#> run   1  2  3
#> 
#> duration of bootstrap calculations:
#> [1] "Wed Nov 19 22:24:57 2025"
#> [1] "Wed Nov 19 22:24:57 2025"
#> 
#> The model-average BMD confidence interval is:
#>   subgroup BMDlower.ma BMDupper.ma
#> 1      all        8.94          10
#> 
#> 
#> 
#> end of analysis for response:  mf_min
#> LN m5-
#> LN m5-
#> Covariate analysis is:FALSE
#> Hill m5-
#> Hill m5-
#> Covariate analysis is:FALSE
#> Inv.Expon. m3-
#> Inv.Expon. m3-
#> Covariate analysis is:FALSE
#> Expon. m5-
#> Expon. m5-
#> Covariate analysis is:FALSE
#> model_averaging
#> Expon. m5-
#> Covariate analysis is:FALSE
#> Handling model averaging case
#>    Selected.Model Response CES   CED CEDL CEDU    AIC Log.Likelihood
#> 1      Expon. m5-   mf_min 0.5  8.45 6.47 10.3 -18.46          14.23
#> 2        Hill m5-   mf_min 0.5  9.03 7.02 10.7 -18.46          14.23
#> 3  Inv.Expon. m3-   mf_min 0.5 4.714  2.9 6.96 -10.18           9.09
#> 4          LN m5-   mf_min 0.5 9.171 7.27 10.7 -18.46          14.23
#> 5 Model averaging   mf_min 0.5 9.065 8.94   10    N/A            N/A
#>                  Var                    a                 d weights
#> 1  0.017888312831692 1.71849444527397e-07  5.61135528663567  0.3316
#> 2 0.0178883278572467 1.71849441477975e-07  6.23608032503139  0.3316
#> 3 0.0274547356721535 1.70137658821088e-07 0.163165375364757  0.0053
#> 4 0.0178883115661859 1.71849497597811e-07  5.95609683921669  0.3316
#> 5                N/A                  N/A               N/A      NA
```

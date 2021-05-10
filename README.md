# Pre-requisites

## vwmvp package

Model fitting, simulation and analysis of model comparison results using the files in this repo require the vwmvp package to be loaded (tar.gz.file here, or see the separate repo [[link](https://github.com/nklange/vwmvp)]. 

## Naming conventions 

Importantly, the package and the analysis files use different shorthand to refer to the models than the manuscript. Consequently, all output from the model fitting in various .rds files and code and further analyses refer to the package shorthand and only re-name to the manuscript names when when finalizing figures and tables. For reference, for the main models we discuss:

Manuscript  | Analysis
------------- | ------------- 
VP(J)A-    | J_RNminus      
VP(J)A+        | J_RNplus
VP(&kappa;)A- | MK_RNminus
VP(&kappa;)A+ | MK_RNplus
VP(&kappa;)F- | MK_FM_RNminus
VP(&kappa;)F+ | MK_FM_RNplus
VP(&kappa;)P- | MK_P_RNminus
VP(&kappa;)P+ | MK_P_RNplus
VP(&kappa;)U- | MK_U_RNminus
VP(&kappa;)U+ | MK_U_RNplus

Additionally, we fitted the VP(&kappa;)A&plusmn; model separately to all set sizes, with &kappa;, &tau; and &kappa;<sub>r</sub (where applicable) were estimated separately for all set sizes and $\alpha$ was ommitted. In the manuscript, these models are identified as **VP, set sizes fitted separately**. In the fitting and prediction files, these models are identified as **VPnosetsize** and **VPplusnosetsize** for the model without and with response noise respectively.

# Data

## Individual data sets

The data has been preprocessed to be in correct format from original sources. No changes were made beyond that.

* **Data/ol17_prepared.rda**: ol17_e1 experiments 1, 2 and 3. Only Experiment 1 is included in the data corpus
* **Data/vdb14_prepared.rda**: data corpus from van den Berg et al. (2014)
* **Data/pratte17_prepared.rda**: pratte17

For the fits of the aggregate data, individual data sets were aggregated within experiments.

## Preparing crossvalidation

**MakeCvSets.R**: splitting individual data sets into separate folds for 10-Fold CV, 5 x 2-Fold CV and LOSsO-CV.

In **CVData/**

Training data  | Test data
------------- | ------------- 
trainforCV10Fold.rds  | testforCV10Fold.rds      
trainforCV5times2.rds       | testforCV5times2.rds
trainforCVLOSZO.rds | testforCVLOSZO.rds

For aggregate fits in LOSSO-CV, the individual-level files were combined.


Out-of-sample predictions for hold-out sets are collected in **Fits/**, based on prediction using best fit parameter estimates for the training sets in **CVPrediction.R**.

# Model fitting

All models were fitted using functions collected in the **vwmvp** package.

## Fitting procedure

The model fitting routine for the fits reported in the manuscript is shown in **ModelFittingRoutine.R**. This uses the top-level fitting function from vwmvp (vwmvp::FitVP). All fits reported in the manuscript were based on setting rep = 20, method = "numint" and startpar = NULL, where starting parameters and data preparation wrapped by the top-level function. For more details on the routine and model implementation, see the vwmvp repo [[link](https://github.com/nklange/vwmvp)].

## Fitting results

All individual and aggregate level fits (all runs for each model and data set) are in files in the **Fits/** folder.

### Individual-level fits

* **FitFull.rds**: fit for comparison by AIC, numerical integration
* **CV10Fold_Trainingfits.rds**: CV training fits, 10 F CV
* **CV5x2Fold_Trainingfits.rds**: CV training fits, 5x2 F CV
* **CVLOSzO_Trainingfits.rds**: CV training fits, LOSsO CV

* **CV10Fold_TestLL.rds**: Hold-out set fits from best fitting training run, 10 F CV
* **CV5x2Fold_TestLL.rds**: Hold-out set fits from best fitting training run, 5x2 F CV
* **CVLOSzO_TestLL.rds**: Hold-out set fits from best fitting training run, LOSsO CV

* **FitFull_sepSS.rds**: VPnosetsize and VPplusnosetsize
* **FitFull_simSS.rds**: negative log likelihood for individual set sizes on the basis of the best fit parameters in **FitFull.rds**

* **GA_FitFull.rds**: Simulation approach in fitting VP models, with GA, for VP(J)A-, VP(J)A+, VP(&kappa;)A- and VP(&kappa;)A+
* **GAnlminb_FitFull.rds**: Simulation + numerical integration for VP(J)A-, VP(J)A+, VP(&kappa;)A- and VP(&kappa;)A+

### Experiment-level (aggregated data) fits

* **FitFull_agg.rds**: fit for comparison by AIC, numerical integration
* **FitFull_AggLOSSO_Training.rds**: LOSSO-CV Training fits for VP models
* **FitFull_AggLOSSO_TestLL.rds**: LOSSO-CV hold-out set fits for VP models)

* **FitFull_Agg_sepSS.rds**: VPnosetsize and VPplusnosetsize
* **FitFull_Agg_simSS.rds**: negative log likelihood for individual set sizes on the basis of the best fit parameters in **FitFull_agg.rds**

# Predictions of behavioral signature pattern

The routine for the predictions of error distributions and summary statistics is in **PredictBehavioralSignature.R**.

The results are stored in **Prediction/** and **SummaryStat/** respectively. While we largely show aggregate experiment-level predictions in the manuscript, we show individual, and averaged individual, predictions in the individual-level graphs available in **Individual Graphs/**

The name of the files, and columns in the files, indicate the content / level of analysis.

* individual: individual best-fit parameters
* id = "Av": averaging individual best-fit parameters
* id = "AvExcl": averaging individual best-fit parameters, excluding participants for all models in the file that show extreme LOSSO-CV out-of-sample deviance for at least one model
*aggregate: best-fit parameters of experiment-level aggregate fit

* "leftout"" column: 1, 2, 3, 4, 5, 6, 7, 8 indicate best-fit parameters were based on fits to data sets without that set size; 0 indicates that the full data set was fitted (which may not include all set sizes by design). For the summary statistics files, 9 indicates observed data.


**Prediction/** for error distributions. [experiment] indicates placeholder for experiment name.

File  | Models | Level | LOSSO-CV
------------- | ------------- |------------- |------------- 
prediction_VP_fullLOSSO_[experiment]_ind.rds | VP | Individual, id = "Av" | Full data set, LOSSO
prediction_nonVP_full_[experiment]_ind.rds | non-VP models | Individual, id = "Av" | Full data set
prediction_idavExclFailed_[experiment]_ind.rds | VP models in RQ2 (no Fisher information) | id = "AvExcl" | Full data set, LOSSO
prediction_sepSS_full_[experiment]_ind.rds       | VPnosetsize,VPplusnosetsize | Individual, id="Av" | Full data set
prediction_VP_fullLOSSO_[experiment]_agg.rds | VP | Individual, id = "Av" | Full data set, LOSSO for VP(&kappa;) models
prediction_nonVP_full_[experiment]_agg.rds | non-VP models | Individual, id = "Av" | Full data set
prediction_sepSS_full_[experiment]_agg.rds       | VPnosetsize,VPplusnosetsize | Aggregate | Full data set

For the figures in the manuscript, resultant summary statistics are calculated in-line in the relevant file (i.e., subject to some sampling error). For the experiment/individual-level files, they are calculated using the routine in the file.

# Manuscript figures/analyses

Figure and Analysis  | File
------------- | ------------- 
Figure 1: Types of model comparison  | **ModelComparisonFigures.R**
Figure 2: Introduction Case Study | **Illustration_VMVP.R**
Figure 3: Theoretical Predictions RQ1 | **RQ1_TheoPred.R**
Speed of fitting | **DiagnoseModels.R**
Figure 4: RQ1 Results | **RNFisher_graph.R** (calls **CompareModels_FisherRN.R** for quantitative fits)
Figure 5: Theoretical Predictions RQ2 | **RQ2_TheoPred.R**
Figure 6: RQ2 Results | **RNCapacity_graph.R**  (calls **CompareModels_LimitRN.R** for quantitative fits)
Figure 7: Set size 1 | **Setsize1_paper.R**
Figure 8: Separate fits to set sizes | **Setsize1_paper.R**

# Supplemental Material figures/analyses

Figure/Analysis  | File
------------- | ------------- 
Figure 1: Additional figures for RQ1  | **CompareModels_FisherRN.R**
Figure 2: Additional figures for RQ2 | **CompareModels_LimitRN.R**
Table 1: Parameter estimates | **DiagnoseModels.R** (*)
Figure 3: Correlation parameters and stability | **DiagnoseModels.R**
Figure 4: Quantitative fits Non-VP models | **NonVP_graph.R**  (calls **CompareModels_NonVP.R**)
Figure 5: Normalized RMSD of Non-VP models | **NonVP_graph.R**  (calls **CompareModels_NonVP.R**)
Figure 6: Shen \& Ma (2019) Factorial Importance Metrics RQ1 |  **RNFisher_graph.R**
Figure 7: Shen \& Ma (2019) Factorial Importance Metrics RQ2 |  **RNCapacity_graph.R**
Figure 8: Shen \& Ma (2019) Factorial Importance Metrics RQ2 as interaction |  **NonVP_graph.R**
Figure 9: Alternate VP variants for limited capacity | **Compare_modelvariants.R**

(*) Note: in vwmvp, the uniform capacity limit is defined as running from 0 - K (rather than 0 - 2K, as defined in the manuscript). Parameter estimates of K in VP(&kappa;)&plusmn; models in fitting files therefore were divided by 2 to reflect the the mean K rather than the max(K). 

## Parameter and model recovery

We ran model/parameter recovery for a limited set of generating parameters. We used two approaches to determine generating parameters: choosing median parameter estimates, and sampling parameters from the multivariate parameter space given by the range of parameter estimates. 

### Generating data sets

In contrast to the predictions of the error distributions which we generated using vwmvp::predict_data(), we generated data sets here using vwmvp::generate_data(). As we show in **vwmvp_predictvsgenerate.R**, given a suitably large sample size (e.g., number of trials), both methods are approximately equivalent. When we initially generated the data sets for recovery, vwmvp::generate_data() and vwmvp::predict_data() both defined the U-capacity limit to run from 0 - $K$ with a mean of $K$/2. In the current version, vwmvp::generate_data() defines it with 0 - 2$K$ with a mean of $K$ (hence the need to adjust the input-$K$ by multiplying it by 2 to generate equivalent predictions in **vwmvp_predictvsgenerate.R**).

The generation of data sets is shown in **Recovery_generateparameters.R**. Detailed explanation for adjusting/determining the values is given in the supplemental material. The data sets are stored in **RecoveryData/** where **median** refers to data sets with median parameter estimates for generating values, and **rmv** referes to samples from a multivariate normal for generating values. The number in the file name refers to the number of trials per set size that were generated.

### Analyzing Recovery

Figure/Analysis  | File
------------- | ------------- 
Table 2: generating parameter values for rmv recovery | **Recovery_generateparameters.R**
Figure 10: Parameter recovery |  **Recovery_analyzePR.R** (calls **Recovery_analyzeMR.R**)
Figure 11: Model recovery |  **Recovery_analyzePR.R** (calls **Recovery_analyzeMR.R**)

# Experiment/Individual-level figures

Experiment-level and individual-level graphs/analyses matching the ones in the manuscript are shown in the .html files in the **IndividualGraphs/** folder. They are based on running the .Rmd files **RQ1RQ2_ExperimentIndividualLevel.Rmd** and **RQ3_ExperimentIndividualLevel.Rmd**.

**IndividualGraphs/SummaryStat/** contains the summary statistics used in these files, **Prediction/** contains the prediction of the error distribution as detailed above. [experiment] indicates placeholder for experiment name.

File  | Models | Level | LOSSO-CV
------------- | ------------- |------------- |------------- 
SumStat_altmod_[experiment]_ind.rds | non-VP | Individual, id= "Av" | full data set
SumStat_altmode_[experiment]_agg.rds | non-VP | Aggregate | full data set
SumStatLOSSO_[experiment]_agg.rds | VP | Aggregate | Full data set (*), LOSSO
SumStat_[experiment]_agg.rds | VP | Aggregate | Full data set
SumStatLOSSO_[experiment]_ind.rds | VP | Individual, id = "Av" | Full data set (*), LOSSO
SumStat_[experiment]_ind.rds | VP | Individual, id="Av" | Full data set
SumStatLOSSO[experiment]_averagepred.rds | VP | summary statistics based on averaging individual error distribution predictions with predictions based on extreme out-of-sample LOSSO-CV removed | Full data set (*), LOSSO
SumStat_[experiment]_averagepred.rds | VP | summary statistics based on averaging individual error distribution predictions
with predictions based on extreme out-of-sample LOSSO-CV removed | Full data set, LOSSO
SumStatLOSSO_[experiment]_idAvExcl.rds | VP | id = "AvExcl" | Full data set, LOSSO

(*) The 'full data set' is exchangeable in both files. Small numerical differences due to simulation noise but based on the same data. In the LOSSO files, VP(J)A&plusmn; models are missing as we do not show the LOSSO-CV predictions of those models at any point and it was not relevant to include them, though simulating the predictions would be trivial.






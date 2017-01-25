# Vortex_OCD

This is the data and code required to obtain the analysis results presented in Vaghi, Luyckx, et al., XXX.Disentangling action from confidence in Obsessive-Compulsive Disorder

The code is written by Matilde Vaghi and Fabrice Luyckx

- The code is a mixture of Matlab and R. 

- Analysis
  	- Vortex_basicstats.m
 	- Vortex_bic.m
 	- Vortex_changePoint.m 
  	- Vortex_ErrorMagnitude.m
  	- Vortex_load.m
  	- Vortex_medication.m 
 	- Vortex_PlotErr.R
 	- Vortex_Regressions.m 

- Data 
  	- Behavioural 
  		- Vortex_ocd_fulldata.mat 

 	- Questionnaires
		- VortexBet_DepOCDMeasures.txt 

	 - RPlotData
		- VortexBet_ErrorMag.txt
		- VortexBet_LR.txt
		- VortexBet_RegConfidence.txt
		- VortexBet_RegDiscrepancy.txt
		- VortexBet_RegLearningRate.txt
  
- Experiment
	- Vortex_main.m 
	- Vortex_functions
		- 

- Functions 
	- diffcirc.m
	- estimateLR.m
	- getChange.m 
	- redBayes_circ.m
	- regr_fitglm.m
	- round.m
	- shift.m

- Plots
	- Plot_Correlation.tiff
	- Plot_ErrorMagnitude.tiff
	- Plot_LearningRate.tiff
	- Plot_RegConfidence.tiff
	- Plot_RegDiscrepancy.tiff
	- Plot_RegLearning.tiff
	- Plot_ModelBehaviour.eps

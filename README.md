# Vaghi, Luyckx et al (2017)

This is the data and code required to obtain the analysis results presented in Vaghi, Luyckx, et al., 2017. Compulsivity reveals a novel dissociation between action and confidence.

The code is written by Matilde Vaghi and Fabrice Luyckx

- The code is a mixture of Matlab and R. 

- Analysis
  	- Vortex_basicstats.m
	- Vortex_basicstats_R1Plots.m
	- Vortex_basicstats_R1Trajectories.m
	- Vortex_basicstats_testInteraction.m
 	- Vortex_bms.m
 	- Vortex_changePoint.m 
	- Vortex_changePoint_R1Subjs.m
	- Vortex_changePoint_testInteraction.m
  	- Vortex_ErrorMagnitude.m
	- Vortex_ErrorMagnitude_R1.m
	- Vortex_ErrorMagnitude_R1_add.m
	- Vortex_hazardFit.m
  	- Vortex_load.m
	- Vortex_lr_med_qst.m
  	- Vortex_medication.m 
 	- Vortex_PlotErre.R
	- Vortex_PlotErretestinteraction.R
 	- Vortex_regressions.m 
	- Vortex_variables.m

- Data 
  	- Behavioural 
  		- Vortex_ocd_fulldata.mat 

 	- Questionnaires
		- VortexBet_DepOCDMeasures.txt 

	 - RPlotData
		- ErrorMagn_20bins.txt
		- ErrorMagn_3bins.txt
		- ErrorMagn_3bins_add.txt
		- VortexBet_ChangePoint.txt
		- VortexBet_ChangePointModel.txt
		- VortexBet_ChangePointSubjs.txt
		- VortexBet_ErrorMag.txt
		- VortexBet_LR.txt
		- VortexBet_LR_CONF_ZSCORES.txt
		- VortexBet_LR_CONF_ZSCORE_TIME.txt
		- VortexBet_RegConfidence.txt
		- VortexBet_RegDiscrepancy.txt
		- VortexBet_RegLearningRate.txt
		- VortexBet_ZVALUES_ChangePoint
		- VortexBet_ZVALUES_LRCF_CPP.txt
  
- Experiment
	- Vortex_main.m 
	- Vortex_functions
		- Instructions
		- Vortex_initialise.m
		- Vortex_instructions.m
		- Vortex_send_email.m
		- Vortex_setup.m
		- Vortex_trial.m
		- circle.m
		- printText.m
		
- Functions 
	- diffcirc.m
	- estimateLR.m
	- getChange.m 
	- redBayes_circ.m
	- regr_fitglm.m
	- round.m
	- shift.m

- Plots
	- BMS_confidence.eps
	- BMS_update.eps
	- Figure1.tiff
	- Figure2.tiff
	- Figure3.tiff
	- Figure4.tiff
	- Hazard_fits.eps
	- Plot_ModelBehaviourAgainstSubject_R1.tif

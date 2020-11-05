README.txt file for Parametric versus rank-based, non-parametric methods 
for differential expression analysis in single-cell RNAseq using a synthetic dataset 
Created 08/2020

The most important document in this appendix is the Guided_Tutorial_of_Analysis 
provided in either .html or .pdf format. This document provides a highly detailed
step-by-step documentation of how this analysis was carried out from the evaluation 
of the data simulators through to the differential expression analysis using a sample
size of 50 cells per condition as an example. Details of various functionalities are 
included as well as the flags that were chosen in the analysis so reading this document 
is highly recommended.

NOTE: Included are also the data simulation R scripts for the other datasets however 
the code is only a repeat of what is outlined in the guided tutorial with necessary 
parameters altered such as names of objects as per the dataset. Also included are the
differential analysis R scripts for the additional sample sizes however similarly these
are repeats of the 50 cell per condition pipeline in the Guided Tutorial with only the 
number of cells parameter changed.  

Also found in this appendix are additional figures that were omitted from the final paper:

- In the Simulator_Analysis_Figures folder you will find the plots for the 
  additional 6 10X datasets that were used to evaluate the simulator package performance

- In the Data_Simulation_Scripts you can also find the R script that was used to produce
  the simulator analysis results plots

- In the DE_Analysis_Figures folder you will find the Excel spreadsheet used to
  calculate all the additional statistic metrics: TP, TN, FP, FN, FDR, TPR, F1 
  and Precision

Happy reading!



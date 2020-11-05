README.txt file for Parametric versus rank-based, non-parametric methods 
for differential expression analysis in single-cell RNAseq using a synthetic dataset 


The most important document in this appendix is the Guided_Tutorial_of_Analysis 
provided in .pdf format. This document provides a highly detailed
step-by-step documentation of how this analysis was carried out from the evaluation 
of the data simulators through to the differential expression analysis using a sample
size of 50 cells per condition as an example. Details of various functionalities are 
included as well as the flags that were chosen in the analysis so reading this document 
is highly recommended.

NOTE: Included is also the data simulation R scripts for the Brain10X dataset however 
the code is only a repeat of what is outlined in the guided tutorial with necessary 
parameters altered such as names of objects as per the dataset. Also included is the
differential analysis R script for the 100 cell per sample size however similarly these
are repeats of the 50 cell per condition pipeline in the Guided Tutorial with only the 
number of cells parameter changed.  

Also found in this repository are:

- My final dissertation paper to give context to the code given

- The script used to calculate all the additional statistic metrics: TP, TN, FP, FN, FDR, TPR, F1 
  and Precision

Happy reading!



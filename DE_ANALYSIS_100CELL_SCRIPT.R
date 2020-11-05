#Simulation of data

# Running scDesign simulation to obtain simulated dataset where we have 2 conditions with 100 cells in each one, there are a total of 10% total genes that are DE (5% UP, 5% DOWN)
scDesign_drosophila <- scDesign::design_data(realcount = drosophila10X, S = rep(1e7, 2), ncell = rep(100, 2), ngroup = 2, 
                                             pUp = 0.05, pDown = 0.05, ncores = 1)

saveRDS(scDesign_drosophila, file = "/datastore/2505621h/scDesign_simulated_dataset_FINAL")

# Merging the count tables from the two conditions into one master count matrix
drosophila <- cbind(scDesign_drosophila$count$count1, scDesign_drosophila$count$count2)

saveRDS(drosophila, file = "/datastore/2505621h/DE_Analysis/droso_object_FINAL")
drosophila <- readRDS(file = "/datastore/2505621h/DE_Analysis/droso_object_FINAL")

# Filter out genes that are all zeros
drosophila <- dataclean(drosophila)

# MAST DE Analysis 

coldata_MAST <- data.frame(
  condition = factor(c(
    rep("0", 100),
    rep("1", 100))))

row.names(coldata_MAST) <- colnames(drosophila)

eset_drosophila <- ExpressionSet(drosophila)
counts = exprs(eset_drosophila)
tpm <- counts*1e6/colSums(counts)
tpm <- log2(tpm+1)

saveRDS(tpm, file = "/datastore/2505621h/100cell_Analysis/drosophila_tpm")

sca <- FromMatrix(tpm,  cData = coldata_MAST)
ngeneson <- apply(exprs(eset_drosophila),2,function(x) mean(x>0))
CD <- colData(sca)
CD$ngeneson <- ngeneson
CD$cngeneson <- CD$ngeneson-mean(ngeneson)
colData(sca) <- CD
testing <- zlm(~ condition + cngeneson, sca = sca)
summaryCond_test <- summary(testing, doLRT='condition1') 
summaryDt <- summaryCond_test$datatable
fcHurdle <- merge(summaryDt[contrast=='condition1' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                  summaryDt[contrast=='condition1' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]

# Renaming & calculating score for compcodeR
fcHurdle$score <- 1-fcHurdle$`Pr(>Chisq)`
names(fcHurdle)[6] <- "FDR"
names(fcHurdle)[2] <- "pvalue"
fcHurdle <- fcHurdle[mixedorder(fcHurdle$primerid),]
fcHurdle <- tibble::column_to_rownames(fcHurdle, var = 'primerid')

saveRDS(fcHurdle, file = "/datastore/2505621h/100cell_Analysis/drosophila_MAST_results")

# SigEMD analysis
data <- dataclean(tpm)

condition_sigEMD <- c(
  rep("0", 100),
  rep("1", 100))

names(condition_sigEMD) <- colnames(data)

results_EMD <- calculate_single(data =  data,condition =  condition_sigEMD, Hur_gene = NULL, binSize=0.2,nperm=100)

emd<- as.data.frame(results_EMD$emdall)
names(emd)[3] <- "adjpvalue"
emd$score <- 1-emd$adjpvalue

saveRDS(emd, file = "/datastore/2505621h/100cell_Analysis/drosophila_SigEMD_results")

# DEsingle Anlaysis

# Setting up group matrix using 100 cells in each group
group <- factor(c(
  rep("1", 100),
  rep("2", 100)))


# Running DESingle analysis
deSingle_drosoph <- DEsingle(drosophila, group = group)

# Identifying significant genes that have a FDR padj value < 0.05
dros_result_groups <- DEtype(results = deSingle_drosoph, threshold = 0.05)
dros_result_groups$score <- 1-dros_result_groups$pvalue.adj.FDR
dros_result_groups <- tibble::rownames_to_column(dros_result_groups, var = "gene")

dros_result_groups <- dros_result_groups[mixedorder(dros_result_groups$gene),]
rownames(dros_result_groups) <- NULL
dros_result_groups <- tibble::column_to_rownames(dros_result_groups, var = 'gene')
names(dros_result_groups)[21] <- "adjpvalue"

saveRDS(dros_result_groups, file = "/datastore/2505621h/100cell_Analysis/deSingle_results_drosophila_FINAL")

# RankStat Analysis

drosophila_rlog <- rlogTransformation(drosophila_sc + 1)
drosophila_rlog <- as.data.frame(drosophila_rlog)

drosophila_rlog <- readRDS(file = "/datastore/2505621h/drosophila_rlog")
drosophila_rlog <- dataclean(drosophila_rlog)
gene_names <- row.names(drosophila)

coldata <- data.frame(
  condition = factor(c(
    rep("0", 100),
    rep("1", 100))))

coldata$condition <- relevel(coldata$condition, ref = "0")

coldata <- as.vector(t(coldata))


result_RS <- RankStat(drosophila_rlog, cl = coldata, method = c("RP", "DRP", "RS", "RD", "R_RD"), alternative = "two.sided", gene.names = gene_names)

names(result)
result_RP <- result_RS$RP_2s
result_RP$score <- 1-result_RP$p.value

saveRDS(result_RP, file = "/datastore/2505621h/100cell_Analysis/drosophila_result_RP_results")

result_DRP <- result_RS$DRP_2s
result_DRP$score <- 1-result_DRP$p.value

saveRDS(result_DRP, file = "/datastore/2505621h/100cell_Analysis/drosophila_result_DRP_results")

result_RSum <- result_RS$RS_2s
result_RSum$score <- 1-result_RSum$p.value

saveRDS(result_RSum, file = "/datastore/2505621h/100cell_Analysis/drosophila_result_RSum_results")

result_RD <- result_RS$RD_2s
result_RD$score <- 1-result_RD$p.value

saveRDS(result_RD, file = "/datastore/2505621h/100cell_Analysis/drosophila_result_RD_results")

result_R_RD <- result_RS$R_RD_2s
result_R_RD$score <- 1-result_R_RD$p.value
saveRDS(result_R_RD, file = "/datastore/2505621h/100cell_Analysis/drosophila_result_R_RD_results")

# Creation of sample annotation table showing which cells belong to each group to create compcodeR object
coldata_compcodeR <- data.frame(
  condition = factor(c(
    rep("0", 100),
    rep("1", 100))))

row.names(coldata_compcodeR) <- colnames(drosophila)

#Creation of parameters
parameters <- list(dataset = "Drosophila",
                   samples.per.cond = 100,
                   n.diffexp = 1760,
                   fraction.upregulated = 0.5,
                   uID = 1234567890)


# Creation of compcodeR object
comp_data_drosophila <- compcodeR::compData(drosophila, sample.annotations = coldata_compcodeR, info.parameters = parameters, variable.annotations = DE_State_annotation)

# Creation of pseudocount compcodeR object for DEseq2 analysis
drosophila_pseudo <- drosophila + 1
comp_data_drosophila_pseudo <- compcodeR::compData(drosophila_pseudo, sample.annotations = coldata_compcodeR, info.parameters = parameters, variable.annotations = DE_State_annotation)


saveRDS(comp_data_drosophila, file = "/datastore/2505621h/100cell_Analysis/comp_data_drosophila.rds")
saveRDS(comp_data_drosophila_pseudo, file = "/datastore/2505621h/100cell_Analysis/comp_data_drosophila_pseudo.rds")


# DEseq2 analysis
comp_code_DEseq2 <- compcodeR::runDiffExp(data.file = "/datastore/2505621h/DE_Analysis/100cell_Analysis/comp_data_drosophila_pseudo.rds", result.extent = "DESeq2",
                               Rmdfunction = "DESeq2.createRmd",
                               output.directory = "/datastore/2505621h/DE_Analysis/100cell_Analysis/", fit.type = "parametric",
                               test = "Wald", beta.prior = TRUE,
                               independent.filtering = TRUE, cooks.cutoff = TRUE,
                               impute.outliers = TRUE)

comp_data_drosophila_DEseq2 <- readRDS(file = "/datastore/2505621h/DE_Analysis/100cell_Analysis/comp_data_drosophila_pseudo_DESeq2.rds")
comp_data_drosophila_DEseq2@method.names$full.name <- "DESeq2"
saveRDS(comp_data_drosophila_DEseq2, file = "/datastore/2505621h/DE_Analysis/100cell_Analysis/comp_data_drosophila_pseudo_DESeq2.rds")

# DEsingle
method_name_DEsingle <- list(short.name = "DEsingle",
                             full.name = "DEsingle") 
comp_data_drosophila_DEsingle <- compcodeR::compData(drosophila, sample.annotations = coldata_compcodeR, info.parameters = parameters, variable.annotations = DE_State_annotation, method.names = method_name_DEsingle, result.table = dros_result_groups)

saveRDS(comp_data_drosophila_DEsingle, file = "/datastore/2505621h/100cell_Analysis/comp_data_drosophila_DEsingle.rds")

# MAST
method_name_MAST <- list(short.name = "MAST",
                         full.name = "MAST") 

comp_data_drosophila_MAST <- compcodeR::compData(drosophila, sample.annotations = coldata_compcodeR, info.parameters = parameters, variable.annotations = DE_State_annotation, method.names = method_name_MAST, result.table = fcHurdle)

saveRDS(comp_data_drosophila_MAST, file = "/datastore/2505621h/100cell_Analysis/comp_data_drosophila_MAST.rds")

# RankStat

# RP
method_name_RP <- list(short.name = "RP",
                       full.name = "Rank Products") 
comp_data_drosophila_RP <- compcodeR::compData(drosophila, sample.annotations = coldata_compcodeR, info.parameters = parameters, variable.annotations = DE_State_annotation, method.names = method_name_RP, result.table = result_RP)

saveRDS(comp_data_drosophila_RP, file = "/datastore/2505621h/100cell_Analysis/comp_data_drosophila_RP.rds")

# RSum
method_name_RSum <- list(short.name = "RS",
                         full.name = "Rank Sums") 
comp_data_drosophila_RSum <- compcodeR::compData(drosophila, sample.annotations = coldata_compcodeR, info.parameters = parameters, variable.annotations = DE_State_annotation, method.names = method_name_RSum, result.table = result_RSum)

saveRDS(comp_data_drosophila_RSum, file = "/datastore/2505621h/100cell_Analysis/comp_data_drosophila_RSum.rds")

# DRP
method_name_DRP <- list(short.name = "DRP",
                        full.name = "Difference of Rank Products")

comp_data_drosophila_DRP <- compcodeR::compData(drosophila, sample.annotations = coldata_compcodeR, info.parameters = parameters, variable.annotations = DE_State_annotation, method.names = method_name_DRP, result.table = result_DRP)

saveRDS(comp_data_drosophila_DRP, file = "/datastore/2505621h/100cell_Analysis/comp_data_drosophila_DRP.rds")

# RD
method_name_RD <- list(short.name = "RD",
                       full.name = "Rank Distances")

comp_data_drosophila_RD <- compcodeR::compData(drosophila, sample.annotations = coldata_compcodeR, info.parameters = parameters, variable.annotations = DE_State_annotation, method.names = method_name_RD, result.table = result_RD)

saveRDS(comp_data_drosophila_RD, file = "/datastore/2505621h/100cell_Analysis/comp_data_drosophila_RD.rds")

# R_RD
method_name_R_RD <- list(short.name = "R_RD",
                         full.name = "Reverse Rank Distances")

comp_data_drosophila_R_RD <- compcodeR::compData(drosophila, sample.annotations = coldata_compcodeR, info.parameters = parameters, variable.annotations = DE_State_annotation, method.names = method_name_R_RD, result.table = result_R_RD)

saveRDS(comp_data_drosophila_R_RD, file = "/datastore/2505621h/100cell_Analysis/comp_data_drosophila_R_RD.rds")

#SigEMD
method_name_SigEMD <- list(short.name = "SigEMD",
                           full.name = "SigEMD")

comp_data_drosophila_SigEMD <- compcodeR::compData(drosophila, sample.annotations = coldata_compcodeR, info.parameters = parameters, variable.annotations = DE_State_annotation, method.names = method_name_SigEMD, result.table = emd)

saveRDS(comp_data_drosophila_SigEMD, file = "/datastore/2505621h/100cell_Analysis/comp_data_drosophila_SigEMD.rds")

## Comparison

file.table <- data.frame(input.files = c("/datastore/2505621h/DE_Analysis/100cell_Analysis/comp_data_drosophila_pseudo_DESeq2.rds",
                                         "/datastore/2505621h/DE_Analysis/100cell_Analysis/comp_data_drosophila_SigEMD.rds",
                                         "/datastore/2505621h/DE_Analysis/100cell_Analysis/comp_data_drosophila_MAST.rds",
                                         "/datastore/2505621h/DE_Analysis/100cell_Analysis/comp_data_drosophila_DEsingle.rds",
                                         "/datastore/2505621h/DE_Analysis/100cell_Analysis/comp_data_drosophila_RP.rds",
                                         "/datastore/2505621h/DE_Analysis/100cell_Analysis/comp_data_drosophila_RSum.rds",
                                         "/datastore/2505621h/DE_Analysis/100cell_Analysis/comp_data_drosophila_DRP.rds",
                                         "/datastore/2505621h/DE_Analysis/100cell_Analysis/comp_data_drosophila_RD.rds",
                                         "/datastore/2505621h/DE_Analysis/100cell_Analysis/comp_data_drosophila_R_RD.rds"),
                         stringsAsFactors = FALSE)

parameters_comp <- list(incl.nbr.samples = NULL, incl.replicates = NULL,
                        incl.dataset = "Drosophila", incl.de.methods = NULL,
                        fdr.threshold = 0.05, tpr.threshold = 0.05,
                        typeI.threshold = 0.05, ma.threshold = 0.05,
                        fdc.maxvar = 1500, overlap.threshold = 0.05,
                        fracsign.threshold = 0.05,
                        comparisons = c("nbrtpfp", "overlap", "auc", "fdr", "tpr", "rocall", "correlation", "sorensen"))
compcodeR::runComparison(file.table = file.table, parameters = parameters_comp, output.directory = "/datastore/2505621h/DE_Analysis/100cell_Analysis/")

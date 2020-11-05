# Calculating tp, tn. fp, fn

# MAST
ground_truth <- readRDS(file = "/datastore/2505621h/100cell_Analysis/ground_truth_list_FINAL")
mast_results <- readRDS(file = "/datastore/2505621h/100cell_Analysis/drosophila_MAST_results")

names(mast_results)
DE <- ground_truth[[1]]
nonDE <- ground_truth[[2]]

mast_results <- tibble::rownames_to_column(mast_results, var = "gene")
sig_mast_results <- subset(mast_results, mast_results$FDR < 0.05)
non_sig_mast_results <- subset(mast_results, mast_results$FDR >= 0.05)

tp_mast <- sig_mast_results[sig_mast_results$gene %in% DE$`unlist(DE_genes)`,]
fp_mast <- sig_mast_results[sig_mast_results$gene %in% nonDE$`unlist(nonDE_genes)`,]
tn_mast <- non_sig_mast_results[non_sig_mast_results$gene %in% nonDE$`unlist(nonDE_genes)`,]
fn_mast <- non_sig_mast_results[non_sig_mast_results$gene %in% DE$`unlist(DE_genes)`,]

# SigEMD
emd_results <- readRDS(file = "/datastore/2505621h/100cell_Analysis/drosophila_SigEMD_results")

emd_results <- tibble::rownames_to_column(emd_results, var = "gene")
sig_emd_results <- subset(emd_results, emd_results$adjpvalue < 0.05)
non_sig_emd_results <- subset(emd_results, emd_results$adjpvalue>= 0.05)

tp_emd <- sig_emd_results[sig_emd_results$gene %in% DE$`unlist(DE_genes)`,]
fp_emd <- sig_emd_results[sig_emd_results$gene %in% nonDE$`unlist(nonDE_genes)`,]
tn_emd <- non_sig_emd_results[non_sig_emd_results$gene %in% nonDE$`unlist(nonDE_genes)`,]
fn_emd <- non_sig_emd_results[non_sig_emd_results$gene %in% DE$`unlist(DE_genes)`,]

# DEsingle
deSingle_results <- readRDS(file = "/datastore/2505621h/100cell_Analysis/deSingle_results_drosophila_FINAL")

deSingle_results <- tibble::rownames_to_column(deSingle_results, var = "gene")
sig_deSingle_results <- subset(deSingle_results, deSingle_results$adjpvalue < 0.05)
non_sig_deSingle_results <- subset(deSingle_results, deSingle_results$adjpvalue>= 0.05)

tp_deSingle <- sig_deSingle_results[sig_deSingle_results$gene %in% DE$`unlist(DE_genes)`,]
fp_deSingle <- sig_deSingle_results[sig_deSingle_results$gene %in% nonDE$`unlist(nonDE_genes)`,]
tn_deSingle <- non_sig_deSingle_results[non_sig_deSingle_results$gene %in% nonDE$`unlist(nonDE_genes)`,]
fn_deSingle <- non_sig_deSingle_results[non_sig_deSingle_results$gene %in% DE$`unlist(DE_genes)`,]

# RankStat

# RankProducts
RP_results <- readRDS(file = "/datastore/2505621h/100cell_Analysis/drosophila_result_RP_results")

RP_results <- tibble::rownames_to_column(RP_results, var = "gene")
sig_RP_results <- subset(RP_results, RP_results$FDR < 0.05)
non_sig_RP_results <- subset(RP_results, RP_results$FDR>= 0.05)

tp_RP <- sig_RP_results[sig_RP_results$gene %in% DE$`unlist(DE_genes)`,]
fp_RP <- sig_RP_results[sig_RP_results$gene %in% nonDE$`unlist(nonDE_genes)`,]
tn_RP <- non_sig_RP_results[non_sig_RP_results$gene %in% nonDE$`unlist(nonDE_genes)`,]
fn_RP <- non_sig_RP_results[non_sig_RP_results$gene %in% DE$`unlist(DE_genes)`,]
# RankSum
RS_results <- readRDS(file = "/datastore/2505621h/100cell_Analysis/drosophila_result_RSum_results")

RS_results <- tibble::rownames_to_column(RS_results, var = "gene")
sig_RS_results <- subset(RS_results, RS_results$FDR < 0.05)
non_sig_RS_results <- subset(RS_results, RS_results$FDR>= 0.05)

tp_RS <- sig_RS_results[sig_RS_results$gene %in% DE$`unlist(DE_genes)`,]
fp_RS <- sig_RS_results[sig_RS_results$gene %in% nonDE$`unlist(nonDE_genes)`,]
tn_RS <- non_sig_RS_results[non_sig_RS_results$gene %in% nonDE$`unlist(nonDE_genes)`,]
fn_RS <- non_sig_RS_results[non_sig_RS_results$gene %in% DE$`unlist(DE_genes)`,]
# RankDistance
RD_results <- readRDS(file = "/datastore/2505621h/100cell_Analysis/drosophila_result_RD_results")

RD_results <- tibble::rownames_to_column(RD_results, var = "gene")
sig_RD_results <- subset(RD_results, RD_results$FDR < 0.05)
non_sig_RD_results <- subset(RD_results, RD_results$FDR>= 0.05)

tp_RD <- sig_RD_results[sig_RD_results$gene %in% DE$`unlist(DE_genes)`,]
fp_RD <- sig_RD_results[sig_RD_results$gene %in% nonDE$`unlist(nonDE_genes)`,]
tn_RD <- non_sig_RD_results[non_sig_RD_results$gene %in% nonDE$`unlist(nonDE_genes)`,]
fn_RD <- non_sig_RD_results[non_sig_RD_results$gene %in% DE$`unlist(DE_genes)`,]
# ReverseRankDistance
R_RD_results <- readRDS(file = "/datastore/2505621h/100cell_Analysis/drosophila_result_R_RD_results")

R_RD_results <- tibble::rownames_to_column(R_RD_results, var = "gene")
sig_R_RD_results <- subset(R_RD_results, R_RD_results$FDR < 0.05)
non_sig_R_RD_results <- subset(R_RD_results, R_RD_results$FDR>= 0.05)

tp_R_RD <- sig_R_RD_results[sig_R_RD_results$gene %in% DE$`unlist(DE_genes)`,]
fp_R_RD <- sig_R_RD_results[sig_R_RD_results$gene %in% nonDE$`unlist(nonDE_genes)`,]
tn_R_RD <- non_sig_R_RD_results[non_sig_R_RD_results$gene %in% nonDE$`unlist(nonDE_genes)`,]
fn_R_RD <- non_sig_R_RD_results[non_sig_R_RD_results$gene %in% DE$`unlist(DE_genes)`,]
# DifferentialRankProducts
DRP_results <- readRDS(file = "/datastore/2505621h/100cell_Analysis/drosophila_result_DRP_results")

DRP_results <- tibble::rownames_to_column(DRP_results, var = "gene")
sig_DRP_results <- subset(DRP_results, DRP_results$FDR < 0.05)
non_sig_DRP_results <- subset(DRP_results, DRP_results$FDR>= 0.05)

tp_DRP <- sig_DRP_results[sig_DRP_results$gene %in% DE$`unlist(DE_genes)`,]
fp_DRP <- sig_DRP_results[sig_DRP_results$gene %in% nonDE$`unlist(nonDE_genes)`,]
tn_DRP <- non_sig_DRP_results[non_sig_DRP_results$gene %in% nonDE$`unlist(nonDE_genes)`,]
fn_DRP <- non_sig_DRP_results[non_sig_DRP_results$gene %in% DE$`unlist(DE_genes)`,]

# DESeq2
comp_DEseq2 <- readRDS(file = "/datastore/2505621h/100cell_Analysis/comp_data_drosophila_pseudo_DESeq2.rds" )
deseq2_results <- comp_DEseq2@result.table
deseq2_results[is.na(deseq2_results)] <- 1

deseq2_results <- tibble::rownames_to_column(deseq2_results, var = "gene")
sig_deseq2_results <- subset(deseq2_results, deseq2_results$adjpvalue < 0.05)
non_sig_deseq2_results <- subset(DRP_results, deseq2_results$pvalue>= 0.05)

tp_deseq2 <- sig_deseq2_results[sig_deseq2_results$gene %in% DE$`unlist(DE_genes)`,]
fp_deseq2 <- sig_deseq2_results[sig_deseq2_results$gene %in% nonDE$`unlist(nonDE_genes)`,]
tn_deseq2 <- non_sig_deseq2_results[non_sig_deseq2_results$gene %in% nonDE$`unlist(nonDE_genes)`,]
fn_deseq2 <- non_sig_deseq2_results[non_sig_deseq2_results$gene %in% DE$`unlist(DE_genes)`,]

# Read data
brain10X <- Read10X(data.dir = "/datastore/2505621h/Simulators/Data/brain10X_filtered_feature_bc_matrix/")

brain10X <- as.matrix(brain10X)

brain10X_sce <- SingleCellExperiment(assays = list(counts = brain10X))

# SPARSim
brain10X_norm <- LogNormalize(brain10X)
brain10X_norm_sce <- SingleCellExperiment(assays = list(counts = brain10X_norm))
brain10X_norm <- as.matrix(brain10X_norm)
brain_condition_A_column_index <- c(1:11843)
brain10X_conditions <- list(cond_A = brain_condition_A_column_index)

brain10X_parameters <- SPARSim_estimate_parameter_from_data(raw_data = brain10X,
                                                            norm_data = brain10X_norm,
                                                            conditions = brain10X_conditions)

brain10X_result <- SPARSim_simulation(dataset_parameter = brain10X_parameters)
brain_sparsim <- brain10X_result$count_matrix
sparsim_brain_sce <- SingleCellExperiment(assays = list(counts = brain_sparsim))

# Splat
splat_parameters_brain <- splatEstimate(brain10X)
splat_sim_brain <- splatSimulate(splat_parameters_brain)

# scDesign
scDesign_brain <- design_data(brain10X, S = 1e+08, ncell = 11843, ngroup = 1)
scDesign_brain_sce <- SingleCellExperiment(assays = list(counts = scDesign_brain))

# SPsimSeq
spsim_brain_sim <- SPsimSeq(n.sim = 1, 
                            s.data = brain10X_sce,
                            n.genes = 31053,
                            batch.config = 1,
                            group.config = 1,
                            tot.samples = 11843,
                            pDE = 0.1,
                            lfc.thrld = 0.25,
                            t.thrld = 1.5,
                            model.zero.prob = TRUE,
                            result.format = "SCE")

spsim_brain_sce <- spsim_brain_sim[[1]]

# Comparison & Difference
comparison_brain <- compareSCEs(list(Real = brain10X_sce,
                                     Splat = splat_sim_brain, 
                                     SPARSim = sparsim_brain_sce,
                                     scDesign = scDesign_brain_sce,
                                     SPsimSeq = spsim_brain_sce))

diff_comp_brain <- diffSCEslog10(list(Real = real_brain_sce,
                                      Splat = splat_sim_brain, 
                                      SPARSim = sparsim_brain_sce,
                                      scDesign = scDesign_brain_sce,
                                      SPsimSeq = spsim_brain_sce), ref = "Real")

summary_brain <- summariseDiff(diff_comp_brain)

saveRDS(summary_brain, file = "/datastore/2505621h/Simulators/RDS/summ_brain")
saveRDS(brain10X, file = "/datastore/2505621h/Simulators/RDS/real_count_matrix_brain")

saveRDS(brain10X_sce, file = "/datastore/2505621h/Simulators/RDS/real_brain_object")
real_brain_sce <- readRDS(file = "/datastore/2505621h/Simulators/RDS/real_brain_object")

saveRDS(sparsim_brain_sce, file = "/datastore/2505621h/Simulators/RDS/SPARSim_brain_object")
sparsim_brain_sce <- readRDS(file = "/datastore/2505621h/Simulators/RDS/SPARSim_brain_object")

saveRDS(splat_sim_brain, file = "/datastore/2505621h/Simulators/RDS/splat_brain_object")
splat_sim_brain <- readRDS(file = "/datastore/2505621h/Simulators/RDS/splat_brain_object")

saveRDS(scDesign_brain_sce, file = "/datastore/2505621h/Simulators/RDS/scDesign_brain_object")
scDesign_brain_sce <- readRDS(file = "/datastore/2505621h/Simulators/RDS/scDesign_brain_object" )

saveRDS(spsim_brain_sce, file = "/datastore/2505621h/Simulators/RDS/SPsimSeq_brain_object")
spsim_brain_sce <- readRDS(file = "/datastore/2505621h/Simulators/RDS/SPsimSeq_brain_object")

saveRDS(comparison_brain, file = "/datastore/2505621h/Simulators/RDS/comparison_brain")
saveRDS(diff_comp_brain, file = "/datastore/2505621h/Simulators/RDS/diff_brain")
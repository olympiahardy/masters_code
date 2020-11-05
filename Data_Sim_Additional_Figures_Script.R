# PBMC
comparison_pbmc <- readRDS(file = "/datastore/2505621h/Simulators/RDS/comparison_pbmc_log10")
diff_comp_pbmc <- readRDS(file = "/datastore/2505621h/Simulators/RDS/diff_pbmc")

# Non-zero plots

mean_plot_pbmc <- comparison_pbmc$Plots$Means
mean_plot_pbmc <- mean_plot_pbmc + scale_y_continuous(limits = c(0, 3)) + ggtitle(label = 'PBMC 10X: Distribution of Mean Expression')

diff_mean_plot_pbmc <- diff_comp_pbmc$Plots$Means
diff_mean_plot_pbmc <- diff_mean_plot_pbmc + ggtitle(label = 'PBMC 10X: Difference in Mean Expression')


variance_plot_pbmc <- comparison_pbmc$Plots$Variances
variance_plot_pbmc <- variance_plot_pbmc + scale_y_continuous(limits = c(0, 20)) + ggtitle(label = 'PBMC 10X: Distribution of Variance')


diff_variance_plot_pbmc <- diff_comp_pbmc$Plots$Variances
diff_variance_plot_pbmc <- diff_variance_plot_pbmc + ggtitle(label = 'PBMC 10X: Difference in Variance')


library_size_pbmc <- comparison_pbmc$Plots$LibrarySizes
library_size_pbmc <- library_size_pbmc + scale_y_continuous(limits = c(0, 30000))+ ggtitle(label = 'PBMC 10X: Distribution of Library Sizes')


diff_library_size_pbmc <- diff_comp_pbmc$Plots$LibrarySizes
diff_library_size_pbmc <- diff_library_size_pbmc + ggtitle(label = 'PBMC 10X: Difference in Library Sizes')


mean_variance_plot_pbmc <- comparison_pbmc$Plots$MeanVar
mean_variance_plot_pbmc <- mean_variance_plot_pbmc + ggtitle(label = 'PBMC 10X: Mean-Variance Relationship')


diff_mean_variance_plot_pbmc <- diff_comp_pbmc$Plots$MeanVar
diff_mean_variance_plot_pbmc <- diff_mean_variance_plot_pbmc + ggtitle(label = 'PBMC 10X: Difference in Mean-Variance Relationship')


plot1_pbmc <- cowplot::plot_grid(mean_plot_pbmc,
                           diff_mean_plot_pbmc,
                           variance_plot_pbmc,
                           diff_variance_plot_pbmc,
                           nrow = 2,
                           ncol = 2,
                           align = "hv")
plot1_pbmc

plot2_pbmc <- cowplot::plot_grid(library_size_pbmc,
                             diff_library_size_pbmc,
                             mean_variance_plot_pbmc,
                             diff_mean_variance_plot_pbmc,
                             nrow = 2,
                             ncol = 2,
                             align = "hv")
plot2_pbmc

# Zero plots 

zero_gene_plot_pbmc <- comparison_pbmc$Plots$ZerosGene
zero_gene_plot_pbmc <- zero_gene_plot_pbmc + scale_y_continuous(name = 'Sparsity [%]', limits = c(75, 100)) + ggtitle(label = 'PBMC 10X: Sparsity by gene')

diff_zero_gene_plot_pbmc <- diff_comp_pbmc$Plots$ZerosGene
diff_zero_gene_plot_pbmc <- diff_zero_gene_plot_pbmc + ggtitle(label = 'PBMC 10X: Difference in sparsity by gene')

zero_cell_plot_pbmc <- comparison_pbmc$Plots$ZerosCell
zero_cell_plot_pbmc <- zero_cell_plot_pbmc + scale_y_continuous(name = 'Sparsity [%]', limits = c(70, 100)) + ggtitle(label = 'PBMC 10X: Sparsity by cell')

diff_zero_cell_plot_pbmc <- diff_comp_pbmc$Plots$ZerosCell
diff_zero_cell_plot_pbmc <- diff_zero_cell_plot_pbmc + ggtitle(label = 'PBMC 10X: Difference in sparsity by cell')

mean_zeros_plot_pbmc <- comparison_pbmc$Plots$MeanZeros
mean_zeros_plot_pbmc <- mean_zeros_plot_pbmc + ggtitle(label = 'PBMC 10X: Mean-Sparsity Relationship')

diff_mean_zeros_plot_pbmc <- diff_comp_pbmc$Plots$MeanZeros
diff_mean_zeros_plot_pbmc <- diff_mean_zeros_plot_pbmc + ggtitle(label = 'PBMC 10X: Difference in Mean-Sparsity Relationship')


plot3_pbmc <- cowplot::plot_grid(zero_gene_plot_pbmc,
                                 diff_zero_gene_plot_pbmc,
                                 zero_cell_plot_pbmc,
                                 diff_zero_cell_plot_pbmc,
                                 mean_zeros_plot_pbmc,
                                 diff_mean_zeros_plot_pbmc,
                                 nrow = 3,
                                 ncol = 2,
                                 align = "hv")

plot3_pbmc

# Brain
diff_comp_brain <- readRDS(file = "/datastore/2505621h/Simulators/RDS/diff_brain")
comparison_brain <- readRDS(file = "/datastore/2505621h/Simulators/RDS/comparison_brain_log10")

# Non-zero Plots

mean_plot_brain <- comparison_brain$Plots$Means
mean_plot_brain <- mean_plot_brain + scale_y_continuous(limits = c(0, 3))+ ggtitle(label = 'Brain 10X: Distribution of Mean Expression')

diff_mean_plot_brain <- diff_comp_brain$Plots$Means
diff_mean_plot_brain <- diff_mean_plot_brain + ggtitle(label = 'Brain 10X: Difference in Mean Expression')

variance_plot_brain <- comparison_brain$Plots$Variances
variance_plot_brain <- variance_plot_brain + scale_y_continuous(limits = c(0, 20))+ ggtitle(label = 'Brain 10X: Distribution of Variance')

diff_variance_plot_brain <- diff_comp_brain$Plots$Variances
diff_variance_plot_brain <- diff_variance_plot_brain + ggtitle(label = 'Brain 10X: Difference in Variance')

library_size_brain <- comparison_brain$Plots$LibrarySizes
library_size_brain <- library_size_brain + scale_y_continuous(limits = c(0, 100000))+ ggtitle(label = 'Brain 10X: Distribution of Library Sizes')

diff_library_size_brain <- diff_comp_brain$Plots$LibrarySizes
diff_library_size_brain <- diff_library_size_brain + ggtitle(label = 'Brain 10X: Difference in Library Sizes')

mean_variance_plot_brain <- comparison_brain$Plots$MeanVar
mean_variance_plot_brain <- mean_variance_plot_brain + ggtitle(label = 'Brain 10X: Mean-Variance Relationship')

diff_mean_variance_plot_brain <- diff_comp_brain$Plots$MeanVar
diff_mean_variance_plot_brain <- diff_mean_variance_plot_brain + ggtitle(label = 'Brain 10X: Difference in Mean-Variance Relationship')

plot1_brain <- cowplot::plot_grid(mean_plot_brain,
                                 diff_mean_plot_brain,
                                 variance_plot_brain,
                                 diff_variance_plot_brain,
                                 nrow = 2,
                                 ncol = 2,
                                 align = "hv")
plot1_brain

plot2_brain <- cowplot::plot_grid(library_size_brain,
                                 diff_library_size_brain,
                                 mean_variance_plot_brain,
                                 diff_mean_variance_plot_brain,
                                 nrow = 2,
                                 ncol = 2,
                                 align = "hv")
plot2_brain

# Zero Plots

zero_gene_plot_brain <- comparison_brain$Plots$ZerosGene
zero_gene_plot_brain <- zero_gene_plot_brain + scale_y_continuous(name = 'Sparsity [%]', limits = c(75, 100)) + ggtitle(label = 'Brain 10X: Sparsity by gene')

diff_zero_gene_plot_brain <- diff_comp_brain$Plots$ZerosGene
diff_zero_gene_plot_brain <- diff_zero_gene_plot_brain + ggtitle(label = 'Brain 10X: Difference in sparsity by gene')

zero_cell_plot_brain <- comparison_brain$Plots$ZerosCell
zero_cell_plot_brain <- zero_cell_plot_brain + scale_y_continuous(name = 'Sparsity [%]', limits = c(70, 100)) + ggtitle(label = 'Brain 10X: Sparsity by cell')

diff_zero_cell_plot_brain <- diff_comp_brain$Plots$ZerosCell
diff_zero_cell_plot_brain <- zero_cell_plot_brain + ggtitle(label = 'Brain 10X: Difference in sparsity by cell')

mean_zeros_plot_brain <- comparison_brain$Plots$MeanZeros
mean_zeros_plot_brain <- mean_zeros_plot_brain + ggtitle(label = 'Brain 10X: Mean-Sparsity Relationship')

diff_mean_zeros_plot_brain <- diff_comp_brain$Plots$MeanZeros
diff_mean_zeros_plot_brain <- diff_mean_zeros_plot_brain + ggtitle(label = 'Brain 10X: Difference in Mean-Sparsity Relationship')

plot3_brain <- cowplot::plot_grid(zero_gene_plot_brain,
                                 diff_zero_gene_plot_brain,
                                 zero_cell_plot_brain,
                                 diff_zero_cell_plot_brain,
                                 mean_zeros_plot_brain,
                                 diff_mean_zeros_plot_brain,
                                 nrow = 3,
                                 ncol = 2,
                                 align = "hv")

plot3_brain

# T-cell
comparison_tcell <- readRDS(file = "/datastore/2505621h/Simulators/RDS/comparison_tcell_log10")
diff_comp_tcell <- readRDS(file = "/datastore/2505621h/Simulators/RDS/diff_tcell")

# Non-zero plots

mean_plot_tcell <- comparison_tcell$Plots$Means
mean_plot_tcell <- mean_plot_tcell + scale_y_continuous(limits = c(0, 3)) + ggtitle(label = 'T-cell 10X: Distribution of Mean Expression')

diff_mean_plot_tcell <- diff_comp_tcell$Plots$Means
diff_mean_plot_tcell <- diff_mean_plot_tcell + ggtitle(label = 'T-cell 10X: Difference in Mean Expression')

variance_plot_tcell <- comparison_tcell$Plots$Variances
variance_plot_tcell <- variance_plot_tcell + scale_y_continuous(limits = c(0, 20))+ ggtitle(label = 'T-cell 10X: Distribution of Variance')

diff_variance_plot_tcell <- diff_comp_tcell$Plots$Variances
diff_variance_plot_tcell <- diff_variance_plot_tcell + ggtitle(label = 'T-cell 10X: Difference in Variance')

library_size_tcell <- comparison_tcell$Plots$LibrarySizes
library_size_tcell <- library_size_tcell + scale_y_continuous(limits = c(0, 30000))+ ggtitle(label = 'T-cell 10X: Distribution of Library Sizes')

diff_library_size_tcell <- diff_comp_tcell$Plots$LibrarySizes
diff_library_size_tcell <- diff_library_size_tcell + ggtitle(label = 'T-cell 10X: Difference in Library Sizes')

mean_variance_plot_tcell <- comparison_tcell$Plots$MeanVar
mean_variance_plot_tcell <- mean_variance_plot_tcell + ggtitle(label = 'T-cell 10X: Mean-Variance Relationship')

diff_mean_variance_plot_tcell <- diff_comp_tcell$Plots$MeanVar
diff_mean_variance_plot_tcell <- diff_mean_variance_plot_tcell + ggtitle(label = 'T-cell 10X: Difference in Mean-Variance Relationship')

plot1_tcell <- cowplot::plot_grid(mean_plot_tcell,
                                  diff_mean_plot_tcell,
                                  variance_plot_tcell,
                                  diff_variance_plot_tcell,
                                  nrow = 2,
                                  ncol = 2,
                                  align = "hv")
plot1_tcell

plot2_tcell <- cowplot::plot_grid(library_size_tcell,
                                  diff_library_size_tcell,
                                  mean_variance_plot_tcell,
                                  diff_mean_variance_plot_tcell,
                                  nrow = 2,
                                  ncol = 2,
                                  align = "hv")
plot2_tcell

# Zero plots

zero_gene_plot_tcell <- comparison_tcell$Plots$ZerosGene
zero_gene_plot_tcell <- zero_gene_plot_tcell + scale_y_continuous(name = 'Sparsity [%]', limits = c(75, 100)) + ggtitle(label = 'T-cell 10X: Sparsity by gene')

diff_zero_gene_plot_tcell <- diff_comp_tcell$Plots$ZerosGene
diff_zero_gene_plot_tcell <- diff_zero_gene_plot_tcell + ggtitle(label = 'T-cell 10X: Difference in sparsity by gene')

zero_cell_plot_tcell <- comparison_tcell$Plots$ZerosCell
zero_cell_plot_tcell <- zero_cell_plot_tcell + scale_y_continuous(name = 'Sparsity [%]', limits = c(70, 100)) + ggtitle(label = 'T-cell 10X: Sparsity by cell')

diff_zero_cell_plot_tcell <- diff_comp_tcell$Plots$ZerosCell
diff_zero_cell_plot_tcell <- diff_zero_cell_plot_tcell + ggtitle(label = 'T-cell 10X: Difference in sparsity by cell')

mean_zeros_plot_tcell <- comparison_tcell$Plots$MeanZeros
mean_zeros_plot_tcell <- mean_zeros_plot_tcell + ggtitle(label = 'T-cell 10X: Mean-Sparsity Relationship')

diff_mean_zeros_plot_tcell <- diff_comp_tcell$Plots$MeanZeros
diff_mean_zeros_plot_tcell <- diff_mean_zeros_plot_tcell + ggtitle(label = 'T-cell 10X: Difference in Mean-Sparsity Relationship')

plot3_tcell <- cowplot::plot_grid(zero_gene_plot_tcell,
                                  diff_zero_gene_plot_tcell,
                                  zero_cell_plot_tcell,
                                  diff_zero_cell_plot_tcell,
                                  mean_zeros_plot_tcell,
                                  diff_mean_zeros_plot_tcell,
                                  nrow = 3,
                                  ncol = 2,
                                  align = "hv")

plot3_tcell

# Sheep
comparison_sheep <- readRDS(file = "/datastore/2505621h/Simulators/RDS/comparison_sheep")
diff_comp_sheep <- readRDS(file = "/datastore/2505621h/Simulators/RDS/diff_sheep")

# Non-zero plots

mean_plot_sheep <- comparison_sheep$Plots$Means
mean_plot_sheep <- mean_plot_sheep + scale_y_continuous(limits = c(0, 3)) + ggtitle(label = 'Sheep 10X: Distribution of Mean Expression')

diff_mean_plot_sheep <- diff_comp_sheep$Plots$Means
diff_mean_plot_sheep <- diff_mean_plot_sheep + ggtitle(label = 'Sheep 10X: Difference in Mean Expression')

variance_plot_sheep <- comparison_sheep$Plots$Variances
variance_plot_sheep <- variance_plot_sheep + scale_y_continuous(limits = c(0, 20)) + ggtitle(label = 'Sheep 10X: Distribution of Variance')

diff_variance_plot_sheep <- diff_comp_sheep$Plots$Variances
diff_variance_plot_sheep <- diff_variance_plot_sheep + ggtitle(label = 'Sheep 10X: Difference in Variance')

library_size_sheep <- comparison_sheep$Plots$LibrarySizes
library_size_sheep <- library_size_sheep + scale_y_continuous(limits = c(0, 50000))+ ggtitle(label = 'Sheep 10X: Distribution of Library Sizes')

diff_library_size_sheep <- diff_comp_sheep$Plots$LibrarySizes
diff_library_size_sheep <- diff_library_size_sheep + ggtitle(label = 'Sheep 10X: Difference in Library Sizes')

mean_variance_plot_sheep <- comparison_sheep$Plots$MeanVar
mean_variance_plot_sheep <- mean_variance_plot_sheep + ggtitle(label = 'Sheep 10X: Mean-Variance Relationship')

diff_mean_variance_plot_sheep <- diff_comp_sheep$Plots$MeanVar
diff_mean_variance_plot_sheep <- diff_mean_variance_plot_sheep + ggtitle(label = 'Sheep 10X: Difference in Mean-Variance Relationship')

plot1_sheep <- cowplot::plot_grid(mean_plot_sheep,
                                  diff_mean_plot_sheep,
                                  variance_plot_sheep,
                                  diff_variance_plot_sheep,
                                  nrow = 2,
                                  ncol = 2,
                                  align = "hv")
plot1_sheep

plot2_sheep <- cowplot::plot_grid(library_size_sheep,
                                  diff_library_size_sheep,
                                  mean_variance_plot_sheep,
                                  diff_mean_variance_plot_sheep,
                                  nrow = 2,
                                  ncol = 2,
                                  align = "hv")
plot2_sheep

# Zero plots

zero_gene_plot_sheep <- comparison_sheep$Plots$ZerosGene
zero_gene_plot_sheep <- zero_gene_plot_sheep + scale_y_continuous(name = 'Sparsity [%]', limits = c(75, 100)) + ggtitle(label = 'Sheep 10X: Sparsity by gene')

diff_zero_gene_plot_sheep <- diff_comp_sheep$Plots$ZerosGene
diff_zero_gene_plot_sheep <- diff_zero_gene_plot_sheep + ggtitle(label = 'Sheep 10X: Difference in sparsity by gene')

zero_cell_plot_sheep <- comparison_sheep$Plots$ZerosCell
zero_cell_plot_sheep <- zero_cell_plot_sheep + scale_y_continuous(name = 'Sparsity [%]', limits = c(70, 100)) + ggtitle(label = 'Sheep 10X: Sparsity by cell')

diff_zero_cell_plot_sheep <- diff_comp_sheep$Plots$ZerosCell
diff_zero_cell_plot_sheep <- diff_zero_cell_plot_sheep + ggtitle(label = 'Sheep 10X: Difference in sparsity by cell')

mean_zeros_plot_sheep <- comparison_sheep$Plots$MeanZeros
mean_zeros_plot_sheep <- mean_zeros_plot_sheep + ggtitle(label = 'Sheep 10X: Mean-Sparsity Relationship')

diff_mean_zeros_plot_sheep <- diff_comp_sheep$Plots$MeanZeros
diff_mean_zeros_plot_sheep <- diff_mean_zeros_plot_sheep + ggtitle(label = 'Sheep 10X: Difference in Mean-Sparsity Relationship')

plot3_sheep <- cowplot::plot_grid(zero_gene_plot_sheep,
                                  diff_zero_gene_plot_sheep,
                                  zero_cell_plot_sheep,
                                  diff_zero_cell_plot_sheep,
                                  mean_zeros_plot_sheep,
                                  diff_mean_zeros_plot_sheep,
                                  nrow = 3,
                                  ncol = 2,
                                  align = "hv")

plot3_sheep

# TGFB1
comparison_tgfb1 <- readRDS(file = "/datastore/2505621h/Simulators/RDS/comparison_tgfb1")
diff_comp_tgfb1 <- readRDS(file = "/datastore/2505621h/Simulators/RDS/diff_tgfb1")

# Non-zero plots

mean_plot_tgfb1 <- comparison_tgfb1$Plots$Means
mean_plot_tgfb1 <- mean_plot_tgfb1 + scale_y_continuous(limits = c(0, 3)) + ggtitle(label = 'TGFB1 10X: Distribution of Mean Expression')

diff_mean_plot_tgfb1 <- diff_comp_tgfb1$Plots$Means
diff_mean_plot_tgfb1 <- diff_mean_plot_tgfb1 + ggtitle(label = 'TGFB1 10X: Difference in Mean Expression')

variance_plot_tgfb1 <- comparison_tgfb1$Plots$Variances
variance_plot_tgfb1 <- variance_plot_tgfb1 + scale_y_continuous(limits = c(0, 20)) + ggtitle(label = 'TGFB1 10X: Distribution of Variance')

diff_variance_plot_tgfb1 <- diff_comp_tgfb1$Plots$Variances
diff_variance_plot_tgfb1 <- diff_variance_plot_tgfb1 + ggtitle(label = 'TGFB1 10X: Difference in Variance')

library_size_tgfb1 <- comparison_tgfb1$Plots$LibrarySizes
library_size_tgfb1 <- library_size_tgfb1 + scale_y_continuous(limits = c(0, 50000))+ ggtitle(label = 'TGFB1 10X: Distribution of Library Sizes')

diff_library_size_tgfb1 <- diff_comp_tgfb1$Plots$LibrarySizes
diff_library_size_tgfb1 <- diff_library_size_tgfb1 + ggtitle(label = 'TGFB1 10X: Difference in Library Sizes')

mean_variance_plot_tgfb1 <- comparison_tgfb1$Plots$MeanVar
mean_variance_plot_tgfb1 <- mean_variance_plot_tgfb1 + ggtitle(label = 'TGFB1 10X: Mean-Variance Relationship')

diff_mean_variance_plot_tgfb1 <- diff_comp_tgfb1$Plots$MeanVar
diff_mean_variance_plot_tgfb1 <- diff_mean_variance_plot_tgfb1 + ggtitle(label = 'TGFB1 10X: Difference in Mean-Variance Relationship')

plot1_tgfb1 <- cowplot::plot_grid(mean_plot_tgfb1,
                                  diff_mean_plot_tgfb1,
                                  variance_plot_tgfb1,
                                  diff_variance_plot_tgfb1,
                                  nrow = 2,
                                  ncol = 2,
                                  align = "hv")
plot1_tgfb1

plot2_tgfb1 <- cowplot::plot_grid(library_size_tgfb1,
                                  diff_library_size_tgfb1,
                                  mean_variance_plot_tgfb1,
                                  diff_mean_variance_plot_tgfb1,
                                  nrow = 2,
                                  ncol = 2,
                                  align = "hv")
plot2_tgfb1

# Zero plots

zero_gene_plot_tgfb1 <- comparison_tgfb1$Plots$ZerosGene
zero_gene_plot_tgfb1 <- zero_gene_plot_tgfb1 + scale_y_continuous(name = 'Sparsity [%]', limits = c(75, 100)) + ggtitle(label = 'TGFB1 10X: Sparsity by gene')

diff_zero_gene_plot_tgfb1 <- diff_comp_tgfb1$Plots$ZerosGene
diff_zero_gene_plot_tgfb1 <- diff_zero_gene_plot_tgfb1 + ggtitle(label = 'TGFB1 10X: Difference in sparsity by gene')

zero_cell_plot_tgfb1 <- comparison_tgfb1$Plots$ZerosCell
zero_cell_plot_tgfb1 <- zero_cell_plot_tgfb1 + scale_y_continuous(name = 'Sparsity [%]', limits = c(70, 100)) + ggtitle(label = 'TGFB1 10X: Sparsity by cell')

diff_zero_cell_plot_tgfb1 <- diff_comp_tgfb1$Plots$ZerosCell
diff_zero_cell_plot_tgfb1 <- diff_zero_cell_plot_tgfb1 + ggtitle(label = 'TGFB1 10X: Difference in sparsity by cell')

mean_zeros_plot_tgfb1 <- comparison_tgfb1$Plots$MeanZeros
mean_zeros_plot_tgfb1 <- mean_zeros_plot_tgfb1 + ggtitle(label = 'TGFB1 10X: Mean-Sparsity Relationship')

diff_mean_zeros_plot_tgfb1 <- diff_comp_tgfb1$Plots$MeanZeros
diff_mean_zeros_plot_tgfb1 <- diff_mean_zeros_plot_tgfb1 + ggtitle(label = 'TGFB1 10X: Difference in Mean-Sparsity Relationship')

plot3_tgfb1 <- cowplot::plot_grid(zero_gene_plot_tgfb1,
                                  diff_zero_gene_plot_tgfb1,
                                  zero_cell_plot_tgfb1,
                                  diff_zero_cell_plot_tgfb1,
                                  mean_zeros_plot_tgfb1,
                                  diff_mean_zeros_plot_tgfb1,
                                  nrow = 3,
                                  ncol = 2,
                                  align = "hv")

plot3_tgfb1

# Jurkat 
diff_comp_jurkat <- readRDS(file = "/datastore/2505621h/Simulators/RDS/diff_jurkat")
comparison_jurkat <- readRDS(file = "/datastore/2505621h/Simulators/RDS/comparison_jurkat_log10")


# Non-zero plots
mean_plot_jurkat <- comparison_jurkat$Plots$Means
mean_plot_jurkat <- mean_plot_jurkat + scale_y_continuous(limits = c(0, 3))+ ggtitle(label = 'Zheng 10X: Distribution of Mean Expression')

diff_mean_plot_jurkat <- diff_comp_jurkat$Plots$Means
diff_mean_plot_jurkat <- diff_mean_plot_jurkat + ggtitle(label = 'Zheng 10X: Difference in Mean Expression')

variance_plot_jurkat <- comparison_jurkat$Plots$Variances
variance_plot_jurkat <- variance_plot_jurkat + scale_y_continuous(limits = c(0, 20))+ ggtitle(label = 'Zheng 10X: Distribution of Variance')

diff_variance_plot_jurkat <- diff_comp_jurkat$Plots$Variances
diff_variance_plot_jurkat <- diff_variance_plot_jurkat + ggtitle(label = 'Zheng 10X: Difference in Variance')

library_size_jurkat <- comparison_jurkat$Plots$LibrarySizes
library_size_jurkat <- library_size_jurkat + scale_y_continuous(limits = c(0, 100000))+ ggtitle(label = 'Zheng 10X: Distribution of Library Sizes')

diff_library_size_jurkat <- diff_comp_jurkat$Plots$LibrarySizes
diff_library_size_jurkat <- diff_library_size_jurkat + ggtitle(label = 'Zheng 10X: Difference in Library Sizes')

mean_variance_plot_jurkat <- comparison_jurkat$Plots$MeanVar
mean_variance_plot_jurkat <- mean_variance_plot_jurkat + ggtitle(label = 'Zheng 10X: Mean-Variance Relationship')

diff_mean_variance_plot_jurkat <- diff_comp_jurkat$Plots$MeanVar
diff_mean_variance_plot_jurkat <- diff_mean_variance_plot_jurkat + ggtitle(label = 'Zheng 10X: Difference in Mean-Variance Relationship')

plot1_jurkat <- cowplot::plot_grid(mean_plot_jurkat,
                                  diff_mean_plot_jurkat,
                                  variance_plot_jurkat,
                                  diff_variance_plot_jurkat,
                                  nrow = 2,
                                  ncol = 2,
                                  align = "hv")
plot1_jurkat

plot2_jurkat <- cowplot::plot_grid(library_size_jurkat,
                                  diff_library_size_jurkat,
                                  mean_variance_plot_jurkat,
                                  diff_mean_variance_plot_jurkat,
                                  nrow = 2,
                                  ncol = 2,
                                  align = "hv")
plot2_jurkat

# Zero plots

zero_gene_plot_jurkat <- comparison_jurkat$Plots$ZerosGene
zero_gene_plot_jurkat <- zero_gene_plot_jurkat + scale_y_continuous(name = 'Sparsity [%]', limits = c(75, 100)) + ggtitle(label = 'Zheng 10X: Sparsity by gene')

diff_zero_gene_plot_jurkat <- diff_comp_jurkat$Plots$ZerosGene
diff_zero_gene_plot_jurkat <- diff_zero_gene_plot_jurkat + ggtitle(label = 'Zheng 10X: Difference in sparsity by gene')

zero_cell_plot_jurkat <- comparison_jurkat$Plots$ZerosCell
zero_cell_plot_jurkat <- zero_cell_plot_jurkat + scale_y_continuous(name = 'Sparsity [%]', limits = c(70, 100)) + ggtitle(label = 'Zheng 10X: Sparsity by cell')

diff_zero_cell_plot_jurkat <- diff_comp_jurkat$Plots$ZerosCell
diff_zero_cell_plot_jurkat <- diff_zero_cell_plot_jurkat + ggtitle(label = 'Zheng 10X: Difference in sparsity by cell')

mean_zeros_plot_jurkat <- comparison_jurkat$Plots$MeanZeros
mean_zeros_plot_jurkat <- mean_zeros_plot_jurkat + ggtitle(label = 'Zheng 10X: Mean-Sparsity Relationship')

diff_mean_zeros_plot_jurkat <- diff_comp_jurkat$Plots$MeanZeros
diff_mean_zeros_plot_jurkat <- diff_mean_zeros_plot_jurkat + ggtitle(label = 'Zheng 10X: Difference in Mean-Sparsity Relationship')

plot3_jurkat <- cowplot::plot_grid(zero_gene_plot_jurkat,
                                  diff_zero_gene_plot_jurkat,
                                  zero_cell_plot_jurkat,
                                  diff_zero_cell_plot_jurkat,
                                  mean_zeros_plot_jurkat,
                                  diff_mean_zeros_plot_jurkat,
                                  nrow = 3,
                                  ncol = 2,
                                  align = "hv")

plot3_jurkat


# Drosophila - included in main paper
comparison_drosoph <- readRDS(file = "/datastore/2505621h/Simulators/RDS/comparison_drosophila")
diff_drosoph <- readRDS(file = "/datastore/2505621h/Simulators/RDS/diff_drosophila")

# Non-zero plots
mean_plot_drosophila <- comparison_drosoph$Plots$Means
mean_plot_drosophila <- mean_plot_drosophila + scale_y_continuous(limits = c(0, 3)) + ggtitle(label = 'Drosophila 10X: Distribution of Mean Expression')

variance_plot_drosophila <- comparison_drosoph$Plots$Variances
variance_plot_drosophila <- variance_plot_drosophila + scale_y_continuous(limits = c(0, 20)) + ggtitle(label = 'Drosophila 10X: Distribution of Variance')

library_size_drosophila <- comparison_drosoph$Plots$LibrarySizes
library_size_drosophila <- library_size_drosophila + scale_y_continuous(limits = c(0, 50000))+ ggtitle(label = 'Drosophila 10X: Distribution of Library Sizes')

mean_variance_plot_drosophila <- comparison_drosoph$Plots$MeanVar
mean_variance_plot_drosophila <- mean_variance_plot_drosophila + ggtitle(label = 'Drosophila 10X: Mean-Variance Relationship')

diff_library_size_drosophila <- diff_drosoph$Plots$LibrarySizes
diff_library_size_drosophila <- diff_library_size_drosophila + ggtitle(label = 'Drosophila 10X: Difference in Library Sizes')

diff_mean_variance_plot_drosophila <- diff_drosoph$Plots$MeanVar
diff_mean_variance_plot_drosophila <- diff_mean_variance_plot_drosophila + ggtitle(label = 'Drosophila 10X: Difference in Mean-Variance Relationship')

diff_mean_plot_drosophila <- diff_drosoph$Plots$Means
diff_mean_plot_drosophila <- diff_mean_plot_drosophila + ggtitle(label = 'Drosophila 10X: Difference in Mean Expression')

diff_variance_plot_drosophila <- diff_drosoph$Plots$Variances
diff_variance_plot_drosophila <- diff_variance_plot_drosophila + ggtitle(label = 'Drosophila 10X: Difference in Variance')


plot1_drosoph <- cowplot::plot_grid(mean_plot_drosophila,
                           diff_mean_plot_drosophila,
                           variance_plot_drosophila,
                           diff_variance_plot_drosophila,
                           nrow = 2,
                           ncol = 2,
                           align = "hv")
plot_drosoph

plot2_drosoph <- cowplot::plot_grid( library_size_drosophila,
                             diff_library_size_drosophila,
                             mean_variance_plot_drosophila,
                             diff_mean_variance_plot_drosophila,
                             nrow = 2,
                             ncol = 2,
                             align = "hv")
plot2_drosoph

# Zero plots

diff_mean_zeros_plot_drosophila <- diff_drosoph$Plots$MeanZeros
diff_mean_zeros_plot_drosophila <- diff_mean_zeros_plot_drosophila + ggtitle(label = 'Drosophila 10X: Difference in Mean-Sparsity Relationship')

diff_zero_cell_plot_drosophila <- diff_drosoph$Plots$ZerosCell
diff_zero_cell_plot_drosophila <- diff_zero_cell_plot_drosophila + ggtitle(label = 'Drosophila 10X: Difference in sparsity by cell')

diff_zero_gene_plot_drosophila <- diff_drosoph$Plots$ZerosGene
diff_zero_gene_plot_drosophila <- diff_zero_gene_plot_drosophila + ggtitle(label = 'Drosophila 10X: Difference in sparsity by gene')
diff_zero_gene_plot_drosophila

zero_gene_plot_drosophila <- comparison_drosoph$Plots$ZerosGene
zero_gene_plot_drosophila <- zero_gene_plot_drosophila + scale_y_continuous(name = 'Sparsity [%]', limits = c(75, 100)) + ggtitle(label = 'Drosophila 10X: Sparsity by gene')

zero_cell_plot_drosophila <- comparison_drosoph$Plots$ZerosCell
zero_cell_plot_drosophila<- zero_cell_plot_drosophila + scale_y_continuous(name = 'Sparsity [%]', limits = c(70, 100)) + ggtitle(label = 'Drosophila 10X: Sparsity by cell')

mean_zeros_plot_drosophila <- comparison_drosoph$Plots$MeanZeros
mean_zeros_plot_drosophila <- mean_zeros_plot_drosophila + ggtitle(label = 'Drosophila 10X: Mean-Sparsity Relationship')


plot3_drosoph <- cowplot::plot_grid(zero_gene_plot_drosophila,
                            diff_zero_gene_plot_drosophila,
                            zero_cell_plot_drosophila,
                            diff_zero_cell_plot_drosophila,
                            mean_zeros_plot_drosophila,
                            diff_mean_zeros_plot_drosophila,
                            nrow = 3,
                            ncol = 2,
                            align = "hv")

plot3_drosoph
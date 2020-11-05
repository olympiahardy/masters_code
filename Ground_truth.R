# Extract list of upregulated DE genes from scDesign object
genesUp <- scDesign_drosophila[[2]]
genesUp <- genesUp[[2]]

# Extract list of downregulated DE genes from scDesign object
genesDown <- scDesign_drosophila[[3]]    
genesDown <- genesDown[[2]]

# Create list of non-DE genes
DE_genes <- as.character(append(genesUp, genesDown))


genes <- row.names(drosophila)
nonDE_genes <- setdiff(genes, DE_genes)

# Create ground truth list containing DE genes and non-DE genes
ground_truth <- list(
  DE = as.data.frame(unlist(DE_genes)), 
  nonDE = as.data.frame(unlist(nonDE_genes))
)

# Creating a binary annotation of DE/non-DE genes 
genes <- as.data.frame(row.names(drosophila))
genes$DE_State <- rep(0, times = length(genes))
genes$DE_State[genes$`row.names(drosophila)` %in% ground_truth$DE$`unlist(DE_genes)`] = 1
rownames(genes) <- genes[,1] 
genes[,1] <- NULL
DE_State_annotation <- genes
names(DE_State_annotation)[1] <- "differential.expression"

#############  Part I:  Clustering  ###############

## Libraries
library(Seurat)  #ver 2.3.4
library(dplyr)  #ver 0.8.1
library(ape) #ver 5.4-1
library(pvclust) #ver 2.2-0
library(RColorBrewer)  #ver 1.1-2
library(goseq)  #ver 1.40.0
library(cowplot)  #ver 1.1.0
library(wgcna)  #ver 1.70-3
library(ggtree)  #v2.2.4
library(treeio)  #v1.12.0
library(tidytree)   #v0.3.3
library(phytools) #v0.7-70
library(mvMORPH)  #v1.1.4
library(gridExtra)  #v2.3
library(data.table)  #v1.13.2
library(parallel)  #v4.0.3
library(future)  #v1.20.1

#import 10x matrix and sample metadata
sl_count_matrix <- read.delim(file = "data/GSE134912_spongilla_10x_count_matrix.txt")
sampleID_metadata <- read.delim(file = "data/sample_ID_metadata.tsv", stringsAsFactors = F)

#gene names are trimmed to 110 characters for purposes of plotting
rownames(sl_count_matrix) <- strtrim(rownames(sl_count_matrix), width = 110)

#initialize seurat object
sl <- CreateSeuratObject(raw.data = sl_count_matrix, project = "spongilla_sc")
sl <- AddMetaData(object = sl, metadata = sampleID_metadata)

#normalization - log counts per ten thousand
sl <- NormalizeData(object = sl, normalization.method = "LogNormalize", scale.factor = 10000)

#Identify variable genes following the method of Macosko et al. 2015
sl <- FindVariableGenes(object = sl, mean.function = ExpMean, dispersion.function = LogVMR, do.plot = F, x.low.cutoff = .025, y.cutoff = 1)

#scale and center
sl <- ScaleData(object = sl)

#principal components analysis of variable genes
sl <- RunPCA(object = sl, pc.genes = sl@var.genes, do.print = FALSE, pcs.print = 1:50, pcs.compute = 100)

#tsne plots
sl <- RunTSNE(object = sl, dims.use = 1:40, do.fast = TRUE)

#manuscript tsne embeddings - included here for full reproducibility
manuscript_tsne_embeddings <- as.matrix(read.delim(file = "./data/tsne_cell_embedding.tsv"))
sl <- SetDimReduction(sl, reduction.type = "manuscript.tsne", slot = "cell.embeddings", new.data = manuscript_tsne_embeddings)
sl <- SetDimReduction(sl, reduction.type = "manuscript.tsne", slot = "key", new.data = "manuscript.tsne")
DimPlot(object = sl, reduction.use = "manuscript.tsne", pt.size = 0.5)

## High resolution Louvain clustering
sl <- FindClusters(object = sl, reduction.type = "pca", dims.use = 1:40, resolution = 10, save.SNN = TRUE, print.output = F)


####### Merging Clusters ############

## The merge_clusters function iteratively merges nearest neighbor clusters 
## if they do not satisfy a given threshold of difference. For spongilla single-cell RNAseq
## we required that for nearest neighbor clusters to be maintained as separate clusters
## there must be at least 20 genes differentially expressed at 2 fold difference either higher
## or lower between the two clusters.

merge_clusters <- function(object, num.genes = 20, fold.difference = 2, adj.p.value = 0.05) {
  # first while loop lumps clusters with 2 or fewer cells with their nearest neighbor cluster in expression space.
  # This initial step is done since differential expression tests on clusters with few cells is not possible.
  if(min(table(object@ident)) < 3){
    counter = 1
    while(counter > 0){
      print("Merging clusters with 2 or fewer cells into nearest neighbor cluster")
      counter = 0 # counter to keep track of clusters being merged
      averages <- AverageExpression(object, show.progress = F)  # calculate avg. expression profiles for each cluster
      merged_clusters <- as.character(object@ident) # create new identity vector that will be updated when merging clusters
      correlations <- cor(averages, averages)  # calculate pairwise pearson correlation matrix between all clusters
      correlations[lower.tri(correlations,diag=TRUE)] <- NA  #convert lower triangle and diagonal to NAs
      correlations <- as.data.frame(as.table(correlations))  #convert to 3-column table
      correlations <- na.omit(correlations)  #Get rid of the NAs from lower triangle (duplicates)
      correlations <- correlations[order(-abs(correlations$Freq)),]  #Sort pairwise correlations highest to lowest
      
      for(i in levels(object@ident)){
        if(length(which(object@ident == i)) < 3){ # test whether cluster has fewer than 3 cells 
          index <- which(correlations == i, arr.ind = T)[1,] # identify nearest neighbor cluster in expression space using correlation table
          closest_cluster <- correlations[index[1],]
          new_cluster = paste(as.character(closest_cluster[1,1]),as.character(closest_cluster[1,2]),sep = ".") # create new merged cluster identity
          merged_clusters[which(merged_clusters == closest_cluster[1,1])] <- new_cluster # assign cells from old cluster to new merged cluster identity
          merged_clusters[which(merged_clusters == closest_cluster[1,2])] <- new_cluster
          counter = counter + 1
        }
      }
      merged_clusters <- as.data.frame(merged_clusters)  # add merged_clusters identity to seurat object 
      row.names(merged_clusters) <- colnames(object@data) 
      object <- AddMetaData(object = object, metadata = merged_clusters, col.name = "merged_clusters")
      object <- SetAllIdent(object = object, id = "merged_clusters")
      print(paste(as.character(counter), " clusters with 1 or 2 cells merged with nearest neighbor cluster"), sep = "")
    }
  }
  
  #LOOP 2 iteratively merges clusters that do not satisfy minimal threshold of differentially expressed genes at 2 fold difference
  counter = 1
  while (counter > 0){
    print(paste("Merging nearest-neighbor clusters with less than ", as.character(num.genes), " genes with log fold difference of ", as.character(fold.difference), ".", sep = ""))
    averages <- AverageExpression(object = object, show.progress = F) # calculate avg. expression profiles for each cluster
    merged_clusters <- as.character(object@ident) # create new identity vector that will be updated when merging clusters
    counter = 0 # reset counter
    correlations <- cor(averages, averages) # calculate pairwise pearson correlation matrix between all clusters
    correlations[lower.tri(correlations,diag=TRUE)] <- NA  #convert lower triangle and diagonal to NAs
    correlations <- as.data.frame(as.table(correlations))  #Turn into a 3-column table
    correlations=na.omit(correlations)  #Get rid of the junk we flagged above
    correlations <- correlations[order(-abs(correlations$Freq)),]    #Sort by highest correlation (whether +ve or -ve)
    tested = c("first") # create vector to keep track of which clusters have been tested for merging in this iteration of the loop
    
    for (i in 1:nrow(correlations)){  #iterate through pairwise correlation table starting with pair of clusters that have highest pearson correlation between their expression vector
      if (!(correlations[i,1] %in% tested) & !(correlations[i,2] %in% tested)){ # Verify that neither cluster has been tested for merging already during this iteration of the loop
        if(length(which(object@ident == as.character(correlations[i,1])))>2 & length(which(object@ident == as.character(correlations[i,2])))>2){  # verify each cluster has more than 2 cells
          tested <- c(tested, as.character(correlations[i,1]), as.character(correlations[i,2])) # Add clusters to the list of clusters that have been tested for merging during this loop iteration
          diff_markers <- FindMarkers(object = object, ident.1 = correlations[i,1], ident.2 = correlations[i,2], print.bar = F, min.pct = 0.25) # conduct differential expression test
          diff_markers <- tibble::rownames_to_column(diff_markers, var = "gene")
          diff_markers_filtered <- diff_markers %>% filter(avg_logFC > log(fold.difference), p_val_adj < adj.p.value) # filter differential marker table for genes satisfying threshold cluster 1 > cluster 2
          diff_markers_filtered2 <- diff_markers %>% filter(avg_logFC < log(1/fold.difference), p_val_adj < adj.p.value) # filter differential marker table for genes satisfying threshold cluster 2 > cluster 1
          
          if(nrow(diff_markers_filtered) < num.genes & nrow(diff_markers_filtered2) < num.genes){ # test whether there are minimum # of genes satisfying differentialexpression test
            new_cluster = paste(correlations[i,1], correlations[i,2], sep = ".") # create new merged cluster identity
            merged_clusters[which(merged_clusters == correlations[i,1])] <- new_cluster # reassign cells from old cluster to new merged cluster identity
            merged_clusters[which(merged_clusters == correlations[i,2])] <- new_cluster # reassign cells from old cluster to new merged cluster identity
            counter = counter + 1
          }
        }
      }
    }
    merged_clusters <- as.data.frame(merged_clusters) # add new cluster identities in metadata slot of seurat object
    row.names(merged_clusters) <- colnames(object@data)
    object <- AddMetaData(object = object, metadata = merged_clusters, col.name = "merged_clusters")
    object <- SetAllIdent(object = object, id = "merged_clusters")
    print(paste(as.character(counter), " clusters merged with nearest neighbor cluster"), sep = "")
  }
  return(object)
}

sl <- merge_clusters(object = sl, num.genes = 20, fold.difference = 2, adj.p.value = 0.05)

DimPlot(object = sl, reduction.use = "manuscript.tsne", pt.size = 0.5, do.label = T)


# rename clusters to reflect those used in manuscript
identities <- levels(sl@ident)
manuscript_clusters <- c("2", "8", "3", "5", "16", "1", "7", "6", "4", "10", "11", "12", "20", "27", "25", "13", "17", "14", "33", "9", "26", "15", "31", "38", "42", "24", "29", "19", "30", "22", "23", "28", "40", "36", "21", "18", "39", "35", "32", "34", "37", "41")
cell_type_assignments <- plyr::mapvalues(sl@ident, identities, manuscript_clusters)
sl <- AddMetaData(object = sl, metadata = cell_type_assignments, col.name = "manuscript_cluster_numbers")
sl <- SetAllIdent(object = sl, id =  "manuscript_cluster_numbers")
DimPlot(object = sl, reduction.use = "manuscript.tsne", pt.size = 0.5, do.label = T)


# identify marker genes
markers <- FindAllMarkers(object = sl, only.pos = TRUE, logfc.thrheshold = 0.25, test.use = "wilcox", min.pct = 0.1)


## Distinguish between clusters with uniquely expressed genes and transitional clusters that lack unique markers
averages <- AverageExpression(object = sl) # calculate avg. expression matrix for each cluster
max_exp_index <- apply(averages, MARGIN = 1, FUN = which.max)  # create index vector of cluster with maximum expression for each gene
max_exp_cluster <- colnames(averages)[max_exp_index] # create vector with the cluster that has highest expression for each gene
names(max_exp_cluster) <- row.names(averages) # name each element of previous vector with gene name

# filter marker table to only contain "unique" markers, i.e. those markers that are most highly expressed in each cluster
unique_markers <- markers %>% 
  filter(p_val_adj < 0.5, avg_logFC > log(2), cluster == max_exp_cluster[gene])

# create table with number of unique markers for each cluster
marker_count <- as.data.frame(table(unique_markers$cluster))

# assign transitional status to cell types with less than 30 unique markers
transitional_assignment <- as.character(sl@ident)
transitional_clusters <- as.character(marker_count$Var1[marker_count$Freq < 30])
for (x in transitional_clusters){
  transitional_assignment[transitional_assignment == x] <- "transitional"
  transitional_assignment <- as.data.frame(transitional_assignment)
  row.names(transitional_assignment) <- names(sl@ident)
}

sl <- AddMetaData(object = sl, metadata = transitional_assignment)
sl <- SetAllIdent(object = sl, id = "transitional_assignment")
DimPlot(object = sl, reduction.use = "manuscript.tsne", pt.size = 0.5, do.label = T)



###### WGCNA gene set analysis ######

# We first identify an expanded set of variable genes to use for identify gene sets.
# NOTE: although not shown here, we comprehensively tested the effect of using different variable gene thresholds for 
# identifying gene sets. These generally yielded similar results.

# Find variable genes
variable_genes2 <- FindVariableGenes(object = sl, do.plot = F, set.var.genes = F, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0)

#prep cell type expression matrix
celltype_table <- read_tsv("./data/Table1_celltype_descriptions_final.tsv", col_types = "ccccc")

#subset seurat object to only include distinct cell type clusters
sl_types <- SubsetData(object = sl, ident.use = celltype_table$`Cluster #`)


#assign cell type names used in the manuscript
celltype_names <- Idents(sl_types_test)
celltype_names2 <- celltype_table$Nomenclature[match(celltype_names, celltype_table$`Cluster #`)]
sl_types_averages <- AddMetaData(object = sl_types, metadata = celltype_names2, col.name = "cell_type")
Idents(sl_types_averages) <- "cell_type"
levels(x = sl_types_averages@idents) <- c("incPin1","incPin2","apnPin1","apnPin2","Lph","basPin","Spc","Met1","Met2","Chb1","Chb2","Cho","Apo","Myp1","Myp2","Amb","Grl","Nrd","Mes1","Mes2","Mes3","Arc","Scl")

#calculate average expression matrices
sl_types_averages <- AverageExpression(object = sl_types, return.seurat = T)
sl_types_average_matrix <- AverageExpression(object = sl_types, return.seurat = F)

#matrix construction for wgcna:
matrix_wgcna <- t(sl_types_averages@data[variable_genes2,])


#Choose soft thresholding:
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(matrix_wgcna, powerVector = powers, verbose = 5)#, corFnc = bicor, corOptions = list(use= 'p', maxPOutliers = 0.05))
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#set soft power:
softPower = 5

#calculate adjacency matrix
adjacency = adjacency(matrix_wgcna, power = softPower)

# Turn adjacency into topological overlap to minimize affects of noise and spurios associations
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)

# Set minimum module size
minModuleSize = 30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                  deepSplit = 4, pamRespectsDendro = FALSE,
                                  minClusterSize = minModuleSize)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
length(table(dynamicColors))
table(dynamicColors)

# Plot the dendrogram and colors underneath
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList = moduleEigengenes(as.matrix(matrix_wgcna), colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
MEDissThres = 0.1
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(as.matrix(matrix_wgcna), dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
length(table(mergedColors))
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

plotDendroAndColors(geneTree, main = "Gene dendrogram with pre- and post-merged module colors", cbind(dynamicColors, mergedColors),c("Dynamic Tree Cut", "Merged dynamic"),dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)




##### Cell Type Phylogeny #####

# In general, we noticed the expression of archeocyte markers in differentiated cell type clusters is highly dependent
# on where the exact boundaries of the differentiated cluster are drawn. Based on this, we remove archaeocyte markers from
# our phylogenetic analysis (NOTE: in practice we get highly similar results whether or not we include archaeocyte genes or not).

#Identify Archeocyte genes
gene_max_clusters_types <- apply(sl_types_average_matrix, MARGIN = 1, FUN = which.max)
gene_max_clusters_types <- colnames(sl_types_average_matrix)[gene_max_clusters_types]
names(gene_max_clusters_types) <- row.names(sl_types_average_matrix)

genes_highest_in_archeocytes <- names(gene_max_clusters_types[which(gene_max_clusters_types == "Ar")])
genes_not_highest_in_archeocytes <- row.names(sl@data)[which(!(row.names(sl@data) %in% genes_highest_in_archeocytes))]

# Next, we identify an expanded set of variable genes to use for phylogeny inference. We tested robustness of the topology
# using a variety of thresholds for variable genes (variable gene set sizes of approximately 2500, 6000, and 13000 genes), \#
# as well as constructing the tree with only transcription factors.

#### Identifying variable gene sets

variable_genes <- FindVariableGenes(object = sl, do.plot = F, set.var.genes = F, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = -1)
variable_genes2 <- FindVariableGenes(object = sl, do.plot = F, set.var.genes = F, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0)
variable_genes3 <- FindVariableGenes(object = sl, do.plot = F, set.var.genes = F, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = .025, y.cutoff = 0.5)


## Tree construction - expanded variable gene set - n=11178 variable gene set
variable_genes_not_highest_in_archeocytes <- variable_genes[which(variable_genes %in% genes_not_highest_in_archeocytes)]
vargenes_noarch_matrix <- t(types_average_norm[variable_genes_not_highest_in_archeocytes,])
vargenes_noarch_matrix_dist <- dist(vargenes_noarch_matrix, method = "euclidean")
vargenes_noarch_matrix_tree <- nj(vargenes_noarch_matrix_dist)
plot.phylo(vargenes_noarch_matrix_tree, type = "unrooted")
plot.phylo(midpoint.root(vargenes_noarch_matrix_tree))

vargenes_noarch_matrix_tree <-midpoint.root(vargenes_noarch_matrix_tree)
vargenes_noarch_matrix_tree_boot <- boot.phylo(vargenes_noarch_matrix_tree, vargenes2_noarch_matrix, FUN = function(xx) nj(dist(xx)), B = 10000)
vargenes_noarch_matrix_tree$node.label <- vargenes_noarch_matrix_tree_boot/100
vargenes_noarch_matrix_tree$node.label <- round(vargenes_noarch_matrix_tree$node.label)

## Tree construction - wgcna variable genes - n=5945 variable gene set
variable_genes2_not_highest_in_archeocytes <- variable_genes2[which(variable_genes2 %in% genes_not_highest_in_archeocytes)]
types_average_norm <- as.matrix(sl_types_averages@data)
vargenes2_noarch_matrix <- t(as.matrix(types_average_norm)[variable_genes2_not_highest_in_archeocytes,])
vargenes2_noarch_matrix_dist <- dist(vargenes2_noarch_matrix, method = "euclidean")
vargenes2_noarch_matrix_tree <- nj(vargenes2_noarch_matrix_dist)
plot.phylo(vargenes2_noarch_matrix_tree, type = "unrooted")
plot.phylo(midpoint.root(vargenes2_noarch_matrix_tree))

vargenes2_noarch_matrix_tree <-midpoint.root(vargenes2_noarch_matrix_tree)
vargenes2_noarch_matrix_tree_boot <- boot.phylo(vargenes2_noarch_matrix_tree, vargenes2_noarch_matrix, FUN = function(xx) nj(dist(xx)), B = 10000)
vargenes2_noarch_matrix_tree$node.label <- vargenes2_noarch_matrix_tree_boot/100
vargenes2_noarch_matrix_tree$node.label <- round(vargenes2_noarch_matrix_tree$node.label)
write.tree(vargenes2_noarch_matrix_tree, file = "./data/vargenes2_noarch_matrix_tree.nwk")

## Tree construction - n=2454 variable gene set
vargenes3_matrix <- t(types_average_norm[variable_genes3,])
vargenes3_matrix_dist <- dist(vargenes3_matrix, method = "euclidean")
vargenes3_matrix_tree <- nj(vargenes3_matrix_dist)
plot.phylo(vargenes3_matrix_tree, type = "unrooted")
plot.phylo(midpoint.root(vargenes3_matrix_tree))

vargenes3_matrix_tree <-midpoint.root(vargenes3_matrix_tree)
vargenes3_matrix_tree_boot <- boot.phylo(vargenes3_matrix_tree, vargenes3_matrix, FUN = function(xx) nj(dist(xx)), B = 10000)
vargenes3_matrix_tree$node.label <- vargenes3_matrix_tree_boot/100
vargenes3_matrix_tree$node.label <- round(vargenes3_matrix_tree$node.label)


######## Modeling gene expression on tree with Ornstein-Uhlenbeck models #######
# This script fits univariate and multivariate Ornstein-Uhlenbeck models to the cell type tree 
# to determine genes that change their expression optima on the tree.


## setup parallelization
plan("multiprocess", workers = 4)
options(future.globals.maxSize= 2000*1024^2)


## Functions
OUtest <- function(gene, exp_matrix, simmap_tree){
  gene_exp <- exp_matrix[gene,]
  try_gene_OU1 <- try(mvOU(simmap_tree, gene_exp, model="OU1", diagnostic=FALSE, 
                           echo=FALSE, param = list(root = F)), silent = T)
  try_gene_OUM <- try(mvOU(simmap_tree, gene_exp, model="OUM", diagnostic=FALSE, 
                           echo=FALSE, param = list(root = F)), silent = T)
  outlist <- fifelse(!is(try_gene_OU1, "try-error") & !is(try_gene_OUM, "try-error"), 
                     return(list(gene_OU1_convergence = try_gene_OU1$convergence, 
                                 gene_OUM_convergence = try_gene_OUM$convergence, 
                                 gene_OU1_hess = try_gene_OU1$hess.values, 
                                 gene_OUM_hess = try_gene_OUM$hess.values, 
                                 lrt_pval = LRT(try_gene_OUM, try_gene_OU1, echo = F)$pval, 
                                 OUM_theta2_max = which.max(as.numeric(try_gene_OUM$theta)) == 2,
                                 AIC_OU1 = AIC(try_gene_OU1), AIC_OUM = AIC(try_gene_OUM))), 
                     return(NA))
}

getnodeID <- function(tree, label){
  tree_tibble <- as_tibble(tree)
  x <- as.integer(tree_tibble$node[which(tree_tibble$label == label)])
  return(x)
}




#### Perform OU tests

## Create expression matrix of counts per ten thousand (CPT) that only contains genes expressed > .03 CPT in at least one cell type
## OU tests will only be performed for these "expressed" genes
gene_max_expression <- apply(sl_types_average_matrix, 1, FUN=max)
sl_types_average_matrix_expressed <- sl_types_average_matrix[gene_max_expression > 0.03,]

## set node and tree
nodes = as_tibble(vargenes2_noarch_matrix_tree) %>%
  filter(!is.na(branch.length)) %>%
  pull(node)

tree = "vargenes2_noarch_matrix_tree"
exp_matrix = sl_types_average_matrix_expressed[,get(tree)$tip.label]
nb_cores = 10

## Loop through each node in the tree
for(i in nodes){
  
  if(isTip(vargenes2_noarch_matrix_tree, i)){
    simmap_tree <- paintBranches(get(tree),edge = i, state = "high", anc.state = "low")
  } else {
    simmap_tree <- paintSubTree(get(tree), node = i, state = "high", anc.state = "low")
  }
  
  print(paste("node ", i, "    ", Sys.time(), sep = ""))
  start_time <- Sys.time()
  OU_genes_out <- mclapply(rownames(exp_matrix), OUtest, exp_matrix = exp_matrix, simmap_tree = simmap_tree, 
                           mc.cores = getOption("mc.cores", nb_cores))
  end_time <- Sys.time()
  print(end_time - start_time)
  
  OU_genes_out2 <- OU_genes_out %>%
    map_dfr(as_tibble) %>%
    rownames_to_column(var = "rownumber") %>%
    mutate(gene = rownames(exp_matrix)[as.integer(rownumber)]) %>%
    add_column(node = i, tree = tree) %>%
    mutate(pval.adj = p.adjust(lrt_pval, method='BH')) %>%
    select(gene, lrt_pval, pval.adj, node, OUM_theta2_max, AIC_OU1, AIC_OUM, gene_OU1_convergence, gene_OUM_convergence, 
           gene_OU1_hess, gene_OUM_hess, tree)
  
  if (exists("OU_genes_out_final")){
    OU_genes_out_final <- bind_rows(OU_genes_out_final, OU_genes_out2)
  } else {
    OU_genes_out_final <- OU_genes_out2
  }
}

# Summarize number of changes along each branch
OU_genes_out_final_node_summary <- OU_genes_out_final %>%
  filter(pval.adj < 0.05) %>%
  group_by(node) %>%
  summarize(num_genes = n()) %>%
  mutate(num_genes_log10 = log10(num_genes))


#### Gene ontology enrichment analysis of OU gene sets (for each clade)

## create table with annotations for each clade in the cell type tree
clade_labels <- tibble(node = c(44, 43, 42, 41, 40, 39, 45, 38, 35, 34, 33, 32, 31, 37, 36, 30, 28, 27, 29, 26, 25, 17, 18, 19, 20, 15, 21, 16, 22, 23, 5, 6, 3, 4, 1, 2, 7, 8, 9, 12, 13, 14, 10, 11),
                       clade_label = c("incPin1/2", "apnPin1/2", "incPin_apnPin", "incPin_apnPin_Lph", "Scp_basPin", "incPin_apnPin_Lph_Scp_basPin", "Met", "Endymocytes", "Chb", "Cho_Apo", "Chb_Cho_Apo", "Myp", "Peptidocytes", "Amb_Grl", "Amoeboid-Neuroid", "Peptidocyte_Amb-Nrd", "Mes1_Mes2", "Mes1_Mes2_Mes3", "Arc_Scl", "Archeocytes_and_relatives", "Peptidocyte_Amb-Nrd_Arc_and_rel", "incPin1", "incPin2", "apnPin1", "apnPin2", "Lph", "Scp", "basPin", "Met1", "Met2", "Chb1", "Chb2", "Cho", "Apo", "Myp1", "Myp2", "Amb", "Grl", "Nrd", "Mes1", "Mes2", "Mes3", "Arc", "Scl"))

## Import and prep gene ontology files
spongilla_genes <- read.delim(file = "./data/2017_12_07_merged_gene_lengths_no_mito_or_ERCCs.txt", stringsAsFactors = F)
rownames(spongilla_genes) <- strtrim(rownames(spongilla_genes), width = 110)
spongilla.go <-  read.delim(file= "./data/2020_08_13_merged_trinity_gene_gene_ontology_assignments_jake_manual_additions2_forsl2_final2.txt", stringsAsFactors = F)
spongilla.go$spongilla_gene <- strtrim(spongilla.go$spongilla_gene, width = 110)
assayed.genes <- rownames(sl_types_average_matrix_expressed)
goterms_catalog <- as_tibble(spongilla.go) 
goterms_catalog <- goterms_catalog %>%
  dplyr::select(go_id, term) %>%
  unique
spongilla_go <- as_tibble(spongilla.go) %>%
  mutate(spongilla_gene2 = gsub("_", "-", spongilla_gene)) %>%
  filter(spongilla_gene %in% rownames(sl_types_average_matrix_expressed))


## Loop through each branch and Run GO analysis
rm(enriched_go_allnodes)
for(i in 1:nrow(clade_labels)){
  
  rm(enriched_go)
  node_number <- as.integer(clade_labels$node[i])
  clade_label <- clade_labels$clade_label[i]
  
  gene_list <- OU_genes_out_final %>%
    dplyr::filter(node == node_number, pval.adj < 0.05, OUM_theta2_max == T) %>%
    dplyr::arrange(pval.adj) %>%
    pull(gene)
  
  gene.vector <- as.integer(assayed.genes%in%gene_list)
  names(gene.vector) <- assayed.genes
  pwf <- nullp(gene.vector, bias.data=spongilla.lengths, plot.fit = F)
  go <- goseq(pwf, gene2cat = spongilla.go, method = "Hypergeometric")
  go$goterm_name <- goterms_catalog$term[match(go$category, goterms_catalog$go_id)]
  
  go$clade_label <- rep(clade_label, nrow(go))
  enriched_go <- as_tibble(go) %>%
    filter(over_represented_pvalue < 0.01)
  
  if(!exists("enriched_go_allnodes")){
    enriched_go_allnodes <- enriched_go
  } else {
    enriched_go_allnodes <- bind_rows(enriched_go_allnodes, enriched_go)
  }
  
  enriched_go_allnodes <- select(enriched_go_allnodes, "category", "over_represented_pvalue", "under_represented_pvalue", "numDEInCat","numInCat","goterm_name","ontology","clade_label")
}





####### Find markers for cell type clusters (excluding transitional clusters) #######

## This was performed using an updated version of Seurat
library(Seurat) #ver 3.2.2
sl_types_v3 <- UpdateSeuratObject(sl_types)


## First identify differentially expressed genes for each distinct cell type
## As background we use all expressed genes.

sl_types_v3_expressed <- subset(sl_types_v3, features = gsub("_", "-", rownames(sl_types_average_matrix_expressed)))
sl_types_markers_expressed <- FindAllMarkers(sl_types_v3_expressed, logfc.threshold = 0, min.pct = 0, features = rownames(sl_types_v3_expressed))
sl_types_markers_expressed_tibble2 <- as_tibble(sl_types_markers_expressed) %>%
  group_by(cluster) %>%
  mutate(BH_fdr = p.adjust(p_val, method = "BH", n = length(rownames(sl_types_v3_expressed)))) %>%
  ungroup()

sl_types_markers_expressed_tibble2$cluster <- ordered(sl_types_markers_expressed_tibble2$cluster, levels = c("incPin1","incPin2","apnPin1","apnPin2","Lph","basPin","Scp","Met1","Met2","Chb1","Chb2","Cho","Apo","Myp1","Myp2","Amb","Grl","Nrd","Mes1","Mes2","Mes3","Arc","Scl"))


## Run GO analysis for differentially expressed gene sets

# update GO files for Seurat version 3 gene names (which replace "_" with "-" in gene names)
assayed.genes2 <- gsub("_", "-", assayed.genes)
spongilla.lengths2 <- spongilla.lengths
names(spongilla.lengths2) <- gsub("_", "-", names(spongilla.lengths))
spongilla.go2 <- spongilla.go
spongilla.go2$spongilla_gene <- gsub("_", "-", spongilla.go$spongilla_gene)

# loop through cell types and determine enriched terms for each
rm(enriched_go_celltypes)
for(i in levels(as.factor(sl_types_markers_expressed_tibble2$cluster))){
  
  rm(enriched_go)
  cell_type <- i
  
  gene_list <- sl_types_markers_expressed_tibble2 %>%
    dplyr::filter(cluster == cell_type, BH_fdr < 0.05, avg_logFC > 0) %>%
    dplyr::arrange(p_val) %>%
    pull(gene)
  
  gene.vector <- as.integer(assayed.genes2%in%gene_list)
  names(gene.vector) <- assayed.genes2
  pwf <- nullp(gene.vector, bias.data=spongilla.lengths2, plot.fit = F)
  go <- goseq(pwf, gene2cat = spongilla.go2, method = "Hypergeometric")
  go$goterm_name <- goterms_catalog$term[match(go$category, goterms_catalog$go_id)]
  
  enriched_go <- as_tibble(go) %>%
    filter(over_represented_pvalue < 0.01)
  
  enriched_go$cell_type <- rep(cell_type, nrow(enriched_go))
  
  if(!exists("enriched_go_celltypes")){
    enriched_go_celltypes <- enriched_go
  } else {
    enriched_go_celltypes <- bind_rows(enriched_go_celltypes, enriched_go)
  }
}

library(RegulationZIndex)

my_data <- readRDS("Res29samp_integration_wHippo_alpha_beta_epsilon_wExpCond8_rmOutlier_wAnno.rds")
group_labels <- my_data@active.ident
counts <- my_data@assays[["RNA"]]@counts
donor.info <- my_data@meta.data[["orig.ident"]]
donors <- c("HI_1", "HI_11", "HI_12", "HI_14", "HI_15", "HI_16", "HI_17", "HI_18", "HI_19", "HI_2", "HI_20", "HI_21", "HI_22", "HI_23", "HI_24",
            "HI_25", "HI_26", "HI_27", "HI_28", "HI_3", "HI_31", "HI_32", "HI_4", "HI_6", "HI_7", "HI10", "HI5", "HI8", "HI9")

sample.info <- cbind(pheno = my_data@meta.data[["expCond1"]], subject = donor.info)
sample.info <- sample.info[!duplicated(sample.info), ]
rownames(sample.info) <- sample.info[, 1]


list.celltype <- c("Immune", "Quiescent stellate", "Activated stellate", "Endothelial", "Ductal", "Acinar", "Epsilon",
                   "Delta", "Gamma", "Beta6", "Beta5", "Beta4", "Beta3", "Beta2", "Beta1", "Alpha")
list.celltype_br <- c("Immune", "QS", "AS", "Endo", "Ductal", "Acinar", "Epsilon",
                      "Delta", "Gamma", "Beta6", "Beta5", "Beta4", "Beta3", "Beta2", "Beta1", "Alpha")

for(i in c(1:6, 8:16)){

  selected <- which(group_labels %in% list.celltype[i])
  testdata <- counts[, selected]

  sub.donor.info <- donor.info[group_labels %in% list.celltype[i]]
  new.clusters <- as.numeric(as.factor(sub.donor.info))

  sample.info <- cbind(new.clusters, sub.donor.info)
  sample.info <- sample.info[!duplicated(sample.info), ]

  cluster.zscore.summary <- NB_inflation_test_var_asymp(testdata, grouping = new.clusters)

  save(sample.info, cluster.zscore.summary, file = paste0("scRNAseq_", list.celltype_br[i], "_score_donors.rds"))
}



library(RegulationZIndex)

my_data <- readRDS("Res29samp_integration_wHippo_alpha_beta_epsilon_wExpCond8_rmOutlier_wAnno.rds")
group_labels <- my_data@active.ident
counts <- my_data@assays[["RNA"]]@counts
donor.info <- my_data@meta.data[["orig.ident"]]
donors <- c("HI_1", "HI_11", "HI_12", "HI_14", "HI_15", "HI_16", "HI_17", "HI_18", "HI_19", "HI_2", "HI_20", "HI_21", "HI_22", "HI_23", "HI_24",
            "HI_25", "HI_26", "HI_27", "HI_28", "HI_3", "HI_31", "HI_32", "HI_4", "HI_6", "HI_7", "HI10", "HI5", "HI8", "HI9")

sample.info <- cbind(pheno = my_data@meta.data[["expCond1"]], subject = donor.info)
sample.info <- sample.info[!duplicated(sample.info), ]
rownames(sample.info) <- sample.info[, 1]

pheno <- cbind(my_data@meta.data[["orig.ident"]], my_data@meta.data[["expCond1"]])
pheno <- pheno[!duplicated(pheno), ]
pheno[pheno[, 2]%in%"OW-ND", 2] <- "O1"
pheno[pheno[, 2]%in%"OB-ND", 2] <- "O2"
rownames(pheno) <- pheno[, 1]


list.celltype_br <- c("Immune", "QS", "AS", "Endo", "Ductal", "Acinar", "Epsilon",
                      "Delta", "Gamma", "Beta6", "Beta5", "Beta4", "Beta3", "Beta2", "Beta1", "Alpha")
files_res <- paste0("scRNAseq_", list.celltype_br, "_score_donors.rds")

hits.list <- NULL
for(i in c(1:6, 8:16)){
  load(files_res[i])
  kkk <- length(cluster.zscore.summary)
  for(j in 1:kkk){
    ids <- which(cluster.zscore.summary[[j]]$Inflated$gene_mean>1)
    hits.list <- c(hits.list, rownames(cluster.zscore.summary[[j]]$Inflated)[ids])
  }
}
aaa <- table(hits.list)
selected <- names(aaa[aaa>=10])

zscore_summary <- NULL
means_summary <- NULL
for(i in c(1:6, 8:16)){
  load(files_res[i])
  kkk <- length(cluster.zscore.summary)
  rownames(sample.info) <- sample.info[, 1]
  all.zscore <- NULL
  all.names <- NULL
  all.means <- NULL
  for(j in 1:kkk){
    all.names <- c(all.names, sample.info[as.character(cluster.zscore.summary[[j]]$groupID), 2])
    all.zscore <- cbind(all.zscore, cluster.zscore.summary[[j]]$Zscores[selected, "zscore"])
    all.means <- cbind(all.means, cluster.zscore.summary[[j]]$Zscores[selected, "gene_mean"])
  }
  colnames(all.zscore) <- paste0(all.names, ":", list.celltype_br[i])
  colnames(all.means) <- paste0(all.names, ":", list.celltype_br[i])
  zscore_summary <- cbind(zscore_summary, all.zscore)
  means_summary <- cbind(means_summary, all.means)
}
rownames(zscore_summary) <- selected
rownames(means_summary) <- selected
save(means_summary, zscore_summary, selected, file = "top_inflated_markers.rds")



library(ComplexHeatmap)
load("top_inflated_markers.rds")
cell_type <- gsub(".*:", "", colnames(zscore_summary))
phenotype <- pheno[gsub(":.*", "", colnames(zscore_summary)), 2]


zscore_summary <- zscore_summary[, order(phenotype)]
cell_type <- cell_type[order(phenotype)]
phenotype <- phenotype[order(phenotype)]

zscore_summary[is.na(zscore_summary)] <- 0
zscore_summary[zscore_summary < 0] <- 0
zscore_summary[zscore_summary > 10] <- 10

ha = HeatmapAnnotation(
  cell_type = cell_type,
  phenotype = phenotype,
  border = TRUE)
##########################################################
########## heatmap of Z-index in selected gene families ##
##########################################################


anno_markers <- read.csv("curated_markers_by_families.csv")[, 1:3]
anno_markers <- anno_markers[!duplicated(anno_markers[, 3]), ]
rownames(anno_markers) <- anno_markers[, "Gene"]

sub_anno <- anno_markers[anno_markers$Family%in%c("Granins", "Synaptotagmins", "ROBO", "Neuregulins", "Nephrins", "Cadherins"), ]

other_anno <- cbind(Family = "Others", Class = "Others", Gene = c("NEGR1", "SEMA5A", "ROR1", "NTNG1"))

selected_anno_markers <- rbind(sub_anno, other_anno)
rownames(selected_anno_markers) <- selected_anno_markers$Gene
annotate_gene_names <- function(selected_anno_markers, gene_list){
  ID_with_annotation <- gene_list[which(gene_list%in%selected_anno_markers[, "Gene"])]
  return(selected_anno_markers[ID_with_annotation, c("Gene", "Family")])
}

test_gene_anno <- annotate_gene_names(selected_anno_markers, selected)


pdf(paste0("scRNAseq_alldonors_celltypes_zscore_all_families_1_slice.pdf"), width=12, height = 8)
annotated_heatmap(zscore_summary, test_gene_anno, sample_annotation = ha, column_split = factor(cell_type), families = unique(test_gene_anno[, "Family"]), plot_name = "zscore")
dev.off()




stress_markers <- read.csv("curated_response_processes_markers.csv")[, 1:3]
stress_markers <- stress_markers[!duplicated(stress_markers[, 3]), ]
rownames(stress_markers) <- stress_markers[, "Gene"]
colnames(stress_markers) <- c("Family", "Class", "Gene")
anno_markers <- read.csv("curated_markers_by_families.csv")[, 1:3]
anno_markers <- anno_markers[!duplicated(anno_markers[, 3]), ]
rownames(anno_markers) <- anno_markers[, "Gene"]
sub2_anno <- anno_markers[anno_markers$Family%in%c("HSPs", "HLA", "Hemoglobins"), ]
other2_anno <- rbind(stress_markers, sub2_anno)
test2_gene_anno <- annotate_gene_names(other2_anno, selected)


pdf(paste0("scRNAseq_alldonors_celltypes_zscore_all_families_2_slice.pdf"), width=12, height = 5.5)
annotated_heatmap(zscore_summary, test2_gene_anno, sample_annotation = ha, column_split = factor(cell_type), families = unique(test2_gene_anno[, "Family"]), plot_name = "zscore")
dev.off()



anno_markers <- read.csv("curated_markers_by_families.csv")[, 1:3]
anno_markers <- anno_markers[!duplicated(anno_markers[, 3]), ]
rownames(anno_markers) <- anno_markers[, "Gene"]
sub3_anno <- anno_markers[anno_markers$Family%in%c("Regenerating islet-derived proteins", "Neuropeptides", "Receptor ligands"), ]

test3_gene_anno <- annotate_gene_names(sub3_anno, selected)


pdf(paste0("scRNAseq_alldonors_celltypes_zscore_all_families_3_slice.pdf"), width=12, height = 4)
annotated_heatmap(zscore_summary, test3_gene_anno, sample_annotation = ha, column_split = factor(cell_type), families = unique(test3_gene_anno[, "Family"]), plot_name = "zscore")
dev.off()


#####################################################################
########## heatmap of log-mean in selected gene families ############
#####################################################################

load("top_inflated_markers.rds")
means_summary[is.na(means_summary)] <- 0
means_summary <- t(apply(means_summary, 1, function(x){log(x+0.01)}))
#means_summary[means_summary > 20] <- 20
means_summary[means_summary >= 5] <- 5
means_summary[means_summary <= -2 ] <- -2
pdf(paste0("scRNAseq_alldonors_celltypes_mean_all_families_1_slice.pdf"), width=12, height = 8)
annotated_heatmap(means_summary, test_gene_anno, sample_annotation = ha, column_split = factor(cell_type), families = unique(test_gene_anno[, "Family"]), plot_name = "mean")
dev.off()

pdf(paste0("scRNAseq_alldonors_celltypes_mean_all_families_2_slice.pdf"), width=12, height = 5.5)
annotated_heatmap(means_summary, test2_gene_anno, sample_annotation = ha, column_split = factor(cell_type), families = unique(test2_gene_anno[, "Family"]), plot_name = "mean")
dev.off()
pdf(paste0("scRNAseq_alldonors_celltypes_mean_all_families_3_slice.pdf"), width=12, height = 4)
annotated_heatmap(means_summary, test3_gene_anno, sample_annotation = ha, column_split = factor(cell_type), families = unique(test3_gene_anno[, "Family"]), plot_name = "mean")
dev.off()





##########################################
########## enrichment results ############
##########################################
load("top_inflated_markers.rds")

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
# Example list of gene names
geneList <- selected

# Perform GO enrichment analysis using default settings
enrich_result <- enrichGO(gene = geneList, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")


pdf(paste0("enrichment.pdf"), width=5, height = 10)
dotplot(enrich_result, showCategory=30) + ggtitle("dotplot for BP")
dev.off()



########################################
########## Correlation plot ############
########################################
annotate_gene_names <- function(selected_anno_markers, gene_list){
  ID_with_annotation <- gene_list[which(gene_list%in%selected_anno_markers[, "Gene"])]
  return(selected_anno_markers[ID_with_annotation, c("Gene", "Family")])
}

anno_markers <- read.csv("curated_markers_by_families.csv")[, 1:3]
anno_markers <- anno_markers[!duplicated(anno_markers[, 3]), ]
rownames(anno_markers) <- anno_markers[, "Gene"]

sub_anno <- anno_markers[anno_markers$Family%in%c("Granins", "Synaptotagmins", "ROBO", "Neuregulins", "Nephrins", "Cadherins"), ]

other_anno <- cbind(Family = "Others", Class = "Others", Gene = c("NEGR1", "SEMA5A", "ROR1", "NTNG1", "SST"))

selected_anno_markers <- rbind(sub_anno, other_anno)
rownames(selected_anno_markers) <- selected_anno_markers$Gene
annotate_gene_names <- function(selected_anno_markers, gene_list){
  ID_with_annotation <- gene_list[which(gene_list%in%selected_anno_markers[, "Gene"])]
  return(selected_anno_markers[ID_with_annotation, c("Gene", "Family")])
}

test_gene_anno <- annotate_gene_names(selected_anno_markers, selected)


stress_markers <- read.csv("curated_response_processes_markers.csv")[, 1:3]
stress_markers <- stress_markers[!duplicated(stress_markers[, 3]), ]
rownames(stress_markers) <- stress_markers[, "Gene"]
colnames(stress_markers) <- c("Family", "Class", "Gene")
anno_markers <- read.csv("curated_markers_by_families.csv")[, 1:3]
anno_markers <- anno_markers[!duplicated(anno_markers[, 3]), ]
rownames(anno_markers) <- anno_markers[, "Gene"]
sub2_anno <- anno_markers[anno_markers$Family%in%c("HSPs", "HLA", "Hemoglobins"), ]
other2_anno <- rbind(stress_markers, sub2_anno)
test2_gene_anno <- annotate_gene_names(other2_anno, selected)

anno_markers <- read.csv("curated_markers_by_families.csv")[, 1:3]
anno_markers <- anno_markers[!duplicated(anno_markers[, 3]), ]
rownames(anno_markers) <- anno_markers[, "Gene"]
sub3_anno <- anno_markers[anno_markers$Family%in%c("Regenerating islet-derived proteins", "Neuropeptides", "Receptor ligands"), ]

test3_gene_anno <- annotate_gene_names(sub3_anno, selected)



library(corrplot)


selected.wave <- counts[unique(c(test_gene_anno[, 1], test2_gene_anno[, 1], test3_gene_anno[, 1])), group_labels%in%"Beta2" & donor.info%in%"HI_25"]
flag_0 <- which(apply(selected.wave, 1, sum) == 0)
Mcor <- cor(as.matrix(t(selected.wave[-flag_0, ])))
pdf(paste0("scRNAseq_HI_25_beta2_cor_allslice_genes.pdf"), width=12, height = 12)
tes_cor <- corrplot(Mcor, order = 'hclust')
dev.off()



selected.wave <- counts[unique(c(test_gene_anno[, 1], test2_gene_anno[, 1], test3_gene_anno[, 1])), group_labels%in%"Beta2" & donor.info%in%"HI_2"]

Mcor <- cor(as.matrix(t(selected.wave)))
pdf(paste0("scRNAseq_HI_2_beta2_cor_allslice_genes.pdf"), width=12, height = 12)
tes_cor <- corrplot(Mcor, order = 'hclust')
dev.off()



selected.wave <- counts[unique(c(test_gene_anno[, 1], test2_gene_anno[, 1], test3_gene_anno[, 1])), group_labels%in%"Beta2" & donor.info%in%"HI_1"]

Mcor <- cor(as.matrix(t(selected.wave)))
pdf(paste0("scRNAseq_HI_1_beta2_cor_allslice_genes.pdf"), width=12, height = 12)
tes_cor <- corrplot(Mcor, order = 'hclust')
dev.off()



selected.wave <- counts[unique(c(test_gene_anno[, 1], test2_gene_anno[, 1], test3_gene_anno[, 1])), group_labels%in%"Beta2"]

Mcor <- cor(as.matrix(t(selected.wave)))
pdf(paste0("scRNAseq_beta2_alldonors_cor_allslice_genes.pdf"), width=12, height = 12)
tes_cor <- corrplot(Mcor, order = 'hclust')
dev.off()


selected.wave <- counts[unique(c(test_gene_anno[, 1], test2_gene_anno[, 1], test3_gene_anno[, 1])), group_labels%in%"Alpha"]

Mcor <- cor(t(as.matrix(selected.wave)))
pdf(paste0("scRNAseq_alpha_alldonors_cor_allslice_genes.pdf"), width=12, height = 12)
tes_cor <- corrplot(Mcor, order = 'hclust')
dev.off()


selected.wave <- counts[unique(c(test_gene_anno[, 1], test2_gene_anno[, 1], test3_gene_anno[, 1])), group_labels%in%"Ductal"]

Mcor <- cor(as.matrix(t(selected.wave)))
pdf(paste0("scRNAseq_ductal_alldonors_cor_allslice_genes.pdf"), width=12, height = 12)
tes_cor <- corrplot(Mcor, order = 'hclust')
dev.off()

pheno <- my_data@meta.data[["expCond1"]]


selected.wave <- counts[unique(c(test_gene_anno[, 1], test2_gene_anno[, 1], test3_gene_anno[, 1])), group_labels%in%"Beta2"&pheno%in%"LN-ND"]

Mcor <- cor(as.matrix(t(selected.wave)))
pdf(paste0("scRNAseq_LN_beta2_cor_allslice_genes.pdf"), width=12, height = 12)
tes_cor <- corrplot(Mcor, order = 'hclust')
dev.off()


selected.wave <- counts[unique(c(test_gene_anno[, 1], test2_gene_anno[, 1], test3_gene_anno[, 1])), group_labels%in%"Beta2"&pheno%in%"OB-ND"]

Mcor <- cor(as.matrix(t(selected.wave)))
pdf(paste0("scRNAseq_OB_beta2_cor_allslice_genes.pdf"), width=12, height = 12)
tes_cor <- corrplot(Mcor, order = 'hclust')
dev.off()


selected.wave <- counts[unique(c(test_gene_anno[, 1], test2_gene_anno[, 1], test3_gene_anno[, 1])), group_labels%in%"Beta2"&pheno%in%"OW-ND"]

Mcor <- cor(as.matrix(t(selected.wave)))
pdf(paste0("scRNAseq_OW_beta2_cor_allslice_genes.pdf"), width=12, height = 12)
tes_cor <- corrplot(Mcor, order = 'hclust')
dev.off()


selected.wave <- counts[unique(c(test_gene_anno[, 1], test2_gene_anno[, 1], test3_gene_anno[, 1])), group_labels%in%"Beta2"&pheno%in%"T2D"]

Mcor <- cor(as.matrix(t(selected.wave)))
pdf(paste0("scRNAseq_T2D_beta2_cor_allslice_genes.pdf"), width=12, height = 12)
tes_cor <- corrplot(Mcor, order = 'hclust')
dev.off()


############################################
########### Differential stability #########
############################################
list.celltype <- c("Immune", "Quiescent stellate", "Activated stellate", "Endothelial", "Ductal", "Acinar", "Epsilon",
                   "Delta", "Gamma", "Beta6", "Beta5", "Beta4", "Beta3", "Beta2", "Beta1", "Alpha")
list.celltype_br <- c("Immune", "QS", "AS", "Endo", "Ductal", "Acinar", "Epsilon",
                      "Delta", "Gamma", "Beta6", "Beta5", "Beta4", "Beta3", "Beta2", "Beta1", "Alpha")

for(i in c(1:6, 8:16)){

  load(paste0("scRNAseq_", list.celltype_br[i], "_score_donors.rds"))
  rownames(sample.info) <- sample.info[,1]
  kkk <- length(cluster.zscore.summary)
  score.mat <- NULL
  sample.list <- NULL
  for(j in 1:kkk) {
    sample.list <- c(sample.list, as.character(cluster.zscore.summary[[j]]$groupID))
    score.mat <- cbind(score.mat, cluster.zscore.summary[[j]]$Zscores[, "zscore"])
  }
  colnames(score.mat) <- sample.info[sample.list, 2]
  rownames(score.mat) <- rownames(cluster.zscore.summary[[j]]$Zscores)
  save(score.mat, file = paste0("scRNAseq_", list.celltype_br[i], "_score_mat.rds"))

}

my_data <- readRDS("Res29samp_integration_wHippo_alpha_beta_epsilon_wExpCond8_rmOutlier_wAnno.rds")
pheno <- cbind(my_data@meta.data[["orig.ident"]], my_data@meta.data[["expCond1"]])
pheno <- pheno[!duplicated(pheno), ]
pheno[pheno[, 2]%in%"OW-ND", 2] <- "O1"
pheno[pheno[, 2]%in%"OB-ND", 2] <- "O2"
rownames(pheno) <- pheno[, 1]

#####################################
########### T2D vs ND ###############
#####################################
all_genes <- rownames(score.mat)
data("otherRNA")
MT_genes <- otherRNA[otherRNA[, "Family"] == "Mitochondiral", "Gene"]
RP_genes <- otherRNA[otherRNA[, "Family"] == "Ribosomal", "Gene"]
IncRNA_genes <- unique(c(otherRNA[otherRNA[, "Family"] == "IncRNA", "Gene"], all_genes[grep("^LINC", all_genes)], all_genes[grep("^A[A-Z][0-9][0-9][0-9][0-9][0-9][0-9]", all_genes)]))
non_coding <- c(MT_genes, RP_genes, IncRNA_genes)

two_group_test <- function(mat, id1, id2){
  summary_tests <- t(apply(mat, 1, function(x){
    tt.res <- t.test(x[id1], x[id2])
    return(c(tt.res$p.value, tt.res$statistic, tt.res$estimate))
  }))
  summary_tests <- as.data.frame(summary_tests)
  colnames(summary_tests) <- c("pvalue", "t", "mean_x", "mean_y")
  summary_tests$qvalue <- p.adjust(summary_tests$pvalue, method = "fdr")
  summary_tests$diff <- summary_tests$mean_x - summary_tests$mean_y
  return(summary_tests)
}

volcano_plot <- function(summary_tests, label_cutoff, xlab, main_text = ""){
  summary_tests$qvalue_log <- -log10(summary_tests$qvalue)
  sub_id <- which(summary_tests$qvalue <= 0.05)
  subsub_id <- which(summary_tests$qvalue <= 0.05 & abs(summary_tests$diff) >= label_cutoff)
  plot(summary_tests$diff, summary_tests$qvalue_log, pch = 20, xlab = xlab, ylab = "-log10(q-value)", main = main_text)
  points(summary_tests$diff[sub_id], summary_tests$qvalue_log[sub_id], pch = 20, col = "brown")
}

for(i in c(1:6, 8:16)){
  load(paste0("scRNAseq_", list.celltype_br[i], "_score_mat.rds"))
  na_flag <- apply(score.mat, 1, function(x){any(is.na(x))})
  sub.score.mat <- score.mat[!na_flag & !rownames(score.mat)%in%non_coding, ]

  t2d_ids <- which(colnames(sub.score.mat)%in%pheno[pheno[, 2] == "T2D", 1])
  nd_ids <- which(colnames(sub.score.mat)%in%pheno[pheno[, 2] != "T2D", 1])

  test_summary <- two_group_test(sub.score.mat, t2d_ids, nd_ids)
  pdf(paste0("", list.celltype_br[i], "volcano1.pdf"), width=12, height = 12)
  volcano_plot(test_summary, label_cutoff = 2, xlab = "T2D - ND")
  dev.off()

  lean_ids <- which(colnames(sub.score.mat)%in%pheno[pheno[, 2] == "LN-ND", 1])
  test_summary <- two_group_test(sub.score.mat, t2d_ids, lean_ids)

  pdf(paste0("", list.celltype_br[i], "volcano2.pdf"), width=12, height = 12)
  volcano_plot(test_summary, label_cutoff = 2, xlab = "T2D - LN")
  dev.off()
}

pdf(paste0("all_volcano1.pdf"), width=15, height = 12)
par(mfrow = c(3, 5))

for(i in c(1:6, 8:16)){
  load(paste0("scRNAseq_", list.celltype_br[i], "_score_mat.rds"))
  na_flag <- apply(score.mat, 1, function(x){any(is.na(x))})
  sub.score.mat <- score.mat[!na_flag & !rownames(score.mat)%in%non_coding, ]

  t2d_ids <- which(colnames(sub.score.mat)%in%pheno[pheno[, 2] == "T2D", 1])
  nd_ids <- which(colnames(sub.score.mat)%in%pheno[pheno[, 2] != "T2D", 1])

  test_summary <- two_group_test(sub.score.mat, t2d_ids, nd_ids)
  volcano_plot(test_summary, label_cutoff = 2, xlab = "T2D - ND", main_text = list.celltype_br[i])

}
dev.off()

pdf(paste0("selected_volcano1.pdf"), width = 18, height = 4)
par(mfrow = c(1, 6))

for(i in c(15, 14, 13, 3, 5, 1)){
  load(paste0("scRNAseq_", list.celltype_br[i], "_score_mat.rds"))
  na_flag <- apply(score.mat, 1, function(x){any(is.na(x))})
  sub.score.mat <- score.mat[!na_flag & !rownames(score.mat)%in%non_coding, ]

  t2d_ids <- which(colnames(sub.score.mat)%in%pheno[pheno[, 2] == "T2D", 1])
  nd_ids <- which(colnames(sub.score.mat)%in%pheno[pheno[, 2] != "T2D", 1])

  test_summary <- two_group_test(sub.score.mat, t2d_ids, nd_ids)
  volcano_plot(test_summary, label_cutoff = 2, xlab = "T2D - ND", main_text = list.celltype_br[i])
}
dev.off()



pdf(paste0("all_volcano2.pdf"), width=15, height = 12)
par(mfrow = c(3, 5))

for(i in c(1:6, 8:16)){
  load(paste0("scRNAseq_", list.celltype_br[i], "_score_mat.rds"))
  na_flag <- apply(score.mat, 1, function(x){any(is.na(x))})
  sub.score.mat <- score.mat[!na_flag & !rownames(score.mat)%in%non_coding, ]

  t2d_ids <- which(colnames(sub.score.mat)%in%pheno[pheno[, 2] == "T2D", 1])
  nd_ids <- which(colnames(sub.score.mat)%in%pheno[pheno[, 2] != "T2D", 1])

  lean_ids <- which(colnames(sub.score.mat)%in%pheno[pheno[, 2] == "LN-ND", 1])
  test_summary <- two_group_test(sub.score.mat, t2d_ids, lean_ids)

  volcano_plot(test_summary, label_cutoff = 2, xlab = "T2D - LN", main_text = list.celltype_br[i])

}
dev.off()


########################################################
########### Differential stability genes ###############
########################################################
i <- 14
load(paste0("scRNAseq_", list.celltype_br[i], "_score_mat.rds"))
na_flag <- apply(score.mat, 1, function(x){any(is.na(x))})
sub.score.mat <- score.mat[!na_flag & !rownames(score.mat)%in%non_coding, ]

t2d_ids <- which(colnames(sub.score.mat)%in%pheno[pheno[, 2] == "T2D", 1])
nd_ids <- which(colnames(sub.score.mat)%in%pheno[pheno[, 2] != "T2D", 1])

test_summary <- two_group_test(sub.score.mat, t2d_ids, nd_ids)
hits <- test_summary[test_summary$qvalue <= 0.05, ]
hits <- hits[order(abs(hits$diff), decreasing = TRUE), ]

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
# Example list of gene names
geneList <- rownames(hits)[hits$diff <= -5]
writeLines(geneList, "down_gene.list")
# Perform GO enrichment analysis using default settings
enrich_result <- enrichGO(gene = geneList, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")

pdf(paste0("DS_enrichment.pdf"), width=5, height = 10)
dotplot(enrich_result, showCategory=30) + ggtitle("dotplot for BP")
dev.off()



################################
########## Figure 4 heatmap ####
################################

genes <- unique(c(test_gene_anno[, 1], test2_gene_anno[, 1], test3_gene_anno[, 1]))

subcounts <- counts[genes, donor.info%in%"HI_25" & group_labels%in%"Beta2"]
library(ComplexHeatmap)
library(viridis)
Colors=viridis(40)


centered <- apply(subcounts, 1, function(x){
  x
})

kkk <- t(centered)
rownames(kkk) <- rownames(subcounts)
pdf(paste0("scRNAseq_HI_25_Beta2_heatmap_genes.pdf"), width=8, height = 10)
Heatmap(as.matrix(kkk), name = "beta HI 25", col=Colors, column_labels = rep("", ncol(subcounts)))
dev.off()

kkk[kkk > 20] <- 20
pdf(paste0("scRNAseq_HI_25_Beta2_heatmap_genes_rm_INS.pdf"), width=8, height = 10)
Heatmap(kkk[!rownames(kkk)%in%c("INS", "IAPP", "GCG", "SST"), ], name = "beta HI 25", col=Colors, column_labels = rep("", ncol(subcounts)))
dev.off()


list1 <- c("NRG3", "DNAJC6", "NEGR1", "SEMA5A", "ROBO2", "ROBO1", "SYT14", "NRG1", "ROR1")
list2 <- c("HLA-B", "HLA-A", "HLA-C", "PCDH17", "HSP90B1", "HSPA5", "SCG3", "PSAP", "CALM1", "PCSK1N", "SCG5", "HSPA8", "HSP90AA1", "HSP90AB1", "VGF", "CHGA", "CHGB", "SCG2")
list3 <- c("JUN", "HSPA6", "DNAJB1", "HSPA1A", "HSPA1B", "HSPB1", "HSPH1")
list4 <- c( "IL32", "REG3A", "REG1A", "REG1B")
combined <- c(list1, list2, list3, list4)
kkk[kkk > 20] <- 20


gene_anno_modules <- data.frame(Gene = c(list1, list2, list3, list4), Family = c(rep("Module 1", length(list1)), rep("Module 4", length(list2)), rep("Module 3", length(list3)), rep("Module 2", length(list4))))
rownames(gene_anno_modules) <- gene_anno_modules$Gene
pdf(paste0("scRNAseq_HI_25_Beta2_heatmap_genes_rm_INS.pdf"), width=8, height = 10)
aaa <- Heatmap(kkk[gene_anno_modules$Gene, ], row_split = factor(gene_anno_modules[, "Family"]), name = "beta HI 25", col=Colors, column_labels = rep("", ncol(subcounts)), cluster_rows = FALSE)
aaa
dev.off()

cell_reorder <- aaa@column_order


centered <- apply(subcounts, 1, function(x){
  x <- log(x+0.1)
  x-mean(x)
})

kkk <- t(centered)
rownames(kkk) <- rownames(subcounts)
Colors=magma(40)
pdf(paste0("scRNAseq_HI_25_Beta2_heatmap_genes_INS.pdf"), width=8, height = 10)
Heatmap(kkk[c("INS", "IAPP", "GCG", "SST"), cell_reorder], name = "beta HI 25", col=Colors, column_labels = rep("", ncol(subcounts)), cluster_columns = FALSE, cluster_rows = FALSE)
dev.off()

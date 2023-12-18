
waveplot_fig3_override <- function(one.cluster.zscore.summary, show_mean_cutoff = 1,
                                   log_xmax = 2, dot_color, hitnames){


  dispersion <- one.cluster.zscore.summary$NBdispersion
  selected_stats <- one.cluster.zscore.summary$Inflated
  xmax <- max(round(selected_stats$gene_mean))+3

  showhits <- selected_stats[hitnames, ]

  score_mat <- one.cluster.zscore.summary$Zscores
  all_genes <- rownames(score_mat)

  selected.IDs <- which(all_genes %in% rownames(selected_stats))
  background_stats <- score_mat[-selected.IDs, ]

  stats_combined <- data.frame(category = c(rep("B", nrow(background_stats)), rep("H", nrow(selected_stats))),
                               rbind(background_stats, selected_stats))

  x_intervals <- seq(0.1, xmax, by = 0.2)
  p_intervals <- rep(0, length(x_intervals))
  logx_intervals <- log10(x_intervals + 0.1)


  upper_lower_intervals <- matrix(0, ncol = 1, nrow = length(x_intervals))
  for(i in 1:length(x_intervals)){
    p_intervals[i] <- determine_k(x_intervals[i])
    upper_lower_intervals[i, 1] <- compute_expected_pi_for_k_NB(x_intervals[i], p_intervals[i], dispersion)
  }
  expected_line <- data.frame(aa = x_intervals, bb = logx_intervals, dd = upper_lower_intervals[, 1])

  stats_combined[, "gene_mean"] <- log10(stats_combined[, "gene_mean"]+0.1)
  showhits[, "gene_mean"] <- log10(showhits[, "gene_mean"]+0.1)
  xmax <- log_xmax #max(stats_combined[, "gene_mean"]) + 0.2
  require(ggplot2)
  g_basic <- ggplot2::ggplot(stats_combined, ggplot2::aes(x = gene_mean,
                                                          y = k_proportion, col = category)) +
    ggplot2::geom_point(size = 2, alpha = 0.5) +
    ggplot2::geom_text(data = showhits, aes(label = rownames(showhits)), vjust = -1, size = 3, color = "brown") +
    ggplot2::geom_point(data = showhits, size = 3, alpha = 0.5, color = "brown") +
    ggplot2::theme(legend.position = "topright") + ggplot2::theme_bw() +
    ggplot2::ylab("K Proportion") + ggplot2::xlab("Gene Mean") + xlim(-1, xmax) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.grid = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   strip.placement = "inside") + ggplot2::geom_line(data = expected_line, aes(bb, dd), color = "#023373", show.legend = FALSE)

  g <- g_basic +  ggplot2::scale_color_manual(values = c("B" = dot_color, "H" ="#4E64A6"),
                                              labels = c("B" = "Homeostatic", "H" = "Inflated"))


  return(g)

}

##########################################
########### spliced vs unspliced #########
##########################################

spliced_counts <- read.csv("spliced.csv")
unspliced_counts <- read.csv("unspliced.csv")
adata_var <- read.csv("adata_var.csv")
adata_obs <- read.csv("adata_obs.csv")

group_labels <- adata_obs$celltype_final
donor.info <- adata_obs$orig.ident
list.celltype <- c("Immune", "Quiescent stellate", "Activated stellate", "Endothelial", "Ductal", "Acinar", "Epsilon",
                   "Delta", "Gamma", "Beta6", "Beta5", "Beta4", "Beta3", "Beta2", "Beta1", "Alpha")
list.celltype_br <- c("Immune", "QS", "AS", "Endo", "Ductal", "Acinar", "Epsilon",
                      "Delta", "Gamma", "Beta6", "Beta5", "Beta4", "Beta3", "Beta2", "Beta1", "Alpha")



load("scRNAseq_Beta2_HI_25_showhits.rds")
list.celltype_br <- c("Immune", "QS", "AS", "Endo", "Ductal", "Acinar", "Epsilon",
                      "Delta", "Gamma", "Beta6", "Beta5", "Beta4", "Beta3", "Beta2", "Beta1", "Alpha")
i <- 14
j <- 17

selected <- which(group_labels %in% list.celltype[i])
testdata <- t(spliced_counts[selected, ])
rownames(testdata) <- adata_var[, 1]
sub.donor.info <- donor.info[group_labels %in% list.celltype[i]]
new.clusters <- as.numeric(as.factor(sub.donor.info))

sample.info <- cbind(new.clusters, sub.donor.info)
sample.info <- sample.info[!duplicated(sample.info), ]
rownames(sample.info) <- sample.info[, 1]

cluster.zscore.summary <- NB_inflation_test_var_asymp(testdata, grouping = new.clusters, min.depth.group = 1000, cell.depth.ranges = c(500, 60000), genemean.cutoff = 0)

wave_one <- waveplot_fig3_override(cluster.zscore.summary[[j]], show_mean_cutoff = 0,
                                   dot_color = "#5C4B51", hitnames = genes)

ggsave(paste0("scRNAseq_unspliced_", list.celltype_br[i], "_",  sample.info[as.character(cluster.zscore.summary[[j]]$groupID), 2], "_wave.pdf"), plot = wave_one, width = 12, height = 3)
save(cluster.zscore.summary, sample.info, file = paste0("scRNAseq_spliced_", list.celltype_br[i], "_all_genes.rds"))


i <- 14
j <- 17

selected <- which(group_labels %in% list.celltype[i])
testdata <- t(unspliced_counts[selected, ])
rownames(testdata) <- adata_var[, 1]
sub.donor.info <- donor.info[group_labels %in% list.celltype[i]]
new.clusters <- as.numeric(as.factor(sub.donor.info))

sample.info <- cbind(new.clusters, sub.donor.info)
sample.info <- sample.info[!duplicated(sample.info), ]
rownames(sample.info) <- sample.info[, 1]

cluster.zscore.summary <- NB_inflation_test_var_asymp(testdata, grouping = new.clusters, min.depth.group = 1000, cell.depth.ranges = c(500, 60000), genemean.cutoff = 0)


wave_one <- waveplot_fig3_override(cluster.zscore.summary[[j]], show_mean_cutoff = 1,
                                   dot_color = "#5C4B51", hitnames = genes)

ggsave(paste0("scRNAseq_spliced_", list.celltype_br[i], "_",  sample.info[as.character(cluster.zscore.summary[[j]]$groupID), 2], "_wave.pdf"), plot = wave_one, width = 12, height = 3)
save(cluster.zscore.summary, sample.info, file = paste0("T2D_spliced_count/scRNAseq_unspliced_", list.celltype_br[i], "_all_genes.rds"))








spliced_unspliced_figure3_pairwise_scatter <- function(summary1, summary2, showgenes){

  score_mat1 <- summary1$Zscores
  score_mat2 <- summary2$Zscores
  colnames(score_mat1) <- paste0("unspliced_", colnames(score_mat1))
  colnames(score_mat2) <- paste0("spliced_", colnames(score_mat2))
  score_mat <- cbind(score_mat1, score_mat2)
  all_genes <- rownames(score_mat)

  data("otherRNA")
  MT_genes <- otherRNA[otherRNA[, "Family"] == "Mitochondiral", "Gene"]
  RP_genes <- otherRNA[otherRNA[, "Family"] == "Ribosomal", "Gene"]
  IncRNA_genes <- otherRNA[otherRNA[, "Family"] == "IncRNA", "Gene"]

  MT.IDs <- which(all_genes%in%MT_genes)
  RP.IDs <- which(all_genes%in%RP_genes)
  IncRNA.IDs <- c(which(all_genes%in%IncRNA_genes), grep("^LINC", all_genes), grep("^A[A-Z][0-9][0-9][0-9][0-9][0-9][0-9]", all_genes))

  background_stats <- score_mat[-c(MT.IDs, RP.IDs, IncRNA.IDs), ]
  showhits <- background_stats[showgenes, ]

  require(ggplot2)
  g1 <- ggplot2::ggplot(background_stats, ggplot2::aes(x = unspliced_gene_mean,
                                                       y = spliced_gene_mean)) +
    ggplot2::geom_point(size = 2, alpha = 0.5) +
    ggplot2::geom_point(data = showhits, size = 3, alpha = 0.5, color = "brown") +
    ggplot2::geom_text(data = showhits, aes(label = rownames(showhits)), vjust = -1, size = 3, color = "brown")  + ggplot2::theme_bw() +
    ggplot2::ylab("Spliced Mean UMI") + ggplot2::xlab("Unspliced Mean UMI") + xlim(0, 20)  + ylim(0, 20) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.grid = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   strip.placement = "inside") + geom_abline(intercept = 0, slope = 1, color="darkblue",
                                                             linetype="dashed", linewidth=1.5)
  g2 <- ggplot2::ggplot(background_stats, ggplot2::aes(x = unspliced_zscore,
                                                       y = spliced_zscore)) +
    ggplot2::geom_point(size = 2, alpha = 0.5) +
    ggplot2::geom_point(data = showhits, size = 3, alpha = 0.5, color = "brown") +
    ggplot2::geom_text(data = showhits, aes(label = rownames(showhits)), vjust = -1, size = 3, color = "brown")  + ggplot2::theme_bw() +
    ggplot2::ylab("Spliced Z-index") + ggplot2::xlab("Unspliced Z-index") +  xlim(-0.5, 40)  + ylim(-0.5, 40) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.grid = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   strip.placement = "inside") + geom_abline(intercept = 0, slope = 1, color="darkblue",
                                                             linetype="dashed", linewidth=1.5)

  allg <- list(g1, g2)
  return(allg)

}

load("scRNAseq_unspliced_Beta2_all_genes.rds")
unspliced.summary <- cluster.zscore.summary
load("scRNAseq_spliced_Beta2_all_genes.rds")
spliced.summary <- cluster.zscore.summary

load("scRNAseq_Beta2_HI_25_showhits.rds")

j <- 17
test_one <- spliced_unspliced_figure3_pairwise_scatter(unspliced.summary[[j]], spliced.summary[[j]], showgenes = genes)
pdf(paste0("scRNAseq_spliced_unspliced_comp_", list.celltype_br[i], "_",  sample.info[as.character(cluster.zscore.summary[[j]]$groupID), 2], "_n.pdf"), width = 12, height = 6)
ggplot2.multiplot(test_one[[1]], test_one[[2]], cols=2)
dev.off()


other_stats_plot_figure2_override <- function(one.cluster.zscore.summary, show_mean_cutoff = 1.5, show_top = 30, log_scale = TRUE, log_xmax = 2){

  dispersion <- one.cluster.zscore.summary$NBdispersion
  selected_stats <- one.cluster.zscore.summary$Inflated

  xmax <-  max(round(selected_stats$gene_mean))+3

  score_mat <- one.cluster.zscore.summary$Zscores
  all_genes <- rownames(score_mat)

  selected.IDs <- which(all_genes %in% rownames(selected_stats))

  add_mat <- one.cluster.zscore.summary$Add
  selected_stats <- cbind(selected_stats, add_mat[rownames(selected_stats), ])
  subselected_stats <- selected_stats[selected_stats$gene_mean >= show_mean_cutoff, ]
  hits_num <- nrow(subselected_stats)
  show_top <- min(hits_num, show_top)
  showhits <- subselected_stats[1:show_top, ]

  data("otherRNA")
  MT_genes <- otherRNA[otherRNA[, "Family"] == "Mitochondiral", "Gene"]
  RP_genes <- otherRNA[otherRNA[, "Family"] == "Ribosomal", "Gene"]
  IncRNA_genes <- otherRNA[otherRNA[, "Family"] == "IncRNA", "Gene"]

  MT.IDs <- which(all_genes%in%MT_genes)
  RP.IDs <- which(all_genes%in%RP_genes)
  IncRNA.IDs <- c(which(all_genes%in%IncRNA_genes), grep("^LINC", all_genes), grep("^A[A-Z][0-9][0-9][0-9][0-9][0-9][0-9]", all_genes))


  background_stats <- cbind(score_mat[-c(selected.IDs, MT.IDs, RP.IDs, IncRNA.IDs), ], add_mat[-c(selected.IDs, MT.IDs, RP.IDs, IncRNA.IDs), ])
  stats_combined <- data.frame(category = c(rep("B", nrow(background_stats)), rep("H", nrow(selected_stats))),
                               rbind(background_stats, selected_stats))

  stats_combined$lognorm_cv2 <- stats_combined$lognorm_cv^2
  showhits$lognorm_cv2 <- showhits$lognorm_cv^2

  x_intervals <- seq(0.1, 100, by = 0.2)
  p_intervals <- rep(0, length(x_intervals))
  logx_intervals <- log10(x_intervals + 0.1)


  upper_lower_intervals <- matrix(0, ncol = 1, nrow = length(x_intervals))
  for(i in 1:length(x_intervals)){
    p_intervals[i] <- determine_k(x_intervals[i])
    upper_lower_intervals[i, 1] <- compute_expected_pi_for_k_NB(x_intervals[i], p_intervals[i], dispersion)
  }
  expected_line <- data.frame(aa = x_intervals, bb = logx_intervals, dd = upper_lower_intervals[, 1])


  if(log_scale == FALSE){
    require(ggplot2)
    g1 <- ggplot2::ggplot(stats_combined, ggplot2::aes(x = count_mean,
                                                       y = count_variance, col = category)) +
      ggplot2::geom_point(size = 2, alpha = 0.5) +
      ggplot2::geom_text(data = showhits, aes(label = rownames(showhits)), vjust = -1, size = 3, color = "brown") +
      ggplot2::theme(legend.position = "topright") + ggplot2::theme_bw() +
      ggplot2::ylab("Variance") + ggplot2::xlab("Gene Mean") + xlim(0, xmax) + ylim(0, 100)+
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     panel.grid = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     strip.placement = "inside")
    g2 <- ggplot2::ggplot(stats_combined, ggplot2::aes(x = count_mean,
                                                       y = count_cv, col = category)) +
      ggplot2::geom_point(size = 2, alpha = 0.5) +
      ggplot2::geom_text(data = showhits, aes(label = rownames(showhits)), vjust = -1, size = 3, color = "brown") +
      ggplot2::theme(legend.position = "topright") + ggplot2::theme_bw() +
      ggplot2::ylab("CV") + ggplot2::xlab("Gene Mean") + xlim(0, xmax) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     panel.grid = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     strip.placement = "inside")

    g3 <- ggplot2::ggplot(stats_combined, ggplot2::aes(x = lognorm_mean,
                                                       y = lognorm_cv2, col = category)) +
      ggplot2::geom_point(size = 2, alpha = 0.5) +
      ggplot2::geom_text(data = showhits, aes(label = rownames(showhits)), vjust = -1, size = 3, color = "brown") +
      ggplot2::theme(legend.position = "topright") + ggplot2::theme_bw() +
      ggplot2::ylab("lognormal CV2") + ggplot2::xlab("Gene Mean") + xlim(0, xmax) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     panel.grid = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     strip.placement = "inside")

    g_basic <- ggplot2::ggplot(stats_combined, ggplot2::aes(x = gene_mean,
                                                            y = k_proportion, col = category)) +
      ggplot2::geom_point(size = 2, alpha = 0.5) +
      ggplot2::geom_text(data = showhits, aes(label = rownames(showhits)), vjust = -1, size = 3, color = "brown") +
      ggplot2::theme(legend.position = "topright") + ggplot2::theme_bw() +
      ggplot2::ylab("K Proportion") + ggplot2::xlab("Gene Mean") + xlim(0, xmax) + ylim(0, 1)+
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     panel.grid = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     strip.placement = "inside") + ggplot2::geom_line(data = expected_line, aes(aa, dd), color = "#023373", show.legend = FALSE)


    g_basic <- g_basic +  ggplot2::scale_color_manual(values = c("B" = "black", "H" ="#F06060"),
                                                      labels = c("B" = "Background", "H" = "Inflated"))
    g1 <- g1 +  ggplot2::scale_color_manual(values = c("B" = "black", "H" ="#F06060"),
                                            labels = c("B" = "Background", "H" = "Inflated"))
    g2 <- g2 +  ggplot2::scale_color_manual(values = c("B" = "black", "H" ="#F06060"),
                                            labels = c("B" = "Background", "H" = "Inflated"))
    g3 <- g3 +  ggplot2::scale_color_manual(values = c("B" = "black", "H" ="#F06060"),
                                            labels = c("B" = "Background", "H" = "Inflated"))

    allg <- list(g_basic, g1, g2, g3)

  } else {

    stats_combined[, "gene_mean"] <- log10(stats_combined[, "gene_mean"]+0.1)
    showhits[, "gene_mean"] <- log10(showhits[, "gene_mean"]+0.1)
    xmax <- log_xmax #max(stats_combined[, "gene_mean"]) + 0.2

    require(ggplot2)
    g_basic <- ggplot2::ggplot(stats_combined, ggplot2::aes(x = gene_mean,
                                                            y = k_proportion, col = category)) +
      ggplot2::geom_point(size = 2, alpha = 0.5) +
      ggplot2::geom_point(data = showhits, size = 3, alpha = 0.5, color = "brown") +
      ggplot2::geom_text(data = showhits, aes(label = rownames(showhits)), vjust = -1, size = 3, color = "brown") +
      ggplot2::theme(legend.position = "topright") + ggplot2::theme_bw() +
      ggplot2::ylab("K Proportion") + ggplot2::xlab("log10 Gene Mean") + xlim(-1, xmax) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     panel.grid = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     strip.placement = "inside") + ggplot2::geom_line(data = expected_line, aes(bb, dd), color = "#023373", show.legend = FALSE)

    g1 <- ggplot2::ggplot(stats_combined, ggplot2::aes(x = gene_mean,
                                                       y = lognorm_variance, col = category)) +
      ggplot2::geom_point(size = 2, alpha = 0.5) +
      ggplot2::geom_point(data = showhits, size = 3, alpha = 0.5, color = "brown") +
      ggplot2::geom_text(data = showhits, aes(label = rownames(showhits)), vjust = -1, size = 3, color = "brown") +
      ggplot2::theme(legend.position = "topright") + ggplot2::theme_bw() +
      ggplot2::ylab("Variance") + ggplot2::xlab("log10 Gene Mean") + xlim(-1, xmax) + ylim(0, max(stats_combined$lognorm_variance))+
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     panel.grid = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     strip.placement = "inside")
    g2 <- ggplot2::ggplot(stats_combined, ggplot2::aes(x = gene_mean,
                                                       y = lognorm_cv, col = category)) +
      ggplot2::geom_point(size = 2, alpha = 0.5) +
      ggplot2::geom_point(data = showhits, size = 3, alpha = 0.5, color = "brown") +
      ggplot2::geom_text(data = showhits, aes(label = rownames(showhits)), vjust = -1, size = 3, color = "brown") +
      ggplot2::theme(legend.position = "topright") + ggplot2::theme_bw() +
      ggplot2::ylab("CV") + ggplot2::xlab("log10 Gene Mean") + xlim(-1, xmax) + ylim(0, 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     panel.grid = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     strip.placement = "inside")

    g3 <- ggplot2::ggplot(stats_combined, ggplot2::aes(x = gene_mean,
                                                       y = lognorm_cv2, col = category)) +
      ggplot2::geom_point(size = 2, alpha = 0.5) +
      ggplot2::geom_point(data = showhits, size = 3, alpha = 0.5, color = "brown") +
      ggplot2::geom_text(data = showhits, aes(label = rownames(showhits)), vjust = -1, size = 3, color = "brown") +
      ggplot2::theme(legend.position = "topright") + ggplot2::theme_bw() +
      ggplot2::ylab("lognormal CV2") + ggplot2::xlab("log10 Gene Mean") + xlim(-1, xmax) + ylim(0, 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     panel.grid = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     strip.placement = "inside")


    g_basic <- g_basic +  ggplot2::scale_color_manual(values = c("B" = "#5C4B51", "H" ="#4E64A6"),
                                                      labels = c("B" = "Background", "H" = "Inflated"))
    g1 <- g1 +  ggplot2::scale_color_manual(values = c("B" = "#5C4B51", "H" ="#4E64A6"),
                                            labels = c("B" = "Background", "H" = "Inflated"))
    g2 <- g2 +  ggplot2::scale_color_manual(values = c("B" = "#5C4B51", "H" ="#4E64A6"),
                                            labels = c("B" = "Background", "H" = "Inflated"))
    g3 <- g3 +  ggplot2::scale_color_manual(values = c("B" = "#5C4B51", "H" ="#4E64A6"),
                                            labels = c("B" = "Background", "H" = "Inflated"))

  }

  allg <- list(g_basic, g1, g2, g3)
  return(allg)

}





other_stats_plot_figure2_pairwise_scatter <- function(one.cluster.zscore.summary, show_mean_cutoff = 0, show_top = 30, log_scale = TRUE){

  dispersion <- one.cluster.zscore.summary$NBdispersion
  selected_stats <- one.cluster.zscore.summary$Inflated

  xmax <- max(round(selected_stats$gene_mean))+3



  score_mat <- one.cluster.zscore.summary$Zscores
  all_genes <- rownames(score_mat)

  selected.IDs <- which(all_genes %in% rownames(selected_stats))

  add_mat <- one.cluster.zscore.summary$Add
  selected_stats <- cbind(selected_stats, add_mat[rownames(selected_stats), ])
  subselected_stats <- selected_stats[selected_stats$gene_mean >= show_mean_cutoff, ]
  hits_num <- nrow(subselected_stats)
  show_top <- min(hits_num, show_top)
  showhits <- subselected_stats[1:show_top, ]

  data("otherRNA")
  MT_genes <- otherRNA[otherRNA[, "Family"] == "Mitochondiral", "Gene"]
  RP_genes <- otherRNA[otherRNA[, "Family"] == "Ribosomal", "Gene"]
  IncRNA_genes <- otherRNA[otherRNA[, "Family"] == "IncRNA", "Gene"]

  MT.IDs <- which(all_genes%in%MT_genes)
  RP.IDs <- which(all_genes%in%RP_genes)
  IncRNA.IDs <- c(which(all_genes%in%IncRNA_genes), grep("^LINC", all_genes), grep("^A[A-Z][0-9][0-9][0-9][0-9][0-9][0-9]", all_genes))


  background_stats <- cbind(score_mat[-c(selected.IDs, MT.IDs, RP.IDs, IncRNA.IDs), ], add_mat[-c(selected.IDs, MT.IDs, RP.IDs, IncRNA.IDs), ])

  stats_combined <- data.frame(category = c(rep("B", nrow(background_stats)), rep("H", nrow(selected_stats))),
                               rbind(background_stats, selected_stats))

  stats_combined$lognorm_cv2 <- stats_combined$lognorm_cv^2
  showhits$lognorm_cv2 <- showhits$lognorm_cv^2


  if(log_scale == TRUE) {

    stats_combined[, "gene_mean"] <- log10(stats_combined[, "gene_mean"]+0.1)
    showhits[, "gene_mean"] <- log10(showhits[, "gene_mean"]+0.1)
    xmax <- max(showhits[, "gene_mean"]) + 0.2

    require(ggplot2)
    g1 <- ggplot2::ggplot(stats_combined, ggplot2::aes(x = lognorm_variance,
                                                       y = k_proportion, col = category)) +
      ggplot2::geom_point(size = 2, alpha = 0.5) +
      ggplot2::geom_point(data = showhits, size = 3, alpha = 0.5, color = "brown") +
      ggplot2::geom_text(data = showhits, aes(label = rownames(showhits)), vjust = -1, size = 3, color = "brown")  + ggplot2::theme_bw() +
      ggplot2::ylab("K Proportion") + ggplot2::xlab("lognormal variance") + xlim(0, 2.2) + ylim(0.18, 1) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     panel.grid = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     strip.placement = "inside")

    g2 <- ggplot2::ggplot(stats_combined, ggplot2::aes(x = lognorm_cv2,
                                                       y = k_proportion, col = category)) +
      ggplot2::geom_point(size = 2, alpha = 0.5) +
      ggplot2::geom_point(data = showhits, size = 3, alpha = 0.5, color = "brown") +
      ggplot2::geom_text(data = showhits, aes(label = rownames(showhits)), vjust = -1, size = 3, color = "brown") + ggplot2::theme_bw() +
      ggplot2::ylab("K Proportion") + ggplot2::xlab("lognormal CV2") + xlim(0, 10) + ylim(0.18, 1)+
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     panel.grid = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     strip.placement = "inside")

    g3 <- ggplot2::ggplot(stats_combined, ggplot2::aes(x = lognorm_variance,
                                                       y = zscore, col = category)) +
      ggplot2::geom_point(size = 2, alpha = 0.5) +
      ggplot2::geom_point(data = showhits, size = 3, alpha = 0.5, color = "brown") +
      ggplot2::geom_text(data = showhits, aes(label = rownames(showhits)), vjust = -1, size = 3, color = "brown")  + ggplot2::theme_bw() +
      ggplot2::ylab("Z Index") + ggplot2::xlab("lognormal variance") + xlim(0, 2.2) + ylim(0, 30) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     panel.grid = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     strip.placement = "inside")

    g4 <- ggplot2::ggplot(stats_combined, ggplot2::aes(x = lognorm_cv2,
                                                       y = zscore, col = category)) +
      ggplot2::geom_point(size = 2, alpha = 0.5) +
      ggplot2::geom_point(data = showhits, size = 3, alpha = 0.5, color = "brown") +
      ggplot2::geom_text(data = showhits, aes(label = rownames(showhits)), vjust = -1, size = 3, color = "brown") + ggplot2::theme_bw() +
      ggplot2::ylab("Z Index") + ggplot2::xlab("lognormal CV2") + xlim(0, 10) + ylim(0, 30)  +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     panel.grid = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     strip.placement = "inside")

    g4 <- g4 +  ggplot2::scale_color_manual(values = c("B" = "#5C4B51", "H" ="#4E64A6"),
                                            labels = c("B" = "Background", "H" = "Inflated"))  + ggplot2::theme(legend.position = "none")
    g1 <- g1 +  ggplot2::scale_color_manual(values = c("B" = "#5C4B51", "H" ="#4E64A6"),
                                            labels = c("B" = "Background", "H" = "Inflated"))  + ggplot2::theme(legend.position = "none")
    g2 <- g2 +  ggplot2::scale_color_manual(values = c("B" = "#5C4B51", "H" ="#4E64A6"),
                                            labels = c("B" = "Background", "H" = "Inflated"))  + ggplot2::theme(legend.position = "none")
    g3 <- g3 +  ggplot2::scale_color_manual(values = c("B" = "#5C4B51", "H" ="#4E64A6"),
                                            labels = c("B" = "Background", "H" = "Inflated")) + ggplot2::theme(legend.position = "none")

  }

  allg <- list(g1, g2, g3, g4)
  return(allg)

}


#######################################
########## other Stats comparison #####
#######################################

library(RegulationZIndex)
library(Seurat)

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

i<-14
j<-16

selected <- which(group_labels %in% list.celltype[i])
testdata <- counts[, selected]

sub.donor.info <- donor.info[group_labels %in% list.celltype[i]]
new.clusters <- as.numeric(as.factor(sub.donor.info))

sample.info <- cbind(new.clusters, sub.donor.info)
sample.info <- sample.info[!duplicated(sample.info), ]
data("otherRNA")
cluster.zscore.summary <- NB_inflation_test_var_asymp_collecting_more_metrics(testdata, grouping = new.clusters, genemean.cutoff = 0, filter.gene = TRUE)


rownames(sample.info) <- sample.info[, 1]
kkk <- length(cluster.zscore.summary)
colnames(cluster.zscore.summary[[j]]$Add) <- c("count_variance", "count_mean", "count_cv", "lognorm_variance", "lognorm_mean", "lognorm_cv")
wave_one <- other_stats_plot_figure2_override(cluster.zscore.summary[[j]], show_mean_cutoff = 1, log_scale = TRUE, show_top = 15)

library(easyGgplot2)
pdf(paste0("scRNAseq_", list.celltype_br[i], "_",  sample.info[as.character(cluster.zscore.summary[[j]]$groupID), 2], "_more_stats_fig2.pdf"), width = 12, height = 9)
ggplot2.multiplot(wave_one[[1]], wave_one[[2]], wave_one[[4]], cols=1)
dev.off()


scatter_all <- other_stats_plot_figure2_pairwise_scatter(cluster.zscore.summary[[j]], show_mean_cutoff = 1, log_scale = TRUE, show_top = 15)
pdf(paste0("scRNAseq_", list.celltype_br[i], "_",  sample.info[as.character(cluster.zscore.summary[[j]]$groupID), 2], "_more_stats_pairwise_fig2.pdf"), width = 12, height = 12)
ggplot2.multiplot(scatter_all[[1]], scatter_all[[2]], scatter_all[[3]], scatter_all[[4]], cols=2)
dev.off()


all_genes <- rownames(cluster.zscore.summary[[j]]$Zscores)

data("otherRNA")
MT_genes <- otherRNA[otherRNA[, "Family"] == "Mitochondiral", "Gene"]
RP_genes <- otherRNA[otherRNA[, "Family"] == "Ribosomal", "Gene"]
IncRNA_genes <- otherRNA[otherRNA[, "Family"] == "IncRNA", "Gene"]

MT.IDs <- which(all_genes%in%MT_genes)
RP.IDs <- which(all_genes%in%RP_genes)
IncRNA.IDs <- c(which(all_genes%in%IncRNA_genes), grep("^LINC", all_genes), grep("^A[A-Z][0-9][0-9][0-9][0-9][0-9][0-9]", all_genes))

coding_genes <- all_genes[-c(IncRNA.IDs, RP.IDs, MT.IDs)]

www_selected <- cbind(cluster.zscore.summary[[j]]$Zscores[coding_genes, ], cluster.zscore.summary[[j]]$Add[coding_genes, ])
www_selected_1 <- www_selected[order(www_selected$lognorm_variance, decreasing = T), ]

high_var_genes <- rownames(www_selected_1[1:20, ])



other_stats_plot_supp_figure2_override <- function(one.cluster.zscore.summary, show_mean_cutoff = 1.5, show_top = 30, log_scale = TRUE, log_xmax = 2, showgenes){

  dispersion <- one.cluster.zscore.summary$NBdispersion
  selected_stats <- one.cluster.zscore.summary$Inflated

  xmax <-  max(round(selected_stats$gene_mean))+3

  score_mat <- one.cluster.zscore.summary$Zscores
  all_genes <- rownames(score_mat)

  selected.IDs <- which(all_genes %in% rownames(selected_stats))

  add_mat <- one.cluster.zscore.summary$Add
  selected_stats <- cbind(selected_stats, add_mat[rownames(selected_stats), ])
  subselected_stats <- selected_stats[selected_stats$gene_mean >= show_mean_cutoff, ]
  hits_num <- nrow(subselected_stats)
  show_top <- min(hits_num, show_top)
  showhits <- subselected_stats[1:show_top, ]

  data("otherRNA")
  MT_genes <- otherRNA[otherRNA[, "Family"] == "Mitochondiral", "Gene"]
  RP_genes <- otherRNA[otherRNA[, "Family"] == "Ribosomal", "Gene"]
  IncRNA_genes <- otherRNA[otherRNA[, "Family"] == "IncRNA", "Gene"]

  MT.IDs <- which(all_genes%in%MT_genes)
  RP.IDs <- which(all_genes%in%RP_genes)
  IncRNA.IDs <- c(which(all_genes%in%IncRNA_genes), grep("^LINC", all_genes), grep("^A[A-Z][0-9][0-9][0-9][0-9][0-9][0-9]", all_genes))


  background_stats <- cbind(score_mat[-c(selected.IDs, MT.IDs, RP.IDs, IncRNA.IDs), ], add_mat[-c(selected.IDs, MT.IDs, RP.IDs, IncRNA.IDs), ])
  stats_combined <- data.frame(category = c(rep("B", nrow(background_stats)), rep("H", nrow(selected_stats))),
                               rbind(background_stats, selected_stats))



  stats_combined$lognorm_cv2 <- stats_combined$lognorm_cv^2
  showhits$lognorm_cv2 <- showhits$lognorm_cv^2

  showmore <- stats_combined[showgenes, ]

  x_intervals <- seq(0.1, 100, by = 0.2)
  p_intervals <- rep(0, length(x_intervals))
  logx_intervals <- log10(x_intervals + 0.1)

  upper_lower_intervals <- matrix(0, ncol = 1, nrow = length(x_intervals))
  for(i in 1:length(x_intervals)){
    p_intervals[i] <- determine_k(x_intervals[i])
    upper_lower_intervals[i, 1] <- compute_expected_pi_for_k_NB(x_intervals[i], p_intervals[i], dispersion)
  }
  expected_line <- data.frame(aa = x_intervals, bb = logx_intervals, dd = upper_lower_intervals[, 1])


  if(log_scale == FALSE){
    require(ggplot2)
    g1 <- ggplot2::ggplot(stats_combined, ggplot2::aes(x = count_mean,
                                                       y = count_variance, col = category)) +
      ggplot2::geom_point(size = 2, alpha = 0.5) +
      ggplot2::geom_text(data = showhits, aes(label = rownames(showhits)), vjust = -1, size = 3, color = "brown") +
      ggplot2::theme(legend.position = "topright") + ggplot2::theme_bw() +
      ggplot2::ylab("Variance") + ggplot2::xlab("Gene Mean") + xlim(0, xmax) + ylim(0, 100)+
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     panel.grid = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     strip.placement = "inside")
    g2 <- ggplot2::ggplot(stats_combined, ggplot2::aes(x = count_mean,
                                                       y = count_cv, col = category)) +
      ggplot2::geom_point(size = 2, alpha = 0.5) +
      ggplot2::geom_text(data = showhits, aes(label = rownames(showhits)), vjust = -1, size = 3, color = "brown") +
      ggplot2::theme(legend.position = "topright") + ggplot2::theme_bw() +
      ggplot2::ylab("CV") + ggplot2::xlab("Gene Mean") + xlim(0, xmax) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     panel.grid = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     strip.placement = "inside")

    g3 <- ggplot2::ggplot(stats_combined, ggplot2::aes(x = lognorm_mean,
                                                       y = lognorm_cv2, col = category)) +
      ggplot2::geom_point(size = 2, alpha = 0.5) +
      ggplot2::geom_text(data = showhits, aes(label = rownames(showhits)), vjust = -1, size = 3, color = "brown") +
      ggplot2::theme(legend.position = "topright") + ggplot2::theme_bw() +
      ggplot2::ylab("lognormal CV2") + ggplot2::xlab("Gene Mean") + xlim(0, xmax) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     panel.grid = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     strip.placement = "inside")

    g_basic <- ggplot2::ggplot(stats_combined, ggplot2::aes(x = gene_mean,
                                                            y = k_proportion, col = category)) +
      ggplot2::geom_point(size = 2, alpha = 0.5) +
      ggplot2::geom_text(data = showhits, aes(label = rownames(showhits)), vjust = -1, size = 3, color = "brown") +
      ggplot2::theme(legend.position = "topright") + ggplot2::theme_bw() +
      ggplot2::ylab("K Proportion") + ggplot2::xlab("Gene Mean") + xlim(0, xmax) + ylim(0, 1)+
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     panel.grid = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     strip.placement = "inside") + ggplot2::geom_line(data = expected_line, aes(aa, dd), color = "#023373", show.legend = FALSE)


    g_basic <- g_basic +  ggplot2::scale_color_manual(values = c("B" = "black", "H" ="#F06060"),
                                                      labels = c("B" = "Background", "H" = "Inflated"))
    g1 <- g1 +  ggplot2::scale_color_manual(values = c("B" = "black", "H" ="#F06060"),
                                            labels = c("B" = "Background", "H" = "Inflated"))
    g2 <- g2 +  ggplot2::scale_color_manual(values = c("B" = "black", "H" ="#F06060"),
                                            labels = c("B" = "Background", "H" = "Inflated"))
    g3 <- g3 +  ggplot2::scale_color_manual(values = c("B" = "black", "H" ="#F06060"),
                                            labels = c("B" = "Background", "H" = "Inflated"))

    allg <- list(g_basic, g1, g2, g3)

  } else {

    stats_combined[, "gene_mean"] <- log10(stats_combined[, "gene_mean"]+0.1)
    showhits[, "gene_mean"] <- log10(showhits[, "gene_mean"]+0.1)
    showmore[, "gene_mean"] <- log10(showmore[, "gene_mean"]+0.1)
    xmax <- log_xmax #max(stats_combined[, "gene_mean"]) + 0.2

    require(ggplot2)
    g_basic <- ggplot2::ggplot(stats_combined, ggplot2::aes(x = gene_mean,
                                                            y = k_proportion, col = category)) +
      ggplot2::geom_point(size = 2, alpha = 0.5) +
      ggplot2::geom_point(data = showhits, size = 3, alpha = 0.5, color = "brown") +
      ggplot2::geom_text(data = showhits, aes(label = rownames(showhits)), vjust = -1, size = 3, color = "brown") +
      ggplot2::geom_point(data = showmore, size = 3, alpha = 0.5, color = "darkgreen") +
      ggplot2::geom_text(data = showmore, aes(label = rownames(showmore)), vjust = -1, size = 3, color = "darkgreen") +
      ggplot2::theme(legend.position = "topright") + ggplot2::theme_bw() +
      ggplot2::ylab("K Proportion") + ggplot2::xlab("log10 Gene Mean") + xlim(-1, xmax) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     panel.grid = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     strip.placement = "inside") + ggplot2::geom_line(data = expected_line, aes(bb, dd), color = "#023373", show.legend = FALSE)

    g1 <- ggplot2::ggplot(stats_combined, ggplot2::aes(x = gene_mean,
                                                       y = lognorm_variance, col = category)) +
      ggplot2::geom_point(size = 2, alpha = 0.5) +
      ggplot2::geom_point(data = showhits, size = 3, alpha = 0.5, color = "brown") +
      ggplot2::geom_text(data = showhits, aes(label = rownames(showhits)), vjust = -1, size = 3, color = "brown") +
      ggplot2::geom_point(data = showmore, size = 3, alpha = 0.5, color = "darkgreen") +
      ggplot2::geom_text(data = showmore, aes(label = rownames(showmore)), vjust = -1, size = 3, color = "darkgreen") +
      ggplot2::theme(legend.position = "topright") + ggplot2::theme_bw() +
      ggplot2::ylab("Variance") + ggplot2::xlab("log10 Gene Mean") + xlim(-1, xmax) + ylim(0, max(stats_combined$lognorm_variance))+
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     panel.grid = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     strip.placement = "inside")
    g2 <- ggplot2::ggplot(stats_combined, ggplot2::aes(x = gene_mean,
                                                       y = lognorm_cv, col = category)) +
      ggplot2::geom_point(size = 2, alpha = 0.5) +
      ggplot2::geom_point(data = showhits, size = 3, alpha = 0.5, color = "brown") +
      ggplot2::geom_text(data = showhits, aes(label = rownames(showhits)), vjust = -1, size = 3, color = "brown") +
      ggplot2::geom_point(data = showmore, size = 3, alpha = 0.5, color = "darkgreen") +
      ggplot2::geom_text(data = showmore, aes(label = rownames(showmore)), vjust = -1, size = 3, color = "darkgreen") +
      ggplot2::theme(legend.position = "topright") + ggplot2::theme_bw() +
      ggplot2::ylab("CV") + ggplot2::xlab("log10 Gene Mean") + xlim(-1, xmax) + ylim(0, 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     panel.grid = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     strip.placement = "inside")

    g3 <- ggplot2::ggplot(stats_combined, ggplot2::aes(x = gene_mean,
                                                       y = lognorm_cv2, col = category)) +
      ggplot2::geom_point(size = 2, alpha = 0.5) +
      ggplot2::geom_point(data = showhits, size = 3, alpha = 0.5, color = "brown") +
      ggplot2::geom_text(data = showhits, aes(label = rownames(showhits)), vjust = -1, size = 3, color = "brown") +
      ggplot2::geom_point(data = showmore, size = 3, alpha = 0.5, color = "darkgreen") +
      ggplot2::geom_text(data = showmore, aes(label = rownames(showmore)), vjust = -1, size = 3, color = "darkgreen") +
      ggplot2::theme(legend.position = "topright") + ggplot2::theme_bw() +
      ggplot2::ylab("lognormal CV2") + ggplot2::xlab("log10 Gene Mean") + xlim(-1, xmax) + ylim(0, 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     panel.grid = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     strip.placement = "inside")


    g_basic <- g_basic +  ggplot2::scale_color_manual(values = c("B" = "#5C4B51", "H" ="#4E64A6"),
                                                      labels = c("B" = "Background", "H" = "Inflated"))
    g1 <- g1 +  ggplot2::scale_color_manual(values = c("B" = "#5C4B51", "H" ="#4E64A6"),
                                            labels = c("B" = "Background", "H" = "Inflated"))
    g2 <- g2 +  ggplot2::scale_color_manual(values = c("B" = "#5C4B51", "H" ="#4E64A6"),
                                            labels = c("B" = "Background", "H" = "Inflated"))
    g3 <- g3 +  ggplot2::scale_color_manual(values = c("B" = "#5C4B51", "H" ="#4E64A6"),
                                            labels = c("B" = "Background", "H" = "Inflated"))

  }

  allg <- list(g_basic, g1, g2, g3)
  return(allg)

}


wave_one <- other_stats_plot_supp_figure2_override(cluster.zscore.summary[[j]], show_mean_cutoff = 1, log_scale = TRUE, show_top = 15, showgenes = high_var_genes)

library(easyGgplot2)
pdf(paste0("scRNAseq_", list.celltype_br[i], "_",  sample.info[as.character(cluster.zscore.summary[[j]]$groupID), 2], "_more_stats_suppfig2_highvar.pdf"), width = 12, height = 9)
ggplot2.multiplot(wave_one[[1]], wave_one[[2]], wave_one[[4]], cols=1)
dev.off()
########################
########## heatmap #####
########################
show_mean_cutoff <- 1
show_top <- 15
one.cluster.zscore.summary <- cluster.zscore.summary[[j]]
selected_stats <- one.cluster.zscore.summary$Inflated
subselected_stats <- selected_stats[selected_stats$gene_mean >= show_mean_cutoff, ]
showhits <- subselected_stats[1:show_top, ]
genes <- rownames(showhits)
save(genes, file = paste0("scRNAseq_", list.celltype_br[i], "_",  sample.info[as.character(cluster.zscore.summary[[j]]$groupID), 2], "_showhits.rds"))
subcounts <- testdata[genes, new.clusters==16]

library(viridis)
Colors=viridis(40)


subcounts[subcounts>20] <- 20
pdf(paste0("scRNAseq_", list.celltype_br[i], "_",  sample.info[as.character(cluster.zscore.summary[[j]]$groupID), 2], "_top15.pdf"), width=8, height = 8)
Heatmap(as.matrix(subcounts), name = "beta HI 25", col=Colors,  column_labels = rep("", ncol(subcounts)))
dev.off()



############################
########## violin_plot #####
############################
gene_violin_plot_override <- function(stats, ylab_text = NULL, xlab_text = NULL){

  require(ggplot2)
  g_basic <- ggplot2::ggplot(stats, ggplot2::aes(x = expression, y = gene)) +
    ggplot2::geom_violin(trim=FALSE, fill = "brown", col = "black") + #+  geom_boxplot(width=0.02) +
    ggplot2::theme(legend.position = "topright") + ggplot2::theme_bw() +
    ggplot2::ylab(ylab_text) + ggplot2::xlab(xlab_text) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                                                                         panel.grid = ggplot2::element_blank(),
                                                                         legend.title = ggplot2::element_blank(),
                                                                         axis.ticks = ggplot2::element_blank(),
                                                                         strip.placement = "inside")
  return(g_basic)
}


library(reshape2)
stats <- melt(as.matrix(subcounts), id = c("cell", "gene"))
colnames(stats) <- c("gene", "cell", "expression")
stats <- stats[stats$gene%in%c("GCG", "IAPP", "RBP4", "PCDH7", "PTPRD", "IER3", "PDE4B"), ]
pdf(paste0("scRNAseq_", list.celltype_br[i], "_",  sample.info[as.character(cluster.zscore.summary[[j]]$groupID), 2], "_top_violin.pdf"), width = 8, height = 8)
gene_violin_plot_override(stats, ylab_text = "Expression", xlab_text = "Gene Name")
dev.off()



############################
########## separate RNA #####
############################

i <- 14
j <- 16

selected <- which(group_labels %in% list.celltype[i])
testdata <- counts[, selected]

sub.donor.info <- donor.info[group_labels %in% list.celltype[i]]
new.clusters <- as.numeric(as.factor(sub.donor.info))

sample.info <- cbind(new.clusters, sub.donor.info)
sample.info <- sample.info[!duplicated(sample.info), ]
data("otherRNA")
cluster.zscore.summary <- NB_inflation_test_var_asymp_collecting_more_metrics(testdata, grouping = new.clusters, genemean.cutoff = 0, filter.gene = FALSE)





waveplot_fig2_override_RNA_hits_only <- function(one.cluster.zscore.summary, show_mean_cutoff = 1, show_top = 15,
                                                 xmax = 6){

  dispersion <- one.cluster.zscore.summary$NBdispersion
  selected_stats <- one.cluster.zscore.summary$Inflated

  score_mat <- one.cluster.zscore.summary$Zscores
  all_genes <- rownames(score_mat)

  selected.IDs <- which(all_genes %in% rownames(selected_stats))


  data("otherRNA")
  MT_genes <- otherRNA[otherRNA[, "Family"] == "Mitochondiral", "Gene"]
  RP_genes <- otherRNA[otherRNA[, "Family"] == "Ribosomal", "Gene"]
  IncRNA_genes <- otherRNA[otherRNA[, "Family"] == "IncRNA", "Gene"]

  MT.IDs <- which(all_genes%in%MT_genes)
  RP.IDs <- which(all_genes%in%RP_genes)
  IncRNA.IDs <- c(which(all_genes%in%IncRNA_genes), grep("^LINC", all_genes), grep("^A[A-Z][0-9][0-9][0-9][0-9][0-9][0-9]", all_genes))

  MT_stats <- score_mat[MT.IDs, ]
  RP_stats <- score_mat[RP.IDs, ]
  IncRNA_stats <- score_mat[IncRNA.IDs, ]

  background_stats <- score_mat[-c(selected.IDs, MT.IDs, RP.IDs, IncRNA.IDs), ]

  stats_combined <- data.frame(category = c(rep("B", nrow(background_stats)), rep("H", nrow(selected_stats)), rep("MT", nrow(MT_stats)),
                                            rep("RP", nrow(RP_stats)), rep("IncRNA", nrow(IncRNA_stats))), rbind(background_stats, selected_stats, MT_stats, RP_stats, IncRNA_stats))

  subselected_stats <- selected_stats[rownames(selected_stats)%in%all_genes[c(MT.IDs, RP.IDs, IncRNA.IDs)], ]
  subselected_stats <- subselected_stats[subselected_stats$gene_mean >= show_mean_cutoff, ]
  hits_num <- nrow(subselected_stats)
  show_top <- min(hits_num, show_top)
  showhits <- subselected_stats[1:show_top, ]



  x_intervals <- seq(0.1, xmax, by = 0.2)
  p_intervals <- rep(0, length(x_intervals))
  logx_intervals <- log10(x_intervals + 0.1)


  upper_lower_intervals <- matrix(0, ncol = 1, nrow = length(x_intervals))
  for(i in 1:length(x_intervals)){
    p_intervals[i] <- determine_k(x_intervals[i])
    upper_lower_intervals[i, 1] <- compute_expected_pi_for_k_NB(x_intervals[i], p_intervals[i], dispersion)
  }
  expected_line <- data.frame(aa = x_intervals, bb = logx_intervals, dd = upper_lower_intervals[, 1])

  require(ggplot2)
  g_basic <- ggplot2::ggplot(stats_combined, ggplot2::aes(x = gene_mean,
                                                          y = k_proportion, col = category)) +
    ggplot2::geom_point(size = 2, alpha = 0.5) +
    ggplot2::geom_text(data = showhits, aes(label = rownames(showhits)), vjust = -1, size = 3, color = "brown") +
    ggplot2::theme(legend.position = "topright") + ggplot2::theme_bw() +
    ggplot2::ylab("K Proportion") + ggplot2::xlab("Gene Mean") + xlim(0, xmax) + ylim(0, 1) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.grid = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   strip.placement = "inside") + ggplot2::geom_line(data = expected_line, aes(aa, dd), color = "#023373", show.legend = FALSE)

  g <- g_basic + ggplot2::scale_color_manual(values = c("B" = "#5C4B51", "IncRNA" = "#0099DD", "MT" = "#C2BB00", "RP" = "#8C1F66", "H" ="#4E64A6"),
                                             labels = c("B" = "Homeostatic", "IncRNA" = "IncRNA", "MT" = "Mt Genes", "RP" = "Ribo Genes", "H" = "Inflated"))

  return(g)

}
rownames(sample.info) <- sample.info[,1]
wave_one <- waveplot_fig2_override_RNA_hits_only(cluster.zscore.summary[[j]], show_mean_cutoff = 1, xmax = 10)
ggsave(paste0("scRNAseq_", list.celltype_br[i], "_",  sample.info[as.character(cluster.zscore.summary[[j]]$groupID), 2], "_otherRNA.pdf"), plot = wave_one, width = 12, height = 3)


hist(cluster.zscore.summary[[j]]$Inflated$k_proportion - cluster.zscore.summary[[j]]$Inflated$expected_pi)
table(cluster.zscore.summary[[j]]$Inflated$k_proportion - cluster.zscore.summary[[j]]$Inflated$expected_pi >= 0.05)
table(cluster.zscore.summary[[j]]$Inflated$k_proportion - cluster.zscore.summary[[j]]$Inflated$expected_pi >= 0.1)

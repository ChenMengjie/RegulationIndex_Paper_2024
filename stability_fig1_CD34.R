colors <- c("#5C4B51", "#8CBEB2", "#CC8D1A", "#F06060")
hits_color <- c("#4E64A6", "#829FD9")
line_colors <- c("#03658C", "#023373")



waveplot_fig1_override <- function(one.cluster.zscore.summary, show_mean_cutoff = 0.5, show_top = 25,
                                   xmax = 6, add_poisson_line = TRUE, dot_color, catenate = TRUE){

  dispersion <- one.cluster.zscore.summary$NBdispersion
  selected_stats <- one.cluster.zscore.summary$Inflated
  selected_stats <- selected_stats[selected_stats$gene_mean >= show_mean_cutoff, ]

  hits_num <- nrow(selected_stats)
  show_top <- min(hits_num, show_top)
  showhits <- selected_stats[1:show_top, ]

  score_mat <- one.cluster.zscore.summary$Zscores
  all_genes <- rownames(score_mat)

  selected.IDs <- which(all_genes %in% rownames(selected_stats))
  background_stats <- score_mat[-selected.IDs, ]

  stats_combined <- data.frame(category = c(rep("B", nrow(background_stats)), rep("H", nrow(selected_stats))),
                               rbind(background_stats, selected_stats))

  x_intervals <- seq(0.1, xmax, by = 0.05)
  p_intervals <- rep(0, length(x_intervals))
  logx_intervals <- log10(x_intervals + 0.1)


  upper_lower_intervals <- matrix(0, ncol = 1, nrow = length(x_intervals))
  for(i in 1:length(x_intervals)){
    p_intervals[i] <- determine_k(x_intervals[i])
    upper_lower_intervals[i, 1] <- compute_expected_pi_for_k_NB(x_intervals[i], p_intervals[i], dispersion)
  }
  expected_line <- data.frame(aa = x_intervals, bb = logx_intervals, cc = p_intervals, dd = upper_lower_intervals[, 1])


  Poisson_line <- data.frame(aap = x_intervals, bbp = logx_intervals, ccp = p_intervals, eep = ppois(p_intervals, x_intervals))

  require(ggplot2)
  g_basic <- ggplot2::ggplot(stats_combined, ggplot2::aes(x = gene_mean,
                                                          y = k_proportion, col = category)) +
    ggplot2::geom_point(size = 2, alpha = 0.5) +
    ggplot2::geom_point(data = showhits, size = 3, alpha = 0.5, color = "brown") +
  #       ggplot2::geom_text(data = showhits, aes(label = rownames(showhits)), vjust = -1, size = 3, color = "brown") +
    ggplot2::theme(legend.position = "topright") + ggplot2::theme_bw() +
    ggplot2::ylab("K Proportion") + ggplot2::xlab("Gene Mean") + xlim(0, xmax) + ylim(0.1, 1) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.grid = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   strip.placement = "inside")

  if(catenate == TRUE){
    g_basic <-  g_basic +  ggplot2::geom_line(data = expected_line, aes(aa, dd), color = "#023373", show.legend = FALSE)

  } else {
    for(jj in 0:max(p_intervals)){
      sub_expected_line <- expected_line[expected_line[, 3] == jj, ]
      g_basic <-  g_basic +  ggplot2::geom_line(data = sub_expected_line, aes(aa, dd), color = "#023373", show.legend = FALSE)
    }
  }



  g <- g_basic +  ggplot2::scale_color_manual(values = c("B" = dot_color, "H" ="brown"),
                                              labels = c("B" = "Homeostatic", "H" = "Inflated"))

  if(add_poisson_line == TRUE){

    if(catenate == TRUE){
      g <-  g +  ggplot2::geom_line(data = Poisson_line, aes(aap, eep), color = "#023373", linetype = "dashed", show.legend = FALSE)

    } else {
      for(jj in 0:max(p_intervals)){
        sub_Poisson_line <- Poisson_line[Poisson_line[, 3] == jj, ]
          g <-  g +  ggplot2::geom_line(data = sub_Poisson_line, aes(aap, eep), color = "#023373", linetype = "dashed", show.legend = FALSE)
      }
    }

  }

  return(g)

}


load("CD34.rds")
counts <- counts[, -which(duplicated(colnames(counts)))]
##############################
########## Wave plots ########
##############################


library(RegulationZIndex)
testdata <- t(counts)
cluster.zscore.summary <- NB_inflation_test_var_asymp(testdata, grouping = final_cluster, min.depth.group = 500, cell.depth.ranges = c(500, 10000), genemean.cutoff = 0.5)

wave_one <- waveplot_fig1_override(cluster.zscore.summary[[1]], show_mean_cutoff = 0.5, xmax = 10, add_poisson_line = TRUE, dot_color = "#CC8D1A")
ggsave(paste0("cRNAseq_CD34_group1_kprop_NB.pdf"), plot = wave_one, width = 12, height = 3)

wave_one <- waveplot_fig1_override(cluster.zscore.summary[[2]], show_mean_cutoff = 0.5, xmax = 10, add_poisson_line = TRUE, dot_color = "#8CBEB2")
ggsave(paste0("scRNAseq_CD34_group2_kprop_NB.pdf"), plot = wave_one, width = 12, height = 3)

wave_one <- waveplot_fig1_override(cluster.zscore.summary[[3]], show_mean_cutoff = 0.5, xmax = 10, add_poisson_line = TRUE, dot_color = "#F06060")
ggsave(paste0("scRNAseq_CD34_group3_kprop_NB.pdf"), plot = wave_one, width = 12, height = 3)


cluster.zscore.summary2 <- NB_inflation_test_var_asymp(testdata, min.depth.group = 500, cell.depth.ranges = c(500, 10000), genemean.cutoff = 0.5)
wave_one <- waveplot_fig1_override(cluster.zscore.summary2, show_mean_cutoff = 0.5, xmax = 10, add_poisson_line = TRUE, dot_color = "#5C4B51")
ggsave(paste0("scRNAseq_CD34_all_kprop_NB.pdf"), plot = wave_one, width = 12, height = 3)



########################
########## UMap ########
########################
final_colors <- final_cluster
final_colors[final_cluster == 2] <- "#CC8D1A"
final_colors[final_cluster == 3] <- "#8CBEB2"
final_colors[final_cluster == 4] <- "#F06060"
#load("/Users/mchen12/Downloads/more_Inflation_analysis/zheng_data/CD34.rds")
#plot(dimred[, 2], dimred[, 1], xlim = c(-8, 10), ylim = c(-5, 4), pch=20, col = final_colors)

pdf(paste0("umap_CD34.pdf"), width= 4, height = 4)
plot(dimred[, 2], dimred[, 1], xlim = c(-8, 10), ylim = c(-5, 4), pch=20, col = final_colors)
dev.off()



########################
########## heatmap #####
########################


genes <- rownames(cluster.zscore.summary2$Inflated)
genes <- c(genes, rownames(cluster.zscore.summary[[3]]$Inflated), rownames(cluster.zscore.summary[[2]]$Inflated), rownames(cluster.zscore.summary[[1]]$Inflated))

### show Poisson curves, estimated NB curves


groups <- c(rep("5", nrow(cluster.zscore.summary2$Inflated)),
            rep("4", nrow(cluster.zscore.summary[[3]]$Inflated)),
            rep("3", nrow(cluster.zscore.summary[[2]]$Inflated)),
            rep("2", nrow(cluster.zscore.summary[[1]]$Inflated)))
genes <- gsub("\\.", "-", genes)
genes[genes=="RP11-367G6-3"] <- "RP11-367G6.3"
genes[genes=="RP11-879F14-2"] <- "RP11-879F14.2"
subcounts <- t(counts[, genes])
library(ComplexHeatmap)
library(viridis)
Colors=viridis(40)
cell_type <- final_cluster
ha = HeatmapAnnotation(
  cell_type = cell_type,
  border = TRUE)

subcounts[subcounts>20] <- 20
pdf(paste0("scRNAseq_CD34_inflated_in_cells.pdf"), width=10, height = 22)
Heatmap(subcounts, name = "CD34+ inflated genes", cluster_row_slices = FALSE, cluster_column_slices = FALSE,
        column_split = factor(cell_type), row_split = factor(groups),
        top_annotation = ha,
        column_title_gp = gpar(fill = c( "#CC8D1A", "#8CBEB2", "#F06060")),
        col=Colors, row_title_gp = gpar(fill = c( "#CC8D1A", "#8CBEB2", "#F06060", "#5C4B51")))
dev.off()

library(Seurat)
subcounts <- t(counts[, genes])
test.normal <- NormalizeData(subcounts)

pdf(paste0("scRNAseq_CD34_inflated_in_cells_normalized.pdf"), width=10, height = 22)
Heatmap(as.matrix(test.normal), name = "CD34+ inflated genes", cluster_row_slices = FALSE, cluster_column_slices = FALSE,
        column_split = factor(cell_type), row_split = factor(groups),
        top_annotation = ha,
        column_title_gp = gpar(fill = c( "#CC8D1A", "#8CBEB2", "#F06060")),
        col=Colors, row_title_gp = gpar(fill = c( "#CC8D1A", "#8CBEB2", "#F06060", "#5C4B51")))
dev.off()

library_size <- apply(counts, 1, sum)
gene_coverage <- apply(counts, 1, function(x){sum(x != 0)})

pdf(paste0("libray_size_CD34.pdf"), width= 4, height = 4)
plot(gene_coverage, library_size, col = final_colors, pch = 20, ylim = c(0, 6000))
dev.off()





##############################
########## Wave plots ########
##############################

load("CD34.rds")
counts <- counts[, -which(duplicated(colnames(counts)))]

library(RegulationZIndex)
testdata <- t(counts)
cluster.zscore.summary <- NB_inflation_test_var_asymp(testdata, grouping = final_cluster, min.depth.group = 500, cell.depth.ranges = c(500, 10000), genemean.cutoff = 0.5)

wave1 <- waveplot_fig1_override(cluster.zscore.summary[[1]], show_mean_cutoff = 0.5, xmax = 10, add_poisson_line = TRUE, dot_color = "#CC8D1A")

wave2 <- waveplot_fig1_override(cluster.zscore.summary[[2]], show_mean_cutoff = 0.5, xmax = 10, add_poisson_line = TRUE, dot_color = "#8CBEB2")

wave3 <- waveplot_fig1_override(cluster.zscore.summary[[3]], show_mean_cutoff = 0.5, xmax = 10, add_poisson_line = TRUE, dot_color = "#F06060")

cluster.zscore.summary2 <- NB_inflation_test_var_asymp(testdata, min.depth.group = 500, cell.depth.ranges = c(500, 10000), genemean.cutoff = 0.5)
wave_one <- waveplot_fig1_override(cluster.zscore.summary2, show_mean_cutoff = 0.5, xmax = 10, add_poisson_line = TRUE, dot_color = "#5C4B51")

library(easyGgplot2)
pdf(paste0("all_CD34_wave_fig2.pdf"), width = 10, height = 12)
ggplot2.multiplot(wave1, wave2, wave3, wave_one, cols=1)
dev.off()






library(RegulationZIndex)
testdata <- t(counts)
cluster.zscore.summary <- NB_inflation_test_var_asymp(testdata, grouping = final_cluster, min.depth.group = 500, cell.depth.ranges = c(500, 10000), genemean.cutoff = 0.5)

wave1 <- waveplot_fig1_override(cluster.zscore.summary[[1]], show_mean_cutoff = 0.5, xmax = 10, add_poisson_line = TRUE, dot_color = "#CC8D1A")

wave2 <- waveplot_fig1_override(cluster.zscore.summary[[2]], show_mean_cutoff = 0.5, xmax = 10, add_poisson_line = TRUE, dot_color = "#8CBEB2")

wave3 <- waveplot_fig1_override(cluster.zscore.summary[[3]], show_mean_cutoff = 0.5, xmax = 10, add_poisson_line = TRUE, dot_color = "#F06060")

cluster.zscore.summary2 <- NB_inflation_test_var_asymp(testdata, min.depth.group = 500, cell.depth.ranges = c(500, 10000), genemean.cutoff = 0.5)
wave_one <- waveplot_fig1_override(cluster.zscore.summary2, show_mean_cutoff = 0.5, show_top = 120, xmax = 10, add_poisson_line = TRUE, dot_color = "#5C4B51")

library(easyGgplot2)
pdf(paste0("all_CD34_wave_fig2_no_text.pdf"), width = 10, height = 12)
ggplot2.multiplot(wave1, wave2, wave3, wave_one, cols=1)
dev.off()

> cluster.zscore.summary[[1]]$NBdispersion
[1] 0.3171141
> cluster.zscore.summary[[2]]$NBdispersion
[1] 0.5259099
> cluster.zscore.summary[[3]]$NBdispersion
[1] 0.163134
cluster.zscore.summary2$NBdispersion
[1] 1.38826



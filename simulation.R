set.seed(12334345)

pdf(file = "simulation_zscore_normality.pdf", width=12, height = 3)
par(mfrow = c(1, 4))

dispersion <- 1
gene_mean <- 10
k <- determine_k(gene_mean)
rrr <- 1/dispersion
rrr_gene_mean <- rrr + gene_mean
p_para <- rrr/rrr_gene_mean
dataX <- matrix(rnbinom(50*100000, size = rrr, prob = p_para), ncol = 1000)
test_typeI <- Estimate_dispersion_zscore_asymp_for_one_cluster(dataX)

sum(test_typeI$score_mat$pvalue<=0.05)/length(test_typeI$score_mat$pvalue)

#test_typeI_2 <- Estimate_dispersion_zscore_for_one_cluster(dataX)
#hist(test_typeI_2$score_mat$pvalue)
#sum(test_typeI_2$score_mat$pvalue<=0.05)/length(test_typeI_2$score_mat$pvalue)

aaa <- test_typeI$score_mat
hist(aaa$zscore, breaks = 50, prob = TRUE, main = "mean = 10, dispersion = 1, N = 1000", xlab = "Z-index")
lines(density(rnorm(10000, mean = 0, sd = 1)), col = "red")


########
dispersion <- 0.5
gene_mean <- 5
k <- determine_k(gene_mean)
rrr <- 1/dispersion
rrr_gene_mean <- rrr + gene_mean
p_para <- rrr/rrr_gene_mean
dataX <- matrix(rnbinom(50*100000, size = rrr, prob = p_para), ncol = 1000)
test_typeI <- Estimate_dispersion_zscore_asymp_for_one_cluster(dataX)
#hist(test_typeI$score_mat$pvalue)
sum(test_typeI$score_mat$pvalue<=0.05)/length(test_typeI$score_mat$pvalue)

#test_typeI_2 <- Estimate_dispersion_zscore_for_one_cluster(dataX)
#hist(test_typeI_2$score_mat$pvalue)
#sum(test_typeI_2$score_mat$pvalue<=0.05)/length(test_typeI_2$score_mat$pvalue)

aaa <- test_typeI$score_mat
hist(aaa$zscore, breaks = 50, prob = TRUE, main = "mean = 5, dispersion = 0.5, N = 1000", xlab = "Z-index")
lines(density(rnorm(10000, mean = 0, sd = 1)), col = "red")

########
dispersion <- 0.1
gene_mean <- 2
k <- determine_k(gene_mean)
rrr <- 1/dispersion
rrr_gene_mean <- rrr + gene_mean
p_para <- rrr/rrr_gene_mean
dataX <- matrix(rnbinom(50*100000, size = rrr, prob = p_para), ncol = 500)
test_typeI <- Estimate_dispersion_zscore_asymp_for_one_cluster(dataX)
#hist(test_typeI$score_mat$pvalue)
sum(test_typeI$score_mat$pvalue<=0.05)/length(test_typeI$score_mat$pvalue)

#test_typeI_2 <- Estimate_dispersion_zscore_for_one_cluster(dataX)
#hist(test_typeI_2$score_mat$pvalue)
#sum(test_typeI_2$score_mat$pvalue<=0.05)/length(test_typeI_2$score_mat$pvalue)

aaa <- test_typeI$score_mat
hist(aaa$zscore, breaks = 50, prob = TRUE, main = "mean = 2, dispersion = 0.1, N = 500", xlab = "Z-index")
lines(density(rnorm(10000, mean = 0, sd = 1)), col = "red")


set.seed(12334345)



dispersion <- 0.5
gene_mean <- 0.1
k <- determine_k(gene_mean)
rrr <- 1/dispersion
rrr_gene_mean <- rrr + gene_mean
p_para <- rrr/rrr_gene_mean
dataX <- matrix(rnbinom(50*100000, size = rrr, prob = p_para), ncol = 500)
test_typeI <- Estimate_dispersion_zscore_asymp_for_one_cluster(dataX)
sum(test_typeI$score_mat$pvalue<=0.05)/length(test_typeI$score_mat$pvalue)

aaa <- test_typeI$score_mat
hist(aaa$zscore, breaks = 50, prob = TRUE, main = "mean = 0.1, dispersion = 0.5, N = 500", xlab = "Z-index")
lines(density(rnorm(10000, mean = 0, sd = 1)), col = "red")
dev.off()




set.seed(12334345)

pdf(file = "simulation_pvalue.pdf", width=12, height = 3)
par(mfrow = c(1, 4))

dispersion <- 1
gene_mean <- 10
k <- determine_k(gene_mean)
rrr <- 1/dispersion
rrr_gene_mean <- rrr + gene_mean
p_para <- rrr/rrr_gene_mean
dataX <- matrix(rnbinom(50*100000, size = rrr, prob = p_para), ncol = 1000)
test_typeI <- Estimate_dispersion_zscore_asymp_for_one_cluster(dataX)
hist(test_typeI$score_mat$pvalue, main = "mean = 10, dispersion = 1, N = 1000", xlab = "p-value")
sum(test_typeI$score_mat$pvalue<=0.05)/length(test_typeI$score_mat$pvalue)


########
dispersion <- 0.5
gene_mean <- 5
k <- determine_k(gene_mean)
rrr <- 1/dispersion
rrr_gene_mean <- rrr + gene_mean
p_para <- rrr/rrr_gene_mean
dataX <- matrix(rnbinom(50*100000, size = rrr, prob = p_para), ncol = 1000)
test_typeI <- Estimate_dispersion_zscore_asymp_for_one_cluster(dataX)
hist(test_typeI$score_mat$pvalue, main = "mean = 5, dispersion = 0.5, N = 1000", xlab = "p-value")


########
dispersion <- 0.1
gene_mean <- 2
k <- determine_k(gene_mean)
rrr <- 1/dispersion
rrr_gene_mean <- rrr + gene_mean
p_para <- rrr/rrr_gene_mean
dataX <- matrix(rnbinom(50*100000, size = rrr, prob = p_para), ncol = 500)
test_typeI <- Estimate_dispersion_zscore_asymp_for_one_cluster(dataX)
hist(test_typeI$score_mat$pvalue, main = "mean = 2, dispersion = 0.1, N = 500", xlab = "p-value", breaks = 10)



set.seed(12334345)

dispersion <- 0.5
gene_mean <- 0.1
k <- determine_k(gene_mean)
rrr <- 1/dispersion
rrr_gene_mean <- rrr + gene_mean
p_para <- rrr/rrr_gene_mean
dataX <- matrix(rnbinom(50*100000, size = rrr, prob = p_para), ncol = 500)
test_typeI <- Estimate_dispersion_zscore_asymp_for_one_cluster(dataX)
hist(test_typeI$score_mat$pvalue, main = "mean = 0.1, dispersion = 0.5, N = 500", xlab = "p-value", breaks = 10)

sum(test_typeI$score_mat$pvalue<=0.05)/length(test_typeI$score_mat$pvalue)

dev.off()




########################
########################
######## Power #########
########################


outlier_val <- c(2, 4, 8, 10)
per_val <- c(0.01, 0.02, 0.05, 0.1)*100




dispersion <- 0.5
gene_mean <- 0.1
k <- determine_k(gene_mean)
rrr <- 1/dispersion
rrr_gene_mean <- rrr + gene_mean
p_para <- rrr/rrr_gene_mean

summary_p_list <- NULL
for(kkkk in 1:20){
for(per in per_val){

  for(outer in outlier_val){

    dataX <- matrix(rnbinom(200*(100-per), size = rrr, prob = p_para), ncol = (100-per))
    dataY <- matrix(rep(outer, 200*per), ncol = per)
    dataXY <- cbind(dataX, dataY)
    test_powerI <- Estimate_dispersion_zscore_asymp_for_one_cluster(dataXY)

    dataZ <- matrix(rnbinom(200*100, size = rrr, prob = p_para), ncol = 100)
    test_typeI <- Estimate_dispersion_zscore_asymp_for_one_cluster(dataZ)
    summary_p_list <- rbind(summary_p_list, c(kkkk, per, outer, sum(test_powerI$score_mat$pvalue<=0.05)/length(test_powerI$score_mat$pvalue), sum(test_typeI$score_mat$pvalue<=0.05)/length(test_typeI$score_mat$pvalue)))

  }

}
}

par(mfrow = c(2, 2))
plot(summary_p_list[summary_p_list[, 2] == 1, 4], summary_p_list[summary_p_list[, 2] == 1, 5], col = as.numeric(as.factor(summary_p_list[summary_p_list[, 2] == 1, 3])), xlab = "Power", ylab = "Type I error", ylim = c(0, 0.1), xlim = c(0, 1), main = "1% outlier", pch = 20)
#legend("bottomleft", pch = 20, col = 1:4, c("Outlier value = 2", "Outlier value = 4", "Outlier value = 6", "Outlier value = 8"))
plot(summary_p_list[summary_p_list[, 2] == 5, 4], summary_p_list[summary_p_list[, 2] == 2, 5], col = as.numeric(as.factor(summary_p_list[summary_p_list[, 2] == 2, 3])), xlab = "Power", ylab = "Type I error", ylim = c(0, 0.1), xlim = c(0, 1), main = "2% outlier", pch = 20)
#legend("bottomleft", pch = 20, col = 1:4, c("Outlier value = 2", "Outlier value = 4", "Outlier value = 6", "Outlier value = 8"))
plot(summary_p_list[summary_p_list[, 2] == 5, 4], summary_p_list[summary_p_list[, 2] == 5, 5], col = as.numeric(as.factor(summary_p_list[summary_p_list[, 2] == 5, 3])), xlab = "Power", ylab = "Type I error", ylim = c(0, 0.1), xlim = c(0, 1), main = "5% outlier", pch = 20)
#legend("bottomleft", pch = 20, col = 1:4, c("Outlier value = 2", "Outlier value = 4", "Outlier value = 6", "Outlier value = 8"))
plot(summary_p_list[summary_p_list[, 2] == 10, 4], summary_p_list[summary_p_list[, 2] == 10, 5], col = as.numeric(as.factor(summary_p_list[summary_p_list[, 2] == 10, 3])), xlab = "Power", ylab = "Type I error", ylim = c(0, 0.1), xlim = c(0, 1), main = "10% outlier", pch = 20)
#legend("bottomleft", pch = 20, col = 1:4, c("Outlier value = 2", "Outlier value = 4", "Outlier value = 6", "Outlier value = 8"))


colnames(summary_p_list) <- c("rep", "per", "outer", "power", "type1")

df <- as.data.frame(summary_p_list)
df$power <- df$power*100
df$type1 <- df$type1*100
df$per <- paste0("Outlier percentage = ", df$per, "%")
df$per = factor(df$per, levels=c('Outlier percentage = 1%','Outlier percentage = 2%','Outlier percentage = 5%','Outlier percentage = 10%'))


library(ggplot2)
library(easyGgplot2)

ff <- ggplot2.scatterplot(data = df, xName = 'power', yName = 'type1', backgroundColor = "white", xtitle="Power", ytitle="Type I error",
                          mainTitle = "Gene mean = 0.1, dispersion = 0.5", removePanelGrid=TRUE, groupName = "outer", removePanelBorder=FALSE, showLegend=TRUE,
                          legendTitle="Outlier value", legendTitleFont=c(15, "bold", "black"), legendTextFont=c(15, "bold", "black"),
                          mainTitleFont = c(15, "bold", "black"), xtitleFont = c(15, "bold", "black"),  ytitleFont = c(15, "bold", "black"),
                          xTickLabelFont = c(15, "bold", "black"), yTickLabelFont = c(15, "bold", "black"), facetingFont = c(15, "bold", "black"),
                          facetingRect =  list(background = NULL, lineType = NULL, lineColor = NULL, lineSize = NULL)) + ylim(0, 10)  + xlim(0, 100) + facet_wrap(~per, ncol=2)  + theme(strip.text.x = element_text(size = 15, colour = "black", face = "bold"))

ggsave("power_mean01.pdf", plot = ff, width = 8, height = 8)
save(summary_p_list, file = "figures_in_man/Supp/power_mean01.rdata")









outlier_val <- c(2, 4, 8, 10)
per_val <- c(0.01, 0.02, 0.05, 0.1)*100


dispersion <- 0.5
gene_mean <- 0.5
k <- determine_k(gene_mean)
rrr <- 1/dispersion
rrr_gene_mean <- rrr + gene_mean
p_para <- rrr/rrr_gene_mean

summary_p_list <- NULL
for(kkkk in 1:20){
  for(per in per_val){

    for(outer in outlier_val){

      dataX <- matrix(rnbinom(200*(100-per), size = rrr, prob = p_para), ncol = (100-per))
      dataY <- matrix(rep(outer, 200*per), ncol = per)
      dataXY <- cbind(dataX, dataY)
      test_powerI <- Estimate_dispersion_zscore_asymp_for_one_cluster(dataXY)

      dataZ <- matrix(rnbinom(200*100, size = rrr, prob = p_para), ncol = 100)
      test_typeI <- Estimate_dispersion_zscore_asymp_for_one_cluster(dataZ)
      summary_p_list <- rbind(summary_p_list, c(kkkk, per, outer, sum(test_powerI$score_mat$pvalue<=0.05)/length(test_powerI$score_mat$pvalue), sum(test_typeI$score_mat$pvalue<=0.05)/length(test_typeI$score_mat$pvalue)))

    }

  }
}


colnames(summary_p_list) <- c("rep", "per", "outer", "power", "type1")

df <- as.data.frame(summary_p_list)
df$power <- df$power*100
df$type1 <- df$type1*100
df$per <- paste0("Outlier percentage = ", df$per, "%")
df$per = factor(df$per, levels=c('Outlier percentage = 1%','Outlier percentage = 2%','Outlier percentage = 5%','Outlier percentage = 10%'))


library(ggplot2)
library(easyGgplot2)

ff <- ggplot2.scatterplot(data = df, xName = 'power', yName = 'type1', backgroundColor = "white", xtitle="Power", ytitle="Type I error",
                          mainTitle = "Gene mean = 0.5, dispersion = 0.5", removePanelGrid=TRUE, groupName = "outer", removePanelBorder=FALSE, showLegend=TRUE,
                          legendTitle="Outlier value", legendTitleFont=c(15, "bold", "black"), legendTextFont=c(15, "bold", "black"),
                          mainTitleFont = c(15, "bold", "black"), xtitleFont = c(15, "bold", "black"),  ytitleFont = c(15, "bold", "black"),
                          xTickLabelFont = c(15, "bold", "black"), yTickLabelFont = c(15, "bold", "black"), facetingFont = c(15, "bold", "black"),
                          facetingRect =  list(background = NULL, lineType = NULL, lineColor = NULL, lineSize = NULL)) + ylim(0, 10)  + xlim(0, 100) + facet_wrap(~per, ncol=2)  + theme(strip.text.x = element_text(size = 15, colour = "black", face = "bold"))

ggsave("power_mean05.pdf", plot = ff, width = 8, height = 8)
save(summary_p_list, file = "power_mean05.rdata")






########### Waveplot



library(RegulationZIndex)
xmax <- 8
dispersion <- 0.25
x_intervals <- seq(0.1, xmax, by = 0.2)


line0 <- compute_expected_pi_for_k_NB(x_intervals, 0, dispersion)
line1 <- compute_expected_pi_for_k_NB(x_intervals, 1, dispersion)
line2 <- compute_expected_pi_for_k_NB(x_intervals, 2, dispersion)
line3 <- compute_expected_pi_for_k_NB(x_intervals, 3, dispersion)
line4 <- compute_expected_pi_for_k_NB(x_intervals, 4, dispersion)



p_intervals <- rep(0, length(x_intervals))
   for(i in 1:length(x_intervals)){
      p_intervals[i] <- determine_k(x_intervals[i])
      upper_lower_intervals[i, 1] <- compute_expected_pi_for_k_NB(x_intervals[i], p_intervals[i], dispersion)
      upper_lower_intervals[i, 2] <- solve_for_pi_NB(x_intervals[i], dispersion, n = n_CI, one_sided = one_sided)
    }
    expected_line <- data.frame(aa = x_intervals, bb = logx_intervals, dd = upper_lower_intervals[, 1])
    CI_line <- data.frame(aa = x_intervals, bb = logx_intervals, ee = upper_lower_intervals[, 2])
    
    Poisson_line <- data.frame(aap = x_intervals, bbp = logx_intervals, eep = ppois(p_intervals, x_intervals))

 
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
                   strip.placement = "inside") + ggplot2::geom_line(data = expected_line, aes(aa, dd), color = "turquoise4", show.legend = FALSE)



par(mfrow = c(1, 4))

all_prop_stats <- NULL
dispersion <- 0.25
rrr <- 1/dispersion
for(gene_mean in x_intervals){
	rrr_gene_mean <- rrr + gene_mean
	p_para <- rrr/rrr_gene_mean
	dataX <- rnbinom(1000, size = rrr, prob = p_para)
	props <- c(length(dataX[dataX<=0]), length(dataX[dataX<=1]), length(dataX[dataX<=2]), length(dataX[dataX<=3]), length(dataX[dataX<=4]))
	all_prop_stats <- rbind(all_prop_stats, props)
}


stats_combined <- data.frame(gene_mean = rep(x_intervals, 5), obs = c(all_prop_stats)/1000, k_proportion = c(line0, line1, line2, line3, line4), category = c(rep("0", 40), rep("1", 40), rep("2", 40), rep("3", 40), rep("4", 40)))
require(ggplot2)
library("viridis")  
  g_basic <- ggplot2::ggplot(stats_combined, ggplot2::aes(x = gene_mean,
                                                     y = k_proportion)) +
      ggplot2::geom_line(aes(color = category)) +
    #  ggplot2::geom_point(aes(x = gene_mean, y = obs), size = 2, alpha = 0.25) +
      ggplot2::theme(legend.position = "topright") + ggplot2::theme_bw() +
     ggplot2::ylab("Proportion") + ggplot2::xlab("Gene Mean") + xlim(0, xmax) + ylim(0, 1) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.grid = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   strip.placement = "inside") # + scale_color_viridis(discrete = TRUE, option = "D")
     
     g_basic  
    
    
   x_intervals <- seq(0.1, xmax, by = 0.01)  
   upper_lower_intervals <- matrix(0, ncol = 2, nrow = length(x_intervals))
 
   p_intervals <- rep(0, length(x_intervals))
   for(i in 1:length(x_intervals)){
      p_intervals[i] <- determine_k(x_intervals[i])
      upper_lower_intervals[i, 1] <- compute_expected_pi_for_k_NB(x_intervals[i], p_intervals[i], dispersion)
      upper_lower_intervals[i, 2] <- solve_for_pi_NB(x_intervals[i], dispersion, n = 100, one_sided = TRUE)
    }
    expected_line <- data.frame(aa = x_intervals, pp = as.character(p_intervals), dd = upper_lower_intervals[, 1])
    CI_line <- data.frame(aa = x_intervals, ee = upper_lower_intervals[, 2])
  
  
   g_basic2 <- ggplot2::ggplot(expected_line, ggplot2::aes(x = aa,
                                                     y = dd)) +
      ggplot2::geom_line(aes(color = pp)) +
     # ggplot2::geom_point(aes(x = gene_mean, y = obs), size = 2, alpha = 0.25) +
      ggplot2::theme(legend.position = "topright") + ggplot2::theme_bw() +
     ggplot2::ylab("Proportion") + ggplot2::xlab("Gene Mean") + xlim(0, xmax) + ylim(0, 1) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.grid = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   strip.placement = "inside") + ggplot2::geom_line(data = CI_line, aes(aa, ee), color = "turquoise4", linetype = "dashed", show.legend = FALSE)
     
     g_basic2  
  
  
     
     
     
    ggplot2::geom_line(data = expected_line, aes(aa, dd), color = "turquoise4", show.legend = FALSE)



##############

library(RegulationZIndex)

dispersion <- 0.25
gene_mean <- 0.5
k <- determine_k(gene_mean)
rrr <- 1/dispersion
rrr_gene_mean <- rrr + gene_mean
p_para <- rrr/rrr_gene_mean


outlier_val <- c(10, 30, 50, 70)
per_val <- c(0.01, 0.02, 0.05, 0.1)*100


all.data <- NULL

 for(per in per_val){
    for(outer in outlier_val){

      dataX <- matrix(rnbinom((100-per), size = rrr, prob = p_para), ncol = (100-per))
      dataY <- matrix(rep(outer, per), ncol = per)
      dataXY <- cbind(dataX, dataY)
      all.data <- rbind(all.data, dataXY)
     
   }
}

Heatmap(t(all.data), cluster_rows = FALSE, cluster_columns = FALSE, col = viridis(40))
summary_stats <- apply(all.data, 1, function(x){
	mean_x <- mean(x)
	de_k <- determine_k(mean_x)
	c(mean_x, length(x[x<=0])/100)
})

obs_points <- as.data.frame(t(summary_stats))
colnames(obs_points) <- c("gene_mean", "obs")

out_info <- data.frame(outlier = rep(c(10, 30, 50, 70), each = 4), percentage = rep(per_val, 4))
obs_points <- cbind(obs_points, out_info)
g_basic2 + ggplot2::geom_point(data = obs_points, aes(x = gene_mean, y = obs), size = 2, alpha = 0.25)


##########
    
    outlier_val <- c(2, 4, 8, 10)
    per_val <- c(0.01, 0.02, 0.05, 0.1)*100
    
    
    dispersion <- 0.5
    gene_mean <- 0.5
    k <- determine_k(gene_mean)
    rrr <- 1/dispersion
    rrr_gene_mean <- rrr + gene_mean
    p_para <- rrr/rrr_gene_mean
    
    wave1 <- list(NULL)
    wave2 <- list(NULL)
    index_per <- 1
    for(per in per_val){
      
      all_gene_mean <- NULL
      all_k_proportion <- NULL
      pvalues <- NULL
      color.bar <- NULL
      for(outer in outlier_val){
        dataX <- matrix(rnbinom(200*(100-per), size = rrr, prob = p_para), ncol = (100-per))
        dataY <- matrix(rep(outer, 200*per), ncol = per)
        dataXY <- cbind(dataX, dataY)
        test_powerI <- Estimate_dispersion_zscore_asymp_for_one_cluster(dataXY)
        all_gene_mean <- c(all_gene_mean, test_powerI$score_mat$gene_mean)
        all_k_proportion <- c(all_k_proportion, test_powerI$score_mat$k_proportion)  
        pvalues <- c(pvalues, test_powerI$score_mat$pvalue)   
        color.bar <- c(color.bar, rep(outer, 200))
      }
      
    wave1[[index_per]] <- waveplot_simplified(all_gene_mean, all_k_proportion, pvalues, 0.5)
    wave2[[index_per]]  <- waveplot_simulation(all_gene_mean, all_k_proportion, 0.5, color.bar = color.bar)
      
    index_per <- index_per + 1
    }
    
    
    library(easyGgplot2)
    pdf("simulation_waveplot.pdf", width = 12, height = 8)
    ggplot2.multiplot(wave2[[1]], wave1[[1]], wave2[[2]], wave1[[2]], wave2[[3]], wave1[[3]], wave2[[4]], wave1[[4]], cols=2)
    dev.off()
    
    
    
    outlier_val <- c(2, 4, 8, 10)
    per_val <- c(0.01, 0.02, 0.05, 0.1)*100
    
    
    dispersion <- 0.5
    gene_mean <- 0.1
    k <- determine_k(gene_mean)
    rrr <- 1/dispersion
    rrr_gene_mean <- rrr + gene_mean
    p_para <- rrr/rrr_gene_mean
    
    wave1 <- list(NULL)
    wave2 <- list(NULL)
    index_per <- 1
    for(per in per_val){
      
      all_gene_mean <- NULL
      all_k_proportion <- NULL
      pvalues <- NULL
      color.bar <- NULL
      for(outer in outlier_val){
        dataX <- matrix(rnbinom(200*(100-per), size = rrr, prob = p_para), ncol = (100-per))
        dataY <- matrix(rep(outer, 200*per), ncol = per)
        dataXY <- cbind(dataX, dataY)
        test_powerI <- Estimate_dispersion_zscore_asymp_for_one_cluster(dataXY)
        all_gene_mean <- c(all_gene_mean, test_powerI$score_mat$gene_mean)
        all_k_proportion <- c(all_k_proportion, test_powerI$score_mat$k_proportion)  
        pvalues <- c(pvalues, test_powerI$score_mat$pvalue)   
        color.bar <- c(color.bar, rep(outer, 200))
      }
      
      wave1[[index_per]] <- waveplot_simplified(all_gene_mean, all_k_proportion, pvalues, 0.5)
      wave2[[index_per]]  <- waveplot_simulation(all_gene_mean, all_k_proportion, 0.5, color.bar = color.bar)
      
      index_per <- index_per + 1
    }
    
    
    library(easyGgplot2)
    pdf("imulation01_waveplot.pdf", width = 12, height = 8)
    ggplot2.multiplot(wave2[[1]], wave1[[1]], wave2[[2]], wave1[[2]], wave2[[3]], wave1[[3]], wave2[[4]], wave1[[4]], cols=2)
    dev.off()
    
    
    
    
    
    


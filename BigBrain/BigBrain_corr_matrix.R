
# Correlations between BigBrain cell density at different depths and MRI markers

markers = c("BSC", "GWC", "CT", "t1t2_ratio_25_GM", "t1t2_ratio_50_GM", "t1t2_ratio_25_WM")
markers_num = seq(1, length(markers))
BB_depths = c("BB25GM", "BB50GM", "BB25WM")
BB_depths_num = seq(1, length(markers))

results_spin = as.data.frame(matrix(ncol=6))
colnames(results_spin) = c("BB_depth", "BB_depth_num", "marker", "marker_num", "pval", "rho")

for (m in 1:length(markers)){
  for (d in 1:length(BB_depths)){
    pval = read.table(paste0("./spin_test/", BB_depths[d],"_vs_", markers[m],"_pval.txt"))
    rho = read.table(paste0("./spin_test/", BB_depths[d],"_vs_", markers[m],"_rho.txt"))
    
    results_spin[nrow(results_spin)+1,] = c(BB_depths[d], BB_depths_num[d], markers[m], markers_num[m], pval[1,1], rho[1,1])
  }
}

#FDR correction
results_spin = results_spin[-1,]
results_spin$pval_fdr = p.adjust(results_spin$pval, method='fdr')
results_spin$pval = as.numeric(results_spin$pval)
results_spin$rho = as.numeric(results_spin$rho)
results_spin$BB_depth = as.factor(results_spin$BB_depth)
results_spin$BB_depth_num = as.factor(results_spin$BB_depth_num)
results_spin$marker = as.factor(results_spin$marker)
results_spin$marker_num = as.factor(results_spin$marker_num)

# col with only sig pvalues and rho
results_spin$pval_sig = NA
results_spin$pval_fdr_sig = NA
results_spin$rho_sig = NA
results_spin$rho_fdr_sig = NA

sig_vars = as.data.frame(matrix(ncol=2))
colnames(sig_vars) = c("BB_depth_num", "marker_num")

sig_vars_fdr = as.data.frame(matrix(ncol=2))
colnames(sig_vars_fdr) = c("BB_depth_num", "marker_num")

for (i in 1:nrow(results_spin)){
  if (results_spin$pval[i] <= 0.05) {
    results_spin$pval_sig[i] = as.numeric(results_spin$pval[i])
    results_spin$rho_sig[i] = as.numeric(results_spin$rho[i])
    
    sig_vars[nrow(sig_vars)+1,] = c(as.numeric(results_spin$BB_depth_num[i]), as.numeric(results_spin$marker_num[i]))
  }
  if (results_spin$pval_fdr[i] <= 0.05) {
    results_spin$pval_fdr_sig[i] = as.numeric(results_spin$pval_fdr[i])
    results_spin$rho_fdr_sig[i] = as.numeric(results_spin$rho[i])
    
    sig_vars_fdr[nrow(sig_vars_fdr)+1,] = c(as.numeric(results_spin$BB_depth_num[i]), as.numeric(results_spin$marker_num[i]))
  }
}

sig_vars = sig_vars[-1,]
sig_vars_fdr = sig_vars_fdr[-1,]


# Heat map of associations
library(ggplot2)

ggplot(results_spin, aes(x=marker_num, y=BB_depth_num)) + 
  geom_tile(aes(fill=as.numeric(rho)), color = "white", lwd = 0.1, linetype = 1) + 
  scale_fill_gradientn(colours=c("#0000FFFF", "#FFFFFFFF", "#FF0000FF"), limits=c(-1,1), breaks=c(-1, -0.5, 0, 0.5, 1), name="Pearson's r", oob = scales::squish) + 
  geom_text(aes(label=paste0("r=",round(as.numeric(rho), 2), '\np=',round(as.numeric(pval), 2))), color='black', size=6) +
  geom_rect(data=sig_vars, mapping=aes(xmin=marker_num-0.5, xmax=marker_num+0.5, ymin=BB_depth_num-0.5, ymax=BB_depth_num+0.5), size=2, fill=NA, color="green3") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.text = element_text(size=12), legend.title = element_text(size=15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + 
  scale_x_discrete(labels=c(markers)) + 
  scale_y_discrete(labels=c(BB_depths)) + 
  ggsave("./all_corr_rho.jpg", height=4.5, width=8)

ggplot(results_spin, aes(x=marker_num, y=BB_depth_num)) + 
  geom_tile(aes(fill=as.numeric(rho)), color = "white", lwd = 0.1, linetype = 1) + 
  scale_fill_gradientn(colours=c("#0000FFFF", "#FFFFFFFF", "#FF0000FF"), limits=c(-1,1), breaks=c(-1, -0.5, 0, 0.5, 1), name="Pearson's r", oob = scales::squish) + 
  geom_text(aes(label=paste0("r=",round(as.numeric(rho), 2), '\np=',round(as.numeric(pval_fdr), 2))), color='black', size=6) +
  geom_rect(data=sig_vars_fdr, mapping=aes(xmin=marker_num-0.5, xmax=marker_num+0.5, ymin=BB_depth_num-0.5, ymax=BB_depth_num+0.5), size=2, fill=NA, color="green3") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.text = element_text(size=12), legend.title = element_text(size=15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + 
  scale_x_discrete(labels=c(markers)) + 
  scale_y_discrete(labels=c(BB_depths)) + 
  ggsave("./all_corr_rho_fdr.jpg", height=4.5, width=8)

# Write results to csv
write.csv(results_spin, "./all_results.csv", row.names=FALSE, col.names=TRUE)


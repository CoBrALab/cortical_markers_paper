
# Correlations between MRI markers and AHBA cell-type data

markers = c("BSC", "GWC", "CT", "t1t2_ratio_25", "t1t2_ratio_50", "t1t2_ratio_25_WM")
markers_num = seq(1, length(markers))
cells = c("astro", "endo", "micro", "neuroex", "neuroin", "oligo", "opc")
cells_num = seq(1, length(cells))

results_spin = as.data.frame(matrix(ncol=6))
colnames(results_spin) = c("marker", "marker_num", "cell", "cell_num", "pval", "rho")

for (m in 1:length(markers)){
  for (c in 1:length(cells)){
    pval = read.table(paste0("./spin_test_vertex_wise/means_", markers[m],"_", cells[c],"_pval.txt"))
    rho = read.table(paste0("./spin_test_vertex_wise/means_", markers[m],"_", cells[c],"_rho.txt"))
    
    results_spin[nrow(results_spin)+1,] = c(markers[m], markers_num[m], cells[c], cells_num[c], pval[1,1], rho[1,1])
  }
}

# Get thresholds for visualization of AHBA cell-type maps
cell_type_thresh = as.data.frame(matrix(ncol=7))
colnames(cell_type_thresh) = c("cell-type", "-1 SD", "+1 SD", "-1.5 SD", "+1.5 SD", "-2 SD", "+2 SD")

for (c in 1:length(cells)){
  print(cells[c])
  vals = as.data.frame(read.table(paste0("./AHBA_cell_type/cell_classes_AHBA_lh152parc_", cells[c],".txt")))
  mean = mean(vals[,1])
  sd = sd(vals[,1])
  
  cell_type_thresh[nrow(cell_type_thresh)+1,] = c(cells[c], (mean-sd), (mean+sd), (mean-(1.5*sd)), (mean+(1.5*sd)), (mean-(2*sd)), (mean+(2*sd)))
}
cell_type_thresh = cell_type_thresh[-1,]
write.csv(cell_type_thresh, './cell_type_thresh.csv', col.names = TRUE, row.names = FALSE)


#FDR correction
results_spin = results_spin[-1,]
results_spin$pval_fdr = p.adjust(results_spin$pval, method='fdr')
results_spin$pval = as.numeric(results_spin$pval)

# col with only sig pvalues and rho
results_spin$pval_sig = NA
results_spin$pval_fdr_sig = NA
results_spin$rho_sig = NA
results_spin$rho_fdr_sig = NA

sig_vars = as.data.frame(matrix(ncol=2))
colnames(sig_vars) = c("cell_num", "marker_num")

for (i in 1:nrow(results_spin)){
  if (results_spin$pval[i] <= 0.05) {
    results_spin$pval_sig[i] = as.numeric(results_spin$pval[i])
    results_spin$rho_sig[i] = as.numeric(results_spin$rho[i])
    
    sig_vars[nrow(sig_vars)+1,] = c(as.numeric(results_spin$cell_num[i]), as.numeric(results_spin$marker_num[i]))
  }
  if (results_spin$pval_fdr[i] <= 0.05) {
    results_spin$pval_fdr_sig[i] = as.numeric(results_spin$pval_fdr[i])
    results_spin$rho_fdr_sig[i] = as.numeric(results_spin$rho[i])
  }
}

sig_vars = sig_vars[-1,]

# Heat map of associations
library(ggplot2)

results_spin$marker = as.factor(results_spin$marker)
results_spin$cell = as.factor(results_spin$cell)

ggplot(results_spin, aes(x=marker_num, y=cell_num)) + 
  geom_tile(aes(fill=as.numeric(rho)), color = "white", lwd = 0.1, linetype = 1) + 
  scale_fill_gradientn(colours=c("#0000FFFF", "#FFFFFFFF", "#FF0000FF"), limits=c(-1,1), breaks=c(-1, 0, 1), name="Pearson's R", oob = scales::squish) + 
  geom_text(aes(label=paste0("r=",round(as.numeric(rho), 2), '\np=',round(as.numeric(pval), 3))), color='black', size=5) +
  geom_rect(data=sig_vars, mapping=aes(xmin=marker_num-0.5, xmax=marker_num+0.5, ymin=cell_num-0.5, ymax=cell_num+0.5), size=2, fill=NA, color="green3") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + 
  scale_x_discrete(labels=c(markers)) + 
  scale_y_discrete(labels=c(cells)) + 
  ggsave("./all_corr_rho.jpg", height=8, width=8)

ggplot(results_spin, aes(x=cell, y=marker, fill=rho_sig)) + 
  geom_tile(color = "white", lwd = 0.1, linetype = 1) + 
  scale_fill_gradientn(colours=c("#0000FFFF", "#FFFFFFFF", "#FF0000FF"), limits=c(-0.6, 0.6)) + 
  geom_text(aes(label=round(as.numeric(rho), 2)), color='white', size=5) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=20),
        axis.text.y = element_text(size=15),
        axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  ggsave("./all_corr_rho_sig.jpg", height=6, width=8)

ggplot(results_spin, aes(x=cell, y=marker, fill=rho_fdr_sig)) + 
  geom_tile(color = "white", lwd = 0.1, linetype = 1) + 
  scale_fill_gradientn(colours=c("#0000FFFF", "#FFFFFFFF", "#FF0000FF"), limits=c(-0.6, 0.6)) + 
  geom_text(aes(label=round(as.numeric(rho), 2)), color='white', size=5) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=20),
        axis.text.y = element_text(size=15),
        axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  ggsave("./all_corr_rho_sig_fdr.jpg", height=6, width=8)

# Wider size for paper figure
ggplot(results_spin, aes(x=cell, y=marker, fill=rho_fdr_sig)) + 
  geom_tile(color = "white", lwd = 0.1, linetype = 1) + 
  scale_fill_gradientn(colours=c("#0000FFFF", "#FFFFFFFF", "#FF0000FF"), limits=c(-0.6, 0.6)) + 
  geom_text(aes(label=round(as.numeric(rho), 2)), color='white', size=6) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=20),
        axis.text.y = element_text(size=15),
        axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  ggsave("./all_corr_rho_sig_fdr_large.jpg", height=6, width=11)

# Correlations AHBA cell-type data between themselves

cells = c("astro", "endo", "micro", "neuroex", "neuroin", "oligo", "opc")

results_spin = as.data.frame(matrix(ncol=4))
colnames(results_spin) = c("cell1", "cell2", "pval", "rho")

for (c1 in 1:(length(cells)-1)){
  for (c2 in 1:(length(cells)-c1)){
    cat(paste0(c1, " ", cells[c1], "\n", c2+c1, " ", cells[c2+c1], "\n\n"))
    pval = read.table(paste0("./spin_test/means_AHBA_", cells[c1],"_", cells[c2+c1],"_pval.txt"))
    rho = read.table(paste0("./spin_test/means_AHBA_", cells[c1],"_", cells[c2+c1],"_rho.txt"))
    
    results_spin[nrow(results_spin)+1,] = c(cells[c1], cells[c2+c1], pval[1,1], rho[1,1])
  }
}

results_spin = results_spin[-1,]
results_spin$pval_fdr = p.adjust(results_spin$pval, method="fdr")

results_spin2 = as.data.frame(results_spin$cell2)
results_spin2$cell2 = results_spin$cell1
results_spin2$pval = results_spin$pval
results_spin2$rho = results_spin$rho
results_spin2$pval_fdr = results_spin$pval_fdr
names(results_spin2) = c("cell1", "cell2", "pval", "rho", "pval_fdr")
results_spin = rbind(results_spin, results_spin2)

# correlation matrix

ggplot(data = results_spin, aes(x=cell1, y=cell2, fill=as.numeric(rho))) + 
  geom_tile(color = "white", size = 0.1) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit=c(-1,1), space="Lab",
                       name="Pearson's r") + 
  theme_minimal() +
  coord_fixed() + 
  geom_text(aes(cell1, cell2, label=paste("r=", round(as.numeric(rho),2), "\np=", round(as.numeric(pval_fdr),2))), color="black", size=2.5) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=10, angle=45, vjust = 0.6),
        axis.text.y = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        legend.text = element_text(size=10)) +
  ggsave("./cells_corr_rho_fdr.jpg")






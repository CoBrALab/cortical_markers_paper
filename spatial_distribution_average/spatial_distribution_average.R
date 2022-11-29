# Across-subjects average of each marker

library(RMINC)
library(readxl)
library(reshape2)

wd = getwd() # IMPORTANT: working directory has to be the same directory that the R file is in
dir.create("./results/", showWarnings = FALSE) # create results folder

# Specify the paths to excel files and vertex files
path_to_csv = "../../master_anon.csv"
path_to_outputs = "../../vertex_files_20mm_anon/"

# Import datasets
all_dataset = read.csv(path_to_csv)

# Remove IDs that don't pass CIVET QC, T1w QC or T2w QC
all_subset = subset(all_dataset, (all_dataset$qc_civet != 0) & (all_dataset$qc_t1w <= 2) & (all_dataset$qc_t2w <= 2), select=c(ID, age, sex))

# Encode vertex files
markers_left = list()
markers_right = list()

# Take inputs that are NOT residualized for curvature
markers_left[[1]] = paste0(path_to_outputs,"BSC/",all_subset$ID,"_model_c_left_log_20mm_rsl.txt")
markers_right[[1]] = paste0(path_to_outputs,"BSC/", all_subset$ID, "_model_c_right_log_20mm_rsl.txt")
markers_left[[2]] = paste0(path_to_outputs,"GWC/",all_subset$ID,"_model_ratio_left_20mm_rsl.txt")
markers_right[[2]] = paste0(path_to_outputs,"GWC/", all_subset$ID, "_model_ratio_right_20mm_rsl.txt")
markers_left[[3]] = paste0(path_to_outputs,"CT/",all_subset$ID,"_native_rms_rsl_tlaplace_20mm_left.txt")
markers_right[[3]] = paste0(path_to_outputs,"CT/", all_subset$ID, "_native_rms_rsl_tlaplace_20mm_right.txt")
markers_left[[4]] = paste0(path_to_outputs,"t1t2_ratio_25_GM/", all_subset$ID, "_25_GM_rsl_left.txt")
markers_right[[4]] = paste0(path_to_outputs,"t1t2_ratio_25_GM/", all_subset$ID, "_25_GM_rsl_right.txt")
markers_left[[5]] = paste0(path_to_outputs,"t1t2_ratio_25_WM/", all_subset$ID, "_25_WM_rsl_left.txt")
markers_right[[5]] = paste0(path_to_outputs,"t1t2_ratio_25_WM/", all_subset$ID, "_25_WM_rsl_right.txt")
markers_left[[6]] = paste0(path_to_outputs,"t1t2_ratio_50_GM/", all_subset$ID, "_50_GM_rsl_left.txt")
markers_right[[6]] = paste0(path_to_outputs,"t1t2_ratio_50_GM/", all_subset$ID, "_50_GM_rsl_right.txt")
markers_left[[7]] = paste0(path_to_outputs,"sigmoid_resid/", all_subset$ID, "_model_resid_left_20mm_rsl.txt")
markers_right[[7]] = paste0(path_to_outputs,"sigmoid_resid/", all_subset$ID, "_model_resid_right_20mm_rsl.txt")

# Marker list (in the right order)
names = c("BSC", "GWC", "CT", "t1t2_ratio_25", "t1t2_ratio_25_WM", "t1t2_ratio_50", "sig_resid")

# Initialize lists
means_left = list()
means_right = list()

row_names = c("-1 SD", "+1 SD", "-1.5 SD", "+1.5 SD", "-2 SD", "+2 SD")
mean_thresholds = as.data.frame(matrix(0, ncol = length(names), nrow = length(row_names)))
colnames(mean_thresholds) = names
thresholds = c(-1, 1, -1.5, 1.5, -2, 2)

# For each marker
for (i in 1:7){
  print(names[i])
  print(head(markers_left[[i]]))

  # Calculate average
  print("Compute vertex-wise means and SD")
  means_left[[i]] = vertexMean(markers_left[[i]])
  means_right[[i]] = vertexMean(markers_right[[i]])

  # Thresholds for visualization
  print("Thresholds for visualization")
  for (t in 1:length(thresholds)){
    mean_thresholds[t,i] = mean(rbind(means_left[[i]], means_right[[i]])) + thresholds[t] * sd(rbind(means_left[[i]], means_right[[i]]))
  }

  # Write results to csv
  write.csv(mean_thresholds, './results/mean_thresholds_all.csv', row.names = TRUE)
  
  print("Write results to csv")
  write.table(means_left[[i]], paste0('./results/mean_all_', names[i],'_left.csv'), sep=",", col.names = TRUE, row.names = FALSE)
  write.table(means_right[[i]], paste0('./results/mean_all_', names[i],'_right.csv'), sep=",", col.names = TRUE, row.names = FALSE)
}
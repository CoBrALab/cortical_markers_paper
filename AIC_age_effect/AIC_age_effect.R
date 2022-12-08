# Akaike transformation criterion (AIC) to evaluate best fit between linear, quadratic and cubic age effects

library(RMINC)
library(ggplot2)
library(readxl)

wd = getwd() # IMPORTANT: working directory has to be the same directory that the R file is in
dir.create("./data_new/", showWarnings = FALSE) # create results folder
dir.create("./spin_test/", showWarnings = FALSE) # create new results folder

# Specify the paths to excel files and vertex files
path_to_csv = "../master_anon.csv"
path_to_outputs = "../vertex_files_20mm_anon/"

# Import dataset
all_dataset = read.csv(path_to_csv)

# Remove IDs that don't pass CIVET QC (0=fail, 2=perfect), T1w QC or T2w QC (1,2 = pass; >2.5 = fail)
all_subset = subset(all_dataset, (all_dataset$qc_civet != 0) & (all_dataset$qc_t1w <= 2) & (all_dataset$qc_t2w <= 2), select=c(ID, age, sex))

# Substract age by the minimum age
all_subset$age = all_subset$age - min(all_subset$age)

# Encode vertex files
markers_left = list()
markers_right = list()

# Take inputs that are residualized for curvature
markers_left[[1]] = paste0(path_to_outputs,"BSC_res_curv/",all_subset$ID,"_BSC_res_curv_left.txt")
markers_right[[1]] = paste0(path_to_outputs,"BSC_res_curv/", all_subset$ID, "_BSC_res_curv_right.txt")
markers_left[[2]] = paste0(path_to_outputs,"GWC_res_curv/",all_subset$ID,"_GWC_res_curv_left.txt")
markers_right[[2]] = paste0(path_to_outputs,"GWC_res_curv/", all_subset$ID, "_GWC_res_curv_right.txt")
markers_left[[3]] = paste0(path_to_outputs,"CT_res_curv/",all_subset$ID,"_CT_res_curv_left.txt")
markers_right[[3]] = paste0(path_to_outputs,"CT_res_curv/", all_subset$ID, "_CT_res_curv_right.txt")
markers_left[[4]] = paste0(path_to_outputs,"t1t2_ratio_25_GM_res_curv/",all_subset$ID,"_t1t2_ratio_25_res_curv_left.txt")
markers_right[[4]] = paste0(path_to_outputs,"t1t2_ratio_25_GM_res_curv/", all_subset$ID, "_t1t2_ratio_25_res_curv_right.txt")
markers_left[[5]] = paste0(path_to_outputs,"t1t2_ratio_25_WM_res_curv/", all_subset$ID, "_t1t2_ratio_25_WM_res_curv_left.txt")
markers_right[[5]] = paste0(path_to_outputs,"t1t2_ratio_25_WM_res_curv/", all_subset$ID, "_t1t2_ratio_25_WM_res_curv_right.txt")
markers_left[[6]] = paste0(path_to_outputs,"t1t2_ratio_50_GM_res_curv/", all_subset$ID, "_t1t2_ratio_50_res_curv_left.txt")
markers_right[[6]] = paste0(path_to_outputs,"t1t2_ratio_50_GM_res_curv/", all_subset$ID, "_t1t2_ratio_50_res_curv_right.txt")

# Marker list (in the right order)
names = c("BSC", "GWC", "CT", "t1t2_ratio_25", "t1t2_ratio_25_WM", "t1t2_ratio_50")

# Initialize lists
lm_age1_left = list()
lm_age1_right = list()
lm_age2_left = list()
lm_age2_right = list()
lm_age3_left = list()
lm_age3_right = list()

results_left = list()
results_right = list()

# Run AIC for the left hemisphere values
for (i in 1:length(names)){
  cat(paste0("\n--------------------------\n", names[i], "\n--------------------------\n"))
  cat(head(markers_left[[i]], 1))
  
  # Fit linear, quadratic and cubic linear models
  cat("\nFit linear models")
  lm_age1_left[[i]] = vertexLm(markers_left[[i]] ~ poly(age,1) + sex, all_subset)
  lm_age2_left[[i]] = vertexLm(markers_left[[i]] ~ poly(age,2) + sex, all_subset)
  lm_age3_left[[i]] = vertexLm(markers_left[[i]] ~ poly(age,3) + sex, all_subset)
  
  # Store results
  cat("\nStore results")
  results_left[[i]] = data.frame(matrix(ncol = 5, nrow = 40962))
  colnames(results_left[[i]]) = c("age1", "age2", "age3", "minAIC", "AIC")
  
  # Run AIC
  cat("\nRun AIC")
  for (v in 1:40962){
    results_left[[i]][,1] = AIC(lm_age1_left[[i]])
    results_left[[i]][,2] = AIC(lm_age2_left[[i]])
    results_left[[i]][,3] = AIC(lm_age3_left[[i]])
    
    results_left[[i]][v,4] = min(results_left[[i]][v,1:3])
    
    if (results_left[[i]][v,1] == results_left[[i]][v,4])  {results_left[[i]][v,5] = as.numeric(1)}
    if (results_left[[i]][v,2] == results_left[[i]][v,4])  {results_left[[i]][v,5] = as.numeric(2)}
    if (results_left[[i]][v,3] == results_left[[i]][v,4])  {results_left[[i]][v,5] = as.numeric(3)}
  }
  
  results_left[[i]]$AIC = as.numeric(results_left[[i]]$AIC)
  
  # Write results to csv
  cat("\nWrite results to csv")
  write.csv(results_left[[i]]$AIC, paste0('./data_new/AIC_', names[i],'_left.csv'), row.names = FALSE)
}

# Run AIC for the right hemisphere values
for (i in 1:length(names)){
  cat(paste0("\n--------------------------\n", names[i], "\n--------------------------\n"))
  cat(head(markers_right[[i]], 1))
  
  # Fit linear, quadratic and cubic linear models
  cat("\nFit linear models")
  lm_age1_right[[i]] = vertexLm(markers_right[[i]] ~ poly(age,1) + sex, all_subset)
  lm_age2_right[[i]] = vertexLm(markers_right[[i]] ~ poly(age,2) + sex, all_subset)
  lm_age3_right[[i]] = vertexLm(markers_right[[i]] ~ poly(age,3) + sex, all_subset)
  
  # Store results
  cat("\nStore results")
  results_right[[i]] = data.frame(matrix(ncol = 5, nrow = 40962))
  colnames(results_right[[i]]) = c("age1", "age2", "age3", "minAIC", "AIC")
  
  # Run AIC
  cat("\nRun AIC")
  for (v in 1:40962){
    results_right[[i]][,1] = AIC(lm_age1_right[[i]])
    results_right[[i]][,2] = AIC(lm_age2_right[[i]])
    results_right[[i]][,3] = AIC(lm_age3_right[[i]])
    
    results_right[[i]][v,4] = min(results_right[[i]][v,1:3])
    
    if (results_right[[i]][v,1] == results_right[[i]][v,4])  {results_right[[i]][v,5] = as.numeric(1)}
    if (results_right[[i]][v,2] == results_right[[i]][v,4])  {results_right[[i]][v,5] = as.numeric(2)}
    if (results_right[[i]][v,3] == results_right[[i]][v,4])  {results_right[[i]][v,5] = as.numeric(3)}
  }
  
  results_right[[i]]$AIC = as.numeric(results_right[[i]]$AIC)
  
  # Write results to csv
  cat("\nWrite results to csv")
  write.csv(results_right[[i]]$AIC, paste0('./data_new/AIC_', names[i],'_right.csv'), row.names = FALSE)
}

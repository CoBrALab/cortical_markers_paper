# BigBrain analysis

library(RMINC)
library(ggplot2)
library(readxl)

wd = getwd()

# Import datasets
#LAM_dataset = read_excel("./LAMQCdatabase_June2020.xlsx", 
#                         sheet = "Scan Info", range = "A1:X185")
ADB_dataset = read_excel("../ADB_master.xlsx", 
                         sheet = "info", range = "A6:AZ74")
HA_dataset = read_excel("../HA_master.xlsx", 
                        sheet = "info", range = "A6:Z112")



# List of IDs used for each measure/dataset

ADB_subset = subset(ADB_dataset, (ADB_dataset$civet_qc_olivier != 0) & (ADB_dataset$qc_aurelie_t1_max <= 2) & (ADB_dataset$qc_aurelie_t2_max <= 2), select=c(ID, age_abs, sex_spelled, T2_sequence_type))
#LAM_subset = subset(LAM_dataset, (LAM_dataset$CIVET_QC_Olivier!=0) & (LAM_dataset$Pre_post_MRI==1) & (LAM_dataset$Used_for_analysis==1) & (LAM_dataset$T1w_Raw_QC != 0), select=c(ID, age_abs, sex_spelled, T2_sequence_type))
HA_subset = subset(HA_dataset, (HA_dataset$civet_qc_olivier != 0) & (HA_dataset$aurelie_QC_t1 <= 2) & (HA_dataset$aurelie_QC_T2 <= 2), select=c(ID, age_abs, sex_spelled, T2_sequence_type))

all_subset = rbind(ADB_subset, HA_subset)

#Spin test files results
all_subset$BSC_rho_file = paste0("./spin_test/single_subjects_BSC/", all_subset$ID, "_BSC_BB_BSC_rho.txt")
all_subset$BSC_pval_file = paste0("./spin_test/single_subjects_BSC/", all_subset$ID, "_BSC_BB_BSC_pval.txt")
all_subset$GWC_rho_file = paste0("./spin_test/single_subjects_GWC/", all_subset$ID, "_GWC_BB_GWC_rho.txt")
all_subset$GWC_pval_file = paste0("./spin_test/single_subjects_GWC/", all_subset$ID, "_GWC_BB_GWC_pval.txt")
all_subset$CT_rho_file = paste0("./spin_test/single_subjects_CT/", all_subset$ID, "_CT_BB_CT_rho.txt")
all_subset$CT_pval_file = paste0("./spin_test/single_subjects_CT/", all_subset$ID, "_CT_BB_CT_pval.txt")
all_subset$t1t2_ratio_25_rho_file = paste0("./spin_test/single_subjects_25/", all_subset$ID, "_t1t2_ratio_25_BB_int_25_rho.txt")
all_subset$t1t2_ratio_25_pval_file = paste0("./spin_test/single_subjects_25/", all_subset$ID, "_t1t2_ratio_25_BB_int_25_pval.txt")
all_subset$t1t2_ratio_25_WM_rho_file = paste0("./spin_test/single_subjects_25_WM/", all_subset$ID, "_t1t2_ratio_25_WM_BB_int_25_WM_rho.txt")
all_subset$t1t2_ratio_25_WM_pval_file = paste0("./spin_test/single_subjects_25_WM/", all_subset$ID, "_t1t2_ratio_25_WM_BB_int_25_WM_pval.txt")

all_subset$BSC_rho = 0
all_subset$BSC_pval = 0
all_subset$GWC_rho = 0
all_subset$GWC_pval = 0
all_subset$CT_rho = 0
all_subset$CT_pval = 0
all_subset$t1t2_ratio_25_rho = 0
all_subset$t1t2_ratio_25_pval = 0
all_subset$t1t2_ratio_25_WM_rho = 0
all_subset$t1t2_ratio_25_WM_pval = 0

for (i in 1:nrow(all_subset)){
  all_subset$BSC_rho[i] = as.numeric(read.table(all_subset$BSC_rho_file[i])[1,1])
  all_subset$BSC_pval[i] = as.numeric(read.table(all_subset$BSC_pval_file[i])[1,1])
  all_subset$GWC_rho[i] = as.numeric(read.table(all_subset$GWC_rho_file[i])[1,1])
  all_subset$GWC_pval[i] = as.numeric(read.table(all_subset$GWC_pval_file[i])[1,1])
  all_subset$CT_rho[i] = as.numeric(read.table(all_subset$CT_rho_file[i])[1,1])
  all_subset$CT_pval[i] = as.numeric(read.table(all_subset$CT_pval_file[i])[1,1])
  all_subset$t1t2_ratio_25_rho[i] = as.numeric(read.table(all_subset$t1t2_ratio_25_rho_file[i])[1,1])
  all_subset$t1t2_ratio_25_pval[i] = as.numeric(read.table(all_subset$t1t2_ratio_25_pval_file[i])[1,1])
  all_subset$t1t2_ratio_25_WM_rho[i] = as.numeric(read.table(all_subset$t1t2_ratio_25_WM_rho_file[i])[1,1])
  all_subset$t1t2_ratio_25_WM_pval[i] = as.numeric(read.table(all_subset$t1t2_ratio_25_WM_pval_file[i])[1,1])
}

summary(all_subset)

##### BSC
#Histogram of rho values
ggplot(all_subset, aes(x = BSC_rho)) + 
  geom_histogram(position = "dodge", binwidth = 0.01) + 
  labs() + xlab("rho") + ylab("#") +
  ggtitle("Spin tests : BigBrain BSC") +
  theme(text=element_text(size=20)) +
  ggsave("./BSC_rho_histogram.png")

#Histogram of p values
ggplot(all_subset, aes(x = BSC_pval)) + 
  geom_histogram(position = "dodge", binwidth = 0.01) + 
  labs() + xlab("pval") + ylab("#") +
  ggtitle("Spin tests : BigBrain BSC") +
  theme(text=element_text(size=20)) +
  ggsave("./BSC_pval_histogram.png")

##### GWC
#Histogram of rho values
ggplot(all_subset, aes(x = GWC_rho)) + 
  geom_histogram(position = "dodge", binwidth = 0.01) + 
  labs() + xlab("rho") + ylab("#") +
  ggtitle("Spin tests : BigBrain GWC") +
  theme(text=element_text(size=20)) +
  ggsave("./GWC_rho_histogram.png")

#Histogram of p values
ggplot(all_subset, aes(x = GWC_pval)) + 
  geom_histogram(position = "dodge", binwidth = 0.01) + 
  labs() + xlab("pval") + ylab("#") +
  ggtitle("Spin tests : BigBrain GWC") +
  theme(text=element_text(size=20)) +
  ggsave("./GWC_pval_histogram.png")

##### CT
#Histogram of rho values
ggplot(all_subset, aes(x = CT_rho)) + 
  geom_histogram(position = "dodge", binwidth = 0.01) + 
  labs() + xlab("rho") + ylab("#") +
  ggtitle("Spin tests : BigBrain CT") +
  theme(text=element_text(size=20)) +
  ggsave("./CT_rho_histogram.png")

#Histogram of p values
ggplot(all_subset, aes(x = CT_pval)) + 
  geom_histogram(position = "dodge", binwidth = 0.001) + 
  labs() + xlab("pval") + ylab("#") +
  ggtitle("Spin tests : BigBrain CT") +
  theme(text=element_text(size=20)) +
  ggsave("./CT_pval_histogram.png")

##### t1t2_ratio_25
#Histogram of rho values
ggplot(all_subset, aes(x = t1t2_ratio_25_rho)) + 
  geom_histogram(position = "dodge", binwidth = 0.01) + 
  labs() + xlab("rho") + ylab("#") +
  ggtitle("Spin tests : BigBrain t1t2_ratio_25") +
  theme(text=element_text(size=20)) +
  ggsave("./t1t2_ratio_25_rho_histogram.png")

#Histogram of p values
ggplot(all_subset, aes(x = t1t2_ratio_25_pval)) + 
  geom_histogram(position = "dodge", binwidth = 0.01) + 
  labs() + xlab("pval") + ylab("#") +
  ggtitle("Spin tests : BigBrain t1t2_ratio_25") +
  theme(text=element_text(size=20)) +
  ggsave("./t1t2_ratio_25_pval_histogram.png")

##### t1t2_ratio_25_WM
#Histogram of rho values
ggplot(all_subset, aes(x = t1t2_ratio_25_WM_rho)) + 
  geom_histogram(position = "dodge", binwidth = 0.01) + 
  labs() + xlab("rho") + ylab("#") +
  ggtitle("Spin tests : BigBrain t1t2_ratio_25_WM") +
  theme(text=element_text(size=20)) +
  ggsave("./t1t2_ratio_25_WM_rho_histogram.png")

#Histogram of p values
ggplot(all_subset, aes(x = t1t2_ratio_25_WM_pval)) + 
  geom_histogram(position = "dodge", binwidth = 0.01) + 
  labs() + xlab("pval") + ylab("#") +
  ggtitle("Spin tests : BigBrain t1t2_ratio_25_WM") +
  theme(text=element_text(size=20)) +
  ggsave("./t1t2_ratio_25_WM_pval_histogram.png")

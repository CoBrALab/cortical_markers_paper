
analysis_directory = '../../R1_descriptive_stats/spin_test/';

%%%%%%%%%%%%%%%%% Means

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spin BSC data
readleft = strcat(analysis_directory, 'mean_all_BSC_left_spin_test.csv');
readright = strcat(analysis_directory, 'mean_all_BSC_right_spin_test.csv');
permno = 500; % how many spins
wsname = sprintf(strcat(analysis_directory, 'means_BSC_rot.mat'));
%SpinPermuCIVET(readleft,readright,permno,wsname)
 
% % Compare BSC vs GWC
% readleft1 = readleft;
% readright1 = readright;
% readleft2 = strcat(analysis_directory, 'mean_all_GWC_left_spin_test.csv');
% readright2 = strcat(analysis_directory, 'mean_all_GWC_right_spin_test.csv');
% real_rho_path = strcat(analysis_directory,'means_BSC_GWC_rho.txt');
% resid_path_left = strcat(analysis_directory,'means_BSC_GWC_resid_left.txt');
% resid_path_right = strcat(analysis_directory,'means_BSC_GWC_resid_right.txt');
% pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
% dlmwrite(strcat(analysis_directory,'means_BSC_GWC_pval.txt'), pval);
% 
% % Compare BSC vs CT
% readleft1 = readleft;
% readright1 = readright;
% readleft2 = strcat(analysis_directory, 'mean_all_CT_left_spin_test.csv');
% readright2 = strcat(analysis_directory, 'mean_all_CT_right_spin_test.csv');
% real_rho_path = strcat(analysis_directory,'means_BSC_CT_rho.txt');
% resid_path_left = strcat(analysis_directory,'means_BSC_CT_resid_left.txt');
% resid_path_right = strcat(analysis_directory,'means_BSC_CT_resid_right.txt');
% pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
% dlmwrite(strcat(analysis_directory,'means_BSC_CT_pval.txt'), pval);
% 
% % Compare BSC vs t1t2_ratio_25
% readleft1 = readleft;
% readright1 = readright;
% readleft2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_left_spin_test.csv');
% readright2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_right_spin_test.csv');
% real_rho_path = strcat(analysis_directory,'means_BSC_t1t2_ratio_25_rho.txt');
% resid_path_left = strcat(analysis_directory,'means_BSC_t1t2_ratio_25_resid_left.txt');
% resid_path_right = strcat(analysis_directory,'means_BSC_t1t2_ratio_25_resid_right.txt');
% pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
% dlmwrite(strcat(analysis_directory,'means_BSC_t1t2_ratio_25_pval.txt'), pval);
% 
% % Compare BSC vs t1t2_ratio_25_WM
% readleft1 = readleft;
% readright1 = readright;
% readleft2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_WM_left_spin_test.csv');
% readright2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_WM_right_spin_test.csv');
% real_rho_path = strcat(analysis_directory,'means_BSC_t1t2_ratio_25_WM_rho.txt');
% resid_path_left = strcat(analysis_directory,'means_BSC_t1t2_ratio_25_WM_resid_left.txt');
% resid_path_right = strcat(analysis_directory,'means_BSC_t1t2_ratio_25_WM_resid_right.txt');
% pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
% dlmwrite(strcat(analysis_directory,'means_BSC_t1t2_ratio_25_WM_pval.txt'), pval);
% 
% % Compare BSC vs t1t2_ratio_50
% readleft1 = readleft;
% readright1 = readright;
% readleft2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_50_left_spin_test.csv');
% readright2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_50_right_spin_test.csv');
% real_rho_path = strcat(analysis_directory,'means_BSC_t1t2_ratio_50_rho.txt');
% resid_path_left = strcat(analysis_directory,'means_BSC_t1t2_ratio_50_resid_left.txt');
% resid_path_right = strcat(analysis_directory,'means_BSC_t1t2_ratio_50_resid_right.txt');
% pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
% dlmwrite(strcat(analysis_directory,'means_BSC_t1t2_ratio_50_pval.txt'), pval);

% Compare BSC vs R1_25
readleft1 = readleft;
readright1 = readright;
x_label = "BSC";
y_label = "GM R1";
readleft2 = strcat(analysis_directory, 'mean_all_R1_25_left_spin_test.csv');
readright2 = strcat(analysis_directory, 'mean_all_R1_25_right_spin_test.csv');
real_rho_path = strcat(analysis_directory,'means_BSC_R1_25_rho.txt');
resid_path_left = strcat(analysis_directory,'means_BSC_R1_25_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_BSC_R1_25_resid_right.txt');
pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_BSC_R1_25_pval.txt'), pval);

% Compare BSC vs R1_25_WM
readleft1 = readleft;
readright1 = readright;
x_label = "BSC";
y_label = "SWM R1";
readleft2 = strcat(analysis_directory, 'mean_all_R1_25_WM_left_spin_test.csv');
readright2 = strcat(analysis_directory, 'mean_all_R1_25_WM_right_spin_test.csv');
real_rho_path = strcat(analysis_directory,'means_BSC_R1_25_WM_rho.txt');
resid_path_left = strcat(analysis_directory,'means_BSC_R1_25_WM_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_BSC_R1_25_WM_resid_right.txt');
pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_BSC_R1_25_WM_pval.txt'), pval);

% Compare BSC vs R1_50
readleft1 = readleft;
readright1 = readright;
x_label = "BSC";
y_label = "GM R1 (50%)";
readleft2 = strcat(analysis_directory, 'mean_all_R1_50_left_spin_test.csv');
readright2 = strcat(analysis_directory, 'mean_all_R1_50_right_spin_test.csv');
real_rho_path = strcat(analysis_directory,'means_BSC_R1_50_rho.txt');
resid_path_left = strcat(analysis_directory,'means_BSC_R1_50_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_BSC_R1_50_resid_right.txt');
pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_BSC_R1_50_pval.txt'), pval);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spin GWC data
readleft = strcat(analysis_directory, 'mean_all_GWC_left_spin_test.csv');
readright = strcat(analysis_directory, 'mean_all_GWC_right_spin_test.csv');
permno = 500; % how many spins
wsname = sprintf(strcat(analysis_directory, 'means_GWC_rot.mat'));
%SpinPermuCIVET(readleft,readright,permno,wsname)

% % Compare GWC vs CT
% readleft1 = readleft;
% readright1 = readright;
% readleft2 = strcat(analysis_directory, 'mean_all_CT_left_spin_test.csv');
% readright2 = strcat(analysis_directory, 'mean_all_CT_right_spin_test.csv');
% real_rho_path = strcat(analysis_directory,'means_GWC_CT_rho.txt');
% resid_path_left = strcat(analysis_directory,'means_GWC_CT_resid_left.txt');
% resid_path_right = strcat(analysis_directory,'means_GWC_CT_resid_right.txt');
% pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
% dlmwrite(strcat(analysis_directory,'means_GWC_CT_pval.txt'), pval);
% 
% % Compare GWC vs t1t2_ratio_25
% readleft1 = readleft;
% readright1 = readright;
% readleft2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_left_spin_test.csv');
% readright2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_right_spin_test.csv');
% real_rho_path = strcat(analysis_directory,'means_GWC_t1t2_ratio_25_rho.txt');
% resid_path_left = strcat(analysis_directory,'means_GWC_t1t2_ratio_25_resid_left.txt');
% resid_path_right = strcat(analysis_directory,'means_GWC_t1t2_ratio_25_resid_right.txt');
% pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
% dlmwrite(strcat(analysis_directory,'means_GWC_t1t2_ratio_25_pval.txt'), pval);
% 
% % Compare GWC vs t1t2_ratio_25_WM
% readleft1 = readleft;
% readright1 = readright;
% readleft2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_WM_left_spin_test.csv');
% readright2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_WM_right_spin_test.csv');
% real_rho_path = strcat(analysis_directory,'means_GWC_t1t2_ratio_25_WM_rho.txt');
% resid_path_left = strcat(analysis_directory,'means_GWC_t1t2_ratio_25_WM_resid_left.txt');
% resid_path_right = strcat(analysis_directory,'means_GWC_t1t2_ratio_25_WM_resid_right.txt');
% pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
% dlmwrite(strcat(analysis_directory,'means_GWC_t1t2_ratio_25_WM_pval.txt'), pval);
% 
% % Compare GWC vs t1t2_ratio_50
% readleft1 = readleft;
% readright1 = readright;
% readleft2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_50_left_spin_test.csv');
% readright2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_50_right_spin_test.csv');
% real_rho_path = strcat(analysis_directory,'means_GWC_t1t2_ratio_50_rho.txt');
% resid_path_left = strcat(analysis_directory,'means_GWC_t1t2_ratio_50_resid_left.txt');
% resid_path_right = strcat(analysis_directory,'means_GWC_t1t2_ratio_50_resid_right.txt');
% pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
% dlmwrite(strcat(analysis_directory,'means_GWC_t1t2_ratio_50_pval.txt'), pval);

% Compare GWC vs R1_25
readleft1 = readleft;
readright1 = readright;
x_label = "GWC";
y_label = "GM R1";
readleft2 = strcat(analysis_directory, 'mean_all_R1_25_left_spin_test.csv');
readright2 = strcat(analysis_directory, 'mean_all_R1_25_right_spin_test.csv');
real_rho_path = strcat(analysis_directory,'means_GWC_R1_25_rho.txt');
resid_path_left = strcat(analysis_directory,'means_GWC_R1_25_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_GWC_R1_25_resid_right.txt');
pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_GWC_R1_25_pval.txt'), pval);

% Compare GWC vs R1_25_WM
readleft1 = readleft;
readright1 = readright;
x_label = "GWC";
y_label = "SWM R1";
readleft2 = strcat(analysis_directory, 'mean_all_R1_25_WM_left_spin_test.csv');
readright2 = strcat(analysis_directory, 'mean_all_R1_25_WM_right_spin_test.csv');
real_rho_path = strcat(analysis_directory,'means_GWC_R1_25_WM_rho.txt');
resid_path_left = strcat(analysis_directory,'means_GWC_R1_25_WM_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_GWC_R1_25_WM_resid_right.txt');
pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_GWC_R1_25_WM_pval.txt'), pval);

% Compare GWC vs R1_50
readleft1 = readleft;
readright1 = readright;
x_label = "GWC";
y_label = "GM R1 (50%)";
readleft2 = strcat(analysis_directory, 'mean_all_R1_50_left_spin_test.csv');
readright2 = strcat(analysis_directory, 'mean_all_R1_50_right_spin_test.csv');
real_rho_path = strcat(analysis_directory,'means_GWC_R1_50_rho.txt');
resid_path_left = strcat(analysis_directory,'means_GWC_R1_50_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_GWC_R1_50_resid_right.txt');
pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_GWC_R1_50_pval.txt'), pval);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spin CT data
readleft = strcat(analysis_directory, 'mean_all_CT_left_spin_test.csv');
readright = strcat(analysis_directory, 'mean_all_CT_right_spin_test.csv');
permno = 500; % how many spins
wsname = sprintf(strcat(analysis_directory, 'means_CT_rot.mat'));
%SpinPermuCIVET(readleft,readright,permno,wsname)

% % Compare CT vs t1t2_ratio_25
% readleft1 = readleft;
% readright1 = readright;
% readleft2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_left_spin_test.csv');
% readright2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_right_spin_test.csv');
% real_rho_path = strcat(analysis_directory,'means_CT_t1t2_ratio_25_rho.txt');
% resid_path_left = strcat(analysis_directory,'means_CT_t1t2_ratio_25_resid_left.txt');
% resid_path_right = strcat(analysis_directory,'means_CT_t1t2_ratio_25_resid_right.txt');
% pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
% dlmwrite(strcat(analysis_directory,'means_CT_t1t2_ratio_25_pval.txt'), pval);
% 
% % Compare CT vs t1t2_ratio_25_WM
% readleft1 = readleft;
% readright1 = readright;
% readleft2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_WM_left_spin_test.csv');
% readright2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_WM_right_spin_test.csv');
% real_rho_path = strcat(analysis_directory,'means_CT_t1t2_ratio_25_WM_rho.txt');
% resid_path_left = strcat(analysis_directory,'means_CT_t1t2_ratio_25_WM_resid_left.txt');
% resid_path_right = strcat(analysis_directory,'means_CT_t1t2_ratio_25_WM_resid_right.txt');
% pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
% dlmwrite(strcat(analysis_directory,'means_CT_t1t2_ratio_25_WM_pval.txt'), pval);
% 
% % Compare CT vs t1t2_ratio_50
% readleft1 = readleft;
% readright1 = readright;
% readleft2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_50_left_spin_test.csv');
% readright2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_50_right_spin_test.csv');
% real_rho_path = strcat(analysis_directory,'means_CT_t1t2_ratio_50_rho.txt');
% resid_path_left = strcat(analysis_directory,'means_CT_t1t2_ratio_50_resid_left.txt');
% resid_path_right = strcat(analysis_directory,'means_CT_t1t2_ratio_50_resid_right.txt');
% pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
% dlmwrite(strcat(analysis_directory,'means_CT_t1t2_ratio_50_pval.txt'), pval);

% Compare CT vs R1_25
readleft1 = readleft;
readright1 = readright;
x_label = "CT";
y_label = "GM R1";
readleft2 = strcat(analysis_directory, 'mean_all_R1_25_left_spin_test.csv');
readright2 = strcat(analysis_directory, 'mean_all_R1_25_right_spin_test.csv');
real_rho_path = strcat(analysis_directory,'means_CT_R1_25_rho.txt');
resid_path_left = strcat(analysis_directory,'means_CT_R1_25_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_CT_R1_25_resid_right.txt');
pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_CT_R1_25_pval.txt'), pval);

% Compare CT vs R1_25_WM
readleft1 = readleft;
readright1 = readright;
x_label = "CT";
y_label = "SWM R1";
readleft2 = strcat(analysis_directory, 'mean_all_R1_25_WM_left_spin_test.csv');
readright2 = strcat(analysis_directory, 'mean_all_R1_25_WM_right_spin_test.csv');
real_rho_path = strcat(analysis_directory,'means_CT_R1_25_WM_rho.txt');
resid_path_left = strcat(analysis_directory,'means_CT_R1_25_WM_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_CT_R1_25_WM_resid_right.txt');
pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_CT_R1_25_WM_pval.txt'), pval);

% Compare CT vs R1_50
readleft1 = readleft;
readright1 = readright;
x_label = "CT";
y_label = "GM R1 (50%)";
readleft2 = strcat(analysis_directory, 'mean_all_R1_50_left_spin_test.csv');
readright2 = strcat(analysis_directory, 'mean_all_R1_50_right_spin_test.csv');
real_rho_path = strcat(analysis_directory,'means_CT_R1_50_rho.txt');
resid_path_left = strcat(analysis_directory,'means_CT_R1_50_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_CT_R1_50_resid_right.txt');
pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_CT_R1_50_pval.txt'), pval);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spin t1t2_ratio_25 data
readleft = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_left_spin_test.csv');
readright = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_right_spin_test.csv');
permno = 500; % how many spins
wsname = sprintf(strcat(analysis_directory, 'means_t1t2_ratio_25_rot.mat'));
%SpinPermuCIVET(readleft,readright,permno,wsname)

% % Compare t1t2_ratio_25 vs t1t2_ratio_25_WM
% readleft1 = readleft;
% readright1 = readright;
% readleft2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_WM_left_spin_test.csv');
% readright2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_WM_right_spin_test.csv');
% real_rho_path = strcat(analysis_directory,'means_t1t2_ratio_25_t1t2_ratio_25_WM_rho.txt');
% resid_path_left = strcat(analysis_directory,'means_t1t2_ratio_25_t1t2_ratio_25_WM_resid_left.txt');
% resid_path_right = strcat(analysis_directory,'means_t1t2_ratio_25_t1t2_ratio_25_WM_resid_right.txt');
% pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
% dlmwrite(strcat(analysis_directory,'means_t1t2_ratio_25_t1t2_ratio_25_WM_pval.txt'), pval);
% 
% % Compare t1t2_ratio_25 vs t1t2_ratio_50
% readleft1 = readleft;
% readright1 = readright;
% readleft2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_50_left_spin_test.csv');
% readright2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_50_right_spin_test.csv');
% real_rho_path = strcat(analysis_directory,'means_t1t2_ratio_25_t1t2_ratio_50_rho.txt');
% resid_path_left = strcat(analysis_directory,'means_t1t2_ratio_25_t1t2_ratio_50_resid_left.txt');
% resid_path_right = strcat(analysis_directory,'means_t1t2_ratio_25_t1t2_ratio_50_resid_right.txt');
% pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
% dlmwrite(strcat(analysis_directory,'means_t1t2_ratio_25_t1t2_ratio_50_pval.txt'), pval);

% Compare t1t2_ratio_25 vs R1_25
readleft1 = readleft;
readright1 = readright;
x_label = "GM T1w/T2w";
y_label = "GM R1";
readleft2 = strcat(analysis_directory, 'mean_all_R1_25_left_spin_test.csv');
readright2 = strcat(analysis_directory, 'mean_all_R1_25_right_spin_test.csv');
real_rho_path = strcat(analysis_directory,'means_t1t2_ratio_25_R1_25_rho.txt');
resid_path_left = strcat(analysis_directory,'means_t1t2_ratio_25_R1_25_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_t1t2_ratio_25_R1_25_resid_right.txt');
pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_t1t2_ratio_25_R1_25_pval.txt'), pval);

% Compare t1t2_ratio_25 vs R1_25_WM
readleft1 = readleft;
readright1 = readright;
x_label = "GM T1w/T2w";
y_label = "SWM R1";
readleft2 = strcat(analysis_directory, 'mean_all_R1_25_WM_left_spin_test.csv');
readright2 = strcat(analysis_directory, 'mean_all_R1_25_WM_right_spin_test.csv');
real_rho_path = strcat(analysis_directory,'means_t1t2_ratio_25_R1_25_WM_rho.txt');
resid_path_left = strcat(analysis_directory,'means_t1t2_ratio_25_R1_25_WM_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_t1t2_ratio_25_R1_25_WM_resid_right.txt');
pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_t1t2_ratio_25_R1_25_WM_pval.txt'), pval);

% Compare t1t2_ratio_25 vs R1_50
readleft1 = readleft;
readright1 = readright;
x_label = "GM T1w/T2w";
y_label = "GM R1 (50%)";
readleft2 = strcat(analysis_directory, 'mean_all_R1_50_left_spin_test.csv');
readright2 = strcat(analysis_directory, 'mean_all_R1_50_right_spin_test.csv');
real_rho_path = strcat(analysis_directory,'means_t1t2_ratio_25_R1_50_rho.txt');
resid_path_left = strcat(analysis_directory,'means_t1t2_ratio_25_R1_50_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_t1t2_ratio_25_R1_50_resid_right.txt');
pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_t1t2_ratio_25_R1_50_pval.txt'), pval);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spin t1t2_ratio_25_WM data
readleft = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_WM_left_spin_test.csv');
readright = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_WM_right_spin_test.csv');
permno = 500; % how many spins
wsname = sprintf(strcat(analysis_directory, 'means_t1t2_ratio_25_WM_rot.mat'));
%SpinPermuCIVET(readleft,readright,permno,wsname)

% % Compare t1t2_ratio_25_WM vs t1t2_ratio_50
% readleft1 = readleft;
% readright1 = readright;
% readleft2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_50_left_spin_test.csv');
% readright2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_50_right_spin_test.csv');
% real_rho_path = strcat(analysis_directory,'means_t1t2_ratio_25_WM_t1t2_ratio_50_rho.txt');
% resid_path_left = strcat(analysis_directory,'means_t1t2_ratio_25_WM_t1t2_ratio_50_resid_left.txt');
% resid_path_right = strcat(analysis_directory,'means_t1t2_ratio_25_WM_t1t2_ratio_50_resid_right.txt');
% pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
% dlmwrite(strcat(analysis_directory,'means_t1t2_ratio_25_WM_t1t2_ratio_50_pval.txt'), pval);

% Compare t1t2_ratio_25_WM vs R1_25
readleft1 = readleft;
readright1 = readright;
x_label = "SWM T1w/T2w";
y_label = "GM R1";
readleft2 = strcat(analysis_directory, 'mean_all_R1_25_left_spin_test.csv');
readright2 = strcat(analysis_directory, 'mean_all_R1_25_right_spin_test.csv');
real_rho_path = strcat(analysis_directory,'means_t1t2_ratio_25_WM_R1_25_rho.txt');
resid_path_left = strcat(analysis_directory,'means_t1t2_ratio_25_WM_R1_25_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_t1t2_ratio_25_WM_R1_25_resid_right.txt');
pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_t1t2_ratio_25_WM_R1_25_pval.txt'), pval);

% Compare t1t2_ratio_25_WM vs R1_25_WM
readleft1 = readleft;
readright1 = readright;
x_label = "SWM T1w/T2w";
y_label = "SWM R1";
readleft2 = strcat(analysis_directory, 'mean_all_R1_25_WM_left_spin_test.csv');
readright2 = strcat(analysis_directory, 'mean_all_R1_25_WM_right_spin_test.csv');
real_rho_path = strcat(analysis_directory,'means_t1t2_ratio_25_WM_R1_25_WM_rho.txt');
resid_path_left = strcat(analysis_directory,'means_t1t2_ratio_25_WM_R1_25_WM_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_t1t2_ratio_25_WM_R1_25_WM_resid_right.txt');
pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_t1t2_ratio_25_WM_R1_25_WM_pval.txt'), pval);

% Compare t1t2_ratio_25_WM vs R1_50
readleft1 = readleft;
readright1 = readright;
x_label = "SWM T1w/T2w";
y_label = "GM R1 (50%)";
readleft2 = strcat(analysis_directory, 'mean_all_R1_50_left_spin_test.csv');
readright2 = strcat(analysis_directory, 'mean_all_R1_50_right_spin_test.csv');
real_rho_path = strcat(analysis_directory,'means_t1t2_ratio_25_WM_R1_50_rho.txt');
resid_path_left = strcat(analysis_directory,'means_t1t2_ratio_25_WM_R1_50_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_t1t2_ratio_25_WM_R1_50_resid_right.txt');
pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_t1t2_ratio_25_WM_R1_50_pval.txt'), pval);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spin t1t2_ratio_50 data
readleft = strcat(analysis_directory, 'mean_all_t1t2_ratio_50_left_spin_test.csv');
readright = strcat(analysis_directory, 'mean_all_t1t2_ratio_50_right_spin_test.csv');
permno = 500; % how many spins
wsname = sprintf(strcat(analysis_directory, 'means_t1t2_ratio_50_rot.mat'));
%SpinPermuCIVET(readleft,readright,permno,wsname)

% Compare t1t2_ratio_50 vs R1_25
readleft1 = readleft;
readright1 = readright;
x_label = "GM T1w/T2w (50%)";
y_label = "GM R1";
readleft2 = strcat(analysis_directory, 'mean_all_R1_25_left_spin_test.csv');
readright2 = strcat(analysis_directory, 'mean_all_R1_25_right_spin_test.csv');
real_rho_path = strcat(analysis_directory,'means_t1t2_ratio_50_R1_25_rho.txt');
resid_path_left = strcat(analysis_directory,'means_t1t2_ratio_50_R1_25_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_t1t2_ratio_50_R1_25_resid_right.txt');
pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_t1t2_ratio_50_R1_25_pval.txt'), pval);

% Compare t1t2_ratio_50 vs R1_25_WM
readleft1 = readleft;
readright1 = readright;
x_label = "GM T1w/T2w (50%)";
y_label = "SWM R1";
readleft2 = strcat(analysis_directory, 'mean_all_R1_25_WM_left_spin_test.csv');
readright2 = strcat(analysis_directory, 'mean_all_R1_25_WM_right_spin_test.csv');
real_rho_path = strcat(analysis_directory,'means_t1t2_ratio_50_R1_25_WM_rho.txt');
resid_path_left = strcat(analysis_directory,'means_t1t2_ratio_50_R1_25_WM_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_t1t2_ratio_50_R1_25_WM_resid_right.txt');
pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_t1t2_ratio_50_R1_25_WM_pval.txt'), pval);

% Compare t1t2_ratio_50 vs R1_50
readleft1 = readleft;
readright1 = readright;
x_label = "GM T1w/T2w (50%)";
y_label = "GM R1 (50%)";
readleft2 = strcat(analysis_directory, 'mean_all_R1_50_left_spin_test.csv');
readright2 = strcat(analysis_directory, 'mean_all_R1_50_right_spin_test.csv');
real_rho_path = strcat(analysis_directory,'means_t1t2_ratio_50_R1_50_rho.txt');
resid_path_left = strcat(analysis_directory,'means_t1t2_ratio_50_R1_50_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_t1t2_ratio_50_R1_50_resid_right.txt');
pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_t1t2_ratio_50_R1_50_pval.txt'), pval);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Spin R1_25 data
% readleft = strcat(analysis_directory, 'mean_all_R1_25_left_spin_test.csv');
% readright = strcat(analysis_directory, 'mean_all_R1_25_right_spin_test.csv');
% permno = 500; % how many spins
% wsname = sprintf(strcat(analysis_directory, 'means_R1_25_rot.mat'));
% SpinPermuCIVET(readleft,readright,permno,wsname)
% 
% % Compare R1_25 vs R1_25_WM
% readleft1 = readleft;
% readright1 = readright;
% readleft2 = strcat(analysis_directory, 'mean_all_R1_25_WM_left_spin_test.csv');
% readright2 = strcat(analysis_directory, 'mean_all_R1_25_WM_right_spin_test.csv');
% real_rho_path = strcat(analysis_directory,'means_R1_25_R1_25_WM_rho.txt');
% resid_path_left = strcat(analysis_directory,'means_R1_25_R1_25_WM_resid_left.txt');
% resid_path_right = strcat(analysis_directory,'means_R1_25_R1_25_WM_resid_right.txt');
% pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
% dlmwrite(strcat(analysis_directory,'means_R1_25_R1_25_WM_pval.txt'), pval);
% 
% % Compare R1_25 vs R1_50
% readleft1 = readleft;
% readright1 = readright;
% readleft2 = strcat(analysis_directory, 'mean_all_R1_50_left_spin_test.csv');
% readright2 = strcat(analysis_directory, 'mean_all_R1_50_right_spin_test.csv');
% real_rho_path = strcat(analysis_directory,'means_R1_25_R1_50_rho.txt');
% resid_path_left = strcat(analysis_directory,'means_R1_25_R1_50_resid_left.txt');
% resid_path_right = strcat(analysis_directory,'means_R1_25_R1_50_resid_right.txt');
% pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
% dlmwrite(strcat(analysis_directory,'means_R1_25_R1_50_pval.txt'), pval);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Spin R1_25_WM data
% readleft = strcat(analysis_directory, 'mean_all_R1_25_WM_left_spin_test.csv');
% readright = strcat(analysis_directory, 'mean_all_R1_25_WM_right_spin_test.csv');
% permno = 500; % how many spins
% wsname = sprintf(strcat(analysis_directory, 'means_R1_25_WM_rot.mat'));
% SpinPermuCIVET(readleft,readright,permno,wsname)
% 
% % Compare R1_25_WM vs R1_50
% readleft1 = readleft;
% readright1 = readright;
% readleft2 = strcat(analysis_directory, 'mean_all_R1_50_left_spin_test.csv');
% readright2 = strcat(analysis_directory, 'mean_all_R1_50_right_spin_test.csv');
% real_rho_path = strcat(analysis_directory,'means_R1_25_WM_R1_50_rho.txt');
% resid_path_left = strcat(analysis_directory,'means_R1_25_WM_R1_50_resid_left.txt');
% resid_path_right = strcat(analysis_directory,'means_R1_25_WM_R1_50_resid_right.txt');
% pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right); 
% dlmwrite(strcat(analysis_directory,'means_R1_25_WM_R1_50_pval.txt'), pval);

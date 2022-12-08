
analysis_directory = '../spatial_distribution_average/data_new/';

%%%%%%%%%%%%%%%%%%%%%%%%%% BSC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spin BSC data
readleft = strcat(analysis_directory, 'mean_all_BSC_left.csv');
readright = strcat(analysis_directory, 'mean_all_BSC_right.csv');
permno = 1000; % how many spins
wsname = sprintf(strcat(analysis_directory, 'means_BSC_rot.mat'));
SpinPermuCIVET(readleft,readright,permno,wsname)

% Compare BSC vs GWC
x_label = 'BSC';
y_label = 'GWC';
readleft1 = readleft;
readright1 = readright;
readleft2 = strcat(analysis_directory, 'mean_all_GWC_left.csv');
readright2 = strcat(analysis_directory, 'mean_all_GWC_right.csv');
real_rho_path = strcat(analysis_directory,'means_BSC_GWC_rho.txt');
resid_path_left = strcat(analysis_directory,'means_BSC_GWC_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_BSC_GWC_resid_right.txt');
pval=pvalvsNull_resid(x_label, y_label, readleft1,readright1,readleft2,readright2,permno,wsname,real_rho_path, resid_path_left, resid_path_right);
dlmwrite(strcat(analysis_directory,'means_BSC_GWC_pval.txt'), pval);

% Compare BSC vs CT
x_label = 'BSC';
y_label = 'CT';
readleft1 = readleft;
readright1 = readright;
readleft2 = strcat(analysis_directory, 'mean_all_CT_left.csv');
readright2 = strcat(analysis_directory, 'mean_all_CT_right.csv');
real_rho_path = strcat(analysis_directory,'means_BSC_CT_rho.txt');
resid_path_left = strcat(analysis_directory,'means_BSC_CT_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_BSC_CT_resid_right.txt');
pval=pvalvsNull_resid(x_label, y_label, readleft1,readright1,readleft2,readright2,permno,wsname,real_rho_path, resid_path_left, resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_BSC_CT_pval.txt'), pval);

% Compare BSC vs t1t2_ratio_25
x_label = 'BSC';
y_label = 'GM T1w/T2w';
readleft1 = readleft;
readright1 = readright;
readleft2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_left.csv');
readright2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_right.csv');
real_rho_path = strcat(analysis_directory,'means_BSC_t1t2_ratio_25_rho.txt');
resid_path_left = strcat(analysis_directory,'means_BSC_t1t2_ratio_25_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_BSC_t1t2_ratio_25_resid_right.txt');
pval=pvalvsNull_resid(x_label, y_label, readleft1,readright1,readleft2,readright2,permno,wsname,real_rho_path, resid_path_left, resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_BSC_t1t2_ratio_25_pval.txt'), pval);

% Compare BSC vs t1t2_ratio_25_WM
x_label = 'BSC';
y_label = 'SWM T1w/T2w';
readleft1 = readleft;
readright1 = readright;
readleft2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_WM_left.csv');
readright2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_WM_right.csv');
real_rho_path = strcat(analysis_directory,'means_BSC_t1t2_ratio_25_WM_rho.txt');
resid_path_left = strcat(analysis_directory,'means_BSC_t1t2_ratio_25_WM_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_BSC_t1t2_ratio_25_WM_resid_right.txt');
pval=pvalvsNull_resid(x_label, y_label, readleft1,readright1,readleft2,readright2,permno,wsname,real_rho_path, resid_path_left, resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_BSC_t1t2_ratio_25_WM_pval.txt'), pval);

% Compare BSC vs t1t2_ratio_50
x_label = 'BSC';
y_label = 'GM T1w/T2w (50%)';
readleft1 = readleft;
readright1 = readright;
readleft2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_50_left.csv');
readright2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_50_right.csv');
real_rho_path = strcat(analysis_directory,'means_BSC_t1t2_ratio_50_rho.txt');
resid_path_left = strcat(analysis_directory,'means_BSC_t1t2_ratio_50_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_BSC_t1t2_ratio_50_resid_right.txt');
pval=pvalvsNull_resid(x_label, y_label, readleft1,readright1,readleft2,readright2,permno,wsname,real_rho_path, resid_path_left, resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_BSC_t1t2_ratio_50_pval.txt'), pval);

%%%%%%%%%%%%%%%%%%%%%%%%%% GWC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spin GWC data
readleft = strcat(analysis_directory, 'mean_all_GWC_left.csv');
readright = strcat(analysis_directory, 'mean_all_GWC_right.csv');
permno = 1000; % how many spins
wsname = sprintf(strcat(analysis_directory, 'means_GWC_rot.mat'));
SpinPermuCIVET(readleft,readright,permno,wsname)

% Compare GWC vs CT
x_label = 'GWC';
y_label = 'CT';
readleft1 = readleft;
readright1 = readright;
readleft2 = strcat(analysis_directory, 'mean_all_CT_left.csv');
readright2 = strcat(analysis_directory, 'mean_all_CT_right.csv');
real_rho_path = strcat(analysis_directory,'means_GWC_CT_rho.txt');
resid_path_left = strcat(analysis_directory,'means_GWC_CT_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_GWC_CT_resid_right.txt');
pval=pvalvsNull_resid(x_label, y_label, readleft1,readright1,readleft2,readright2,permno,wsname,real_rho_path, resid_path_left, resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_GWC_CT_pval.txt'), pval);

% Compare GWC vs t1t2_ratio_25
x_label = 'GWC';
y_label = 'GM T1w/T2w';
readleft1 = readleft;
readright1 = readright;
readleft2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_left.csv');
readright2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_right.csv');
real_rho_path = strcat(analysis_directory,'means_GWC_t1t2_ratio_25_rho.txt');
resid_path_left = strcat(analysis_directory,'means_GWC_t1t2_ratio_25_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_GWC_t1t2_ratio_25_resid_right.txt');
pval=pvalvsNull_resid(x_label, y_label, readleft1,readright1,readleft2,readright2,permno,wsname,real_rho_path, resid_path_left, resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_GWC_t1t2_ratio_25_pval.txt'), pval);

% Compare GWC vs t1t2_ratio_25_WM
x_label = 'GWC';
y_label = 'SWM T1w/T2w';
readleft1 = readleft;
readright1 = readright;
readleft2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_WM_left.csv');
readright2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_WM_right.csv');
real_rho_path = strcat(analysis_directory,'means_GWC_t1t2_ratio_25_WM_rho.txt');
resid_path_left = strcat(analysis_directory,'means_GWC_t1t2_ratio_25_WM_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_GWC_t1t2_ratio_25_WM_resid_right.txt');
pval=pvalvsNull_resid(x_label, y_label, readleft1,readright1,readleft2,readright2,permno,wsname,real_rho_path, resid_path_left, resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_GWC_t1t2_ratio_25_WM_pval.txt'), pval);

% Compare GWC vs t1t2_ratio_50
x_label = 'GWC';
y_label = 'GM T1w/T2w (50%)';
readleft1 = readleft;
readright1 = readright;
readleft2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_50_left.csv');
readright2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_50_right.csv');
real_rho_path = strcat(analysis_directory,'means_GWC_t1t2_ratio_50_rho.txt');
resid_path_left = strcat(analysis_directory,'means_GWC_t1t2_ratio_50_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_GWC_t1t2_ratio_50_resid_right.txt');
pval=pvalvsNull_resid(x_label, y_label, readleft1,readright1,readleft2,readright2,permno,wsname,real_rho_path, resid_path_left, resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_GWC_t1t2_ratio_50_pval.txt'), pval);

%%%%%%%%%%%%%%%%%%%%%%%%%% CT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spin CT data
readleft = strcat(analysis_directory, 'mean_all_CT_left.csv');
readright = strcat(analysis_directory, 'mean_all_CT_right.csv');
permno = 1000; % how many spins
wsname = sprintf(strcat(analysis_directory, 'means_CT_rot.mat'));
%SpinPermuCIVET(readleft,readright,permno,wsname)

% Compare CT vs t1t2_ratio_25
x_label = 'CT';
y_label = 'GM T1w/T2w';
readleft1 = readleft;
readright1 = readright;
readleft2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_left.csv');
readright2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_right.csv');
real_rho_path = strcat(analysis_directory,'means_CT_t1t2_ratio_25_rho.txt');
resid_path_left = strcat(analysis_directory,'means_CT_t1t2_ratio_25_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_CT_t1t2_ratio_25_resid_right.txt');
pval=pvalvsNull_resid(x_label, y_label, readleft1,readright1,readleft2,readright2,permno,wsname,real_rho_path, resid_path_left, resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_CT_t1t2_ratio_25_pval.txt'), pval);

% Compare CT vs t1t2_ratio_25_WM
x_label = 'CT';
y_label = 'SWM T1w/T2w';
readleft1 = readleft;
readright1 = readright;
readleft2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_WM_left.csv');
readright2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_WM_right.csv');
real_rho_path = strcat(analysis_directory,'means_CT_t1t2_ratio_25_WM_rho.txt');
resid_path_left = strcat(analysis_directory,'means_CT_t1t2_ratio_25_WM_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_CT_t1t2_ratio_25_WM_resid_right.txt');
pval=pvalvsNull_resid(x_label, y_label, readleft1,readright1,readleft2,readright2,permno,wsname,real_rho_path, resid_path_left, resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_CT_t1t2_ratio_25_WM_pval.txt'), pval);

% Compare CT vs t1t2_ratio_50
x_label = 'CT';
y_label = 'GM T1w/T2w (50%)';
readleft1 = readleft;
readright1 = readright;
readleft2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_50_left.csv');
readright2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_50_right.csv');
real_rho_path = strcat(analysis_directory,'means_CT_t1t2_ratio_50_rho.txt');
resid_path_left = strcat(analysis_directory,'means_CT_t1t2_ratio_50_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_CT_t1t2_ratio_50_resid_right.txt');
pval=pvalvsNull_resid(x_label, y_label, readleft1,readright1,readleft2,readright2,permno,wsname,real_rho_path, resid_path_left, resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_CT_t1t2_ratio_50_pval.txt'), pval);

%%%%%%%%%%%%%%%%%%%%%%%%%% t1t2_ratio_25 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spin t1t2_ratio_25 data
readleft = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_left.csv');
readright = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_right.csv');
permno = 1000; % how many spins
wsname = sprintf(strcat(analysis_directory, 'means_t1t2_ratio_25_rot.mat'));
SpinPermuCIVET(readleft,readright,permno,wsname)

% Compare t1t2_ratio_25 vs t1t2_ratio_25_WM
x_label = 'GM T1w/T2w';
y_label = 'SWM T1w/T2w';
readleft1 = readleft;
readright1 = readright;
readleft2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_WM_left.csv');
readright2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_WM_right.csv');
real_rho_path = strcat(analysis_directory,'means_t1t2_ratio_25_t1t2_ratio_25_WM_rho.txt');
resid_path_left = strcat(analysis_directory,'means_t1t2_ratio_25_t1t2_ratio_25_WM_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_t1t2_ratio_25_t1t2_ratio_25_WM_resid_right.txt');
pval=pvalvsNull_resid(x_label, y_label, readleft1,readright1,readleft2,readright2,permno,wsname,real_rho_path, resid_path_left, resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_t1t2_ratio_25_t1t2_ratio_25_WM_pval.txt'), pval);

% Compare t1t2_ratio_25 vs t1t2_ratio_50
x_label = 'GM T1w/T2w';
y_label = 'GM T1w/T2w (50%)';
readleft1 = readleft;
readright1 = readright;
readleft2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_50_left.csv');
readright2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_50_right.csv');
real_rho_path = strcat(analysis_directory,'means_t1t2_ratio_25_t1t2_ratio_50_rho.txt');
resid_path_left = strcat(analysis_directory,'means_t1t2_ratio_25_t1t2_ratio_50_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_t1t2_ratio_25_t1t2_ratio_50_resid_right.txt');
pval=pvalvsNull_resid(x_label, y_label, readleft1,readright1,readleft2,readright2,permno,wsname,real_rho_path, resid_path_left, resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_t1t2_ratio_25_t1t2_ratio_50_pval.txt'), pval);

%%%%%%%%%%%%%%%%%%%%%%%%%% t1t2_ratio_25_WM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spin t1t2_ratio_25_WM data
readleft = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_WM_left.csv');
readright = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_WM_right.csv');
permno = 1000; % how many spins
wsname = sprintf(strcat(analysis_directory, 'means_t1t2_ratio_25_WM_rot.mat'));
SpinPermuCIVET(readleft,readright,permno,wsname)

% Compare t1t2_ratio_25_WM vs t1t2_ratio_50
x_label = 'SWM T1w/T2w';
y_label = 'GM T1w/T2w (50%)';
readleft1 = readleft;
readright1 = readright;
readleft2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_50_left.csv');
readright2 = strcat(analysis_directory, 'mean_all_t1t2_ratio_50_right.csv');
real_rho_path = strcat(analysis_directory,'means_t1t2_ratio_25_WM_t1t2_ratio_50_rho.txt');
resid_path_left = strcat(analysis_directory,'means_t1t2_ratio_25_WM_t1t2_ratio_50_resid_left.txt');
resid_path_right = strcat(analysis_directory,'means_t1t2_ratio_25_WM_t1t2_ratio_50_resid_right.txt');
pval=pvalvsNull_resid(x_label, y_label, readleft1,readright1,readleft2,readright2,permno,wsname,real_rho_path, resid_path_left, resid_path_right); 
dlmwrite(strcat(analysis_directory,'means_t1t2_ratio_25_WM_t1t2_ratio_50_pval.txt'), pval);


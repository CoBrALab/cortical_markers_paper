
analysis_directory = '../../AHBA/data_paper/';

%%%%%%%%%%%%%%%%%%%%%%%%%% BSC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spin BSC data
readleft = strcat(analysis_directory, 'mean_all_BSC_left_spin_test_DKT308_vertex.csv');
permno = 1000; % how many spins
wsname = sprintf(strcat(analysis_directory, 'means_BSC_rot.mat'));
SpinPermuCIVET_left(readleft,permno,wsname)

cells = ["astro" "endo" "micro" "neuroex" "neuroin" "oligo" "opc"];

% Compare BSC vs cells
x_label = 'BSC';
for i = 1:length(cells)
    disp(cells(i))
    y_label = cells(i);
    readleft1 = readleft;
    readleft2 = strcat(analysis_directory, 'cell_classes_AHBA_lh152parc_', cells(i),'_vertex.csv');
    real_rho_path = strcat(analysis_directory,'means_BSC_', cells(i),'_rho.txt');
    resid_path_left = strcat(analysis_directory,'means_BSC_', cells(i),'_resid_left.txt');
    pval=pvalvsNull_resid_left(x_label, y_label, readleft1,readleft2,permno,wsname,real_rho_path, resid_path_left);
    dlmwrite(strcat(analysis_directory,'means_BSC_', cells(i),'_pval.txt'), pval);
end

%%%%%%%%%%%%%%%%%%%%%%%%%% GWC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spin GWC data
readleft = strcat(analysis_directory, 'mean_all_GWC_left_spin_test_DKT308_vertex.csv');
permno = 1000; % how many spins
wsname = sprintf(strcat(analysis_directory, 'means_GWC_rot.mat'));
SpinPermuCIVET_left(readleft,permno,wsname)

cells = ["astro" "endo" "micro" "neuroex" "neuroin" "oligo" "opc"];

% Compare GWC vs cells
x_label = 'GWC';
for i = 1:length(cells)
    disp(cells(i))
    y_label = cells(i);
    readleft1 = readleft;
    readleft2 = strcat(analysis_directory, 'cell_classes_AHBA_lh152parc_', cells(i),'_vertex.csv');
    real_rho_path = strcat(analysis_directory,'means_GWC_', cells(i),'_rho.txt');
    resid_path_left = strcat(analysis_directory,'means_GWC_', cells(i),'_resid_left.txt');
    pval=pvalvsNull_resid_left(x_label, y_label, readleft1,readleft2,permno,wsname,real_rho_path, resid_path_left);
    dlmwrite(strcat(analysis_directory,'means_GWC_', cells(i),'_pval.txt'), pval);
end

%%%%%%%%%%%%%%%%%%%%%%%%%% CT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spin CT data
readleft = strcat(analysis_directory, 'mean_all_CT_left_spin_test_DKT308_vertex.csv');
permno = 1000; % how many spins
wsname = sprintf(strcat(analysis_directory, 'means_CT_rot.mat'));
SpinPermuCIVET_left(readleft,permno,wsname)

cells = ["astro" "endo" "micro" "neuroex" "neuroin" "oligo" "opc"];

% Compare CT vs cells
x_label = 'CT';
for i = 1:length(cells)
    disp(cells(i))
    y_label = cells(i);
    readleft1 = readleft;
    readleft2 = strcat(analysis_directory, 'cell_classes_AHBA_lh152parc_', cells(i),'_vertex.csv');
    real_rho_path = strcat(analysis_directory,'means_CT_', cells(i),'_rho.txt');
    resid_path_left = strcat(analysis_directory,'means_CT_', cells(i),'_resid_left.txt');
    pval=pvalvsNull_resid_left(x_label, y_label, readleft1,readleft2,permno,wsname,real_rho_path, resid_path_left);
    dlmwrite(strcat(analysis_directory,'means_CT_', cells(i),'_pval.txt'), pval);
end

%%%%%%%%%%%%%%%%%%%%%%%%%% t1t2_ratio_25 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spin t1t2_ratio_25 data
readleft = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_left_spin_test_DKT308_vertex.csv');
permno = 1000; % how many spins
wsname = sprintf(strcat(analysis_directory, 'means_t1t2_ratio_25_rot.mat'));
SpinPermuCIVET_left(readleft,permno,wsname)

cells = ["astro" "endo" "micro" "neuroex" "neuroin" "oligo" "opc"];

% Compare t1t2_ratio_25 vs cells
x_label = 't1t2_ratio_25';
for i = 1:length(cells)
    disp(cells(i))
    y_label = cells(i);
    readleft1 = readleft;
    readleft2 = strcat(analysis_directory, 'cell_classes_AHBA_lh152parc_', cells(i),'_vertex.csv');
    real_rho_path = strcat(analysis_directory,'means_t1t2_ratio_25_', cells(i),'_rho.txt');
    resid_path_left = strcat(analysis_directory,'means_t1t2_ratio_25_', cells(i),'_resid_left.txt');
    pval=pvalvsNull_resid_left(x_label, y_label, readleft1,readleft2,permno,wsname,real_rho_path, resid_path_left);
    dlmwrite(strcat(analysis_directory,'means_t1t2_ratio_25_', cells(i),'_pval.txt'), pval);
end

%%%%%%%%%%%%%%%%%%%%%%%%%% t1t2_ratio_50 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spin t1t2_ratio_50 data
readleft = strcat(analysis_directory, 'mean_all_t1t2_ratio_50_left_spin_test_DKT308_vertex.csv');
permno = 1000; % how many spins
wsname = sprintf(strcat(analysis_directory, 'means_t1t2_ratio_50_rot.mat'));
SpinPermuCIVET_left(readleft,permno,wsname)

cells = ["astro" "endo" "micro" "neuroex" "neuroin" "oligo" "opc"];

% Compare t1t2_ratio_50 vs cells
x_label = 't1t2_ratio_50';
for i = 1:length(cells)
    disp(cells(i))
    y_label = cells(i);
    readleft1 = readleft;
    readleft2 = strcat(analysis_directory, 'cell_classes_AHBA_lh152parc_', cells(i),'_vertex.csv');
    real_rho_path = strcat(analysis_directory,'means_t1t2_ratio_50_', cells(i),'_rho.txt');
    resid_path_left = strcat(analysis_directory,'means_t1t2_ratio_50_', cells(i),'_resid_left.txt');
    pval=pvalvsNull_resid_left(x_label, y_label, readleft1,readleft2,permno,wsname,real_rho_path, resid_path_left);
    dlmwrite(strcat(analysis_directory,'means_t1t2_ratio_50_', cells(i),'_pval.txt'), pval);
end

%%%%%%%%%%%%%%%%%%%%%%%%%% t1t2_ratio_25_WM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spin t1t2_ratio_25_WM data
readleft = strcat(analysis_directory, 'mean_all_t1t2_ratio_25_WM_left_spin_test_DKT308_vertex.csv');
permno = 1000; % how many spins
wsname = sprintf(strcat(analysis_directory, 'means_t1t2_ratio_25_WM_rot.mat'));
SpinPermuCIVET_left(readleft,permno,wsname)

cells = ["astro" "endo" "micro" "neuroex" "neuroin" "oligo" "opc"];

% Compare t1t2_ratio_25_WM vs cells
x_label = 't1t2_ratio_25_WM';
for i = 1:length(cells)
    disp(cells(i))
    y_label = cells(i);
    readleft1 = readleft;
    readleft2 = strcat(analysis_directory, 'cell_classes_AHBA_lh152parc_', cells(i),'_vertex.csv');
    real_rho_path = strcat(analysis_directory,'means_t1t2_ratio_25_WM_', cells(i),'_rho.txt');
    resid_path_left = strcat(analysis_directory,'means_t1t2_ratio_25_WM_', cells(i),'_resid_left.txt');
    pval=pvalvsNull_resid_left(x_label, y_label, readleft1,readleft2,permno,wsname,real_rho_path, resid_path_left);
    dlmwrite(strcat(analysis_directory,'means_t1t2_ratio_25_WM_', cells(i),'_pval.txt'), pval);
end

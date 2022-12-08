
analysis_directory = '../../descriptive_stats_AHBA/';

%%%%%%%%%%%%%%%%%%%%%%%%%% MEANS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spin BSC data
readleft = strcat(analysis_directory, '/results_parc/mean_all_BSC_left_spin_test_DKT308_vertex.csv');
permno = 1000; % how many spins
wsname = sprintf(strcat(analysis_directory, '/spin_test_vertex_wise/means_BSC_rot.mat'));
%SpinPermuCIVET_left(readleft,permno,wsname)

cells = ["astro" "endo" "micro" "neuroex" "neuroin" "oligo" "opc"];
%cells = ["endo"];

% Compare BSC vs cells
x_label = 'BSC';
for i = 1:length(cells)
    disp(cells(i))
    y_label = cells(i);
    readleft1 = readleft;
    readleft2 = strcat(analysis_directory, 'AHBA_cell_type/cell_classes_AHBA_lh152parc_', cells(i),'_vertex.csv');
    real_rho_path = strcat(analysis_directory,'spin_test_vertex_wise/means_BSC_', cells(i),'_rho.txt');
    resid_path_left = strcat(analysis_directory,'spin_test_vertex_wise/means_BSC_', cells(i),'_resid_left.txt');
    pval=pvalvsNull_resid_left(x_label, y_label, readleft1,readleft2,permno,wsname,real_rho_path, resid_path_left);
    dlmwrite(strcat(analysis_directory,'spin_test_vertex_wise/means_BSC_', cells(i),'_pval.txt'), pval);
end

% Spin GWC data
readleft = strcat(analysis_directory, '/results_parc/mean_all_GWC_left_spin_test_DKT308_vertex.csv');
permno = 1000; % how many spins
wsname = sprintf(strcat(analysis_directory, '/spin_test_vertex_wise/means_GWC_rot.mat'));
SpinPermuCIVET_left(readleft,permno,wsname)

cells = ["astro" "endo" "micro" "neuroex" "neuroin" "oligo" "opc"];
%cells = ["endo"];

% Compare GWC vs cells
x_label = 'GWC';
for i = 1:length(cells)
    disp(cells(i))
    y_label = cells(i);
    readleft1 = readleft;
    readleft2 = strcat(analysis_directory, 'AHBA_cell_type/cell_classes_AHBA_lh152parc_', cells(i),'_vertex.csv');
    real_rho_path = strcat(analysis_directory,'spin_test_vertex_wise/means_GWC_', cells(i),'_rho.txt');
    resid_path_left = strcat(analysis_directory,'spin_test_vertex_wise/means_GWC_', cells(i),'_resid_left.txt');
    pval=pvalvsNull_resid_left(x_label, y_label, readleft1,readleft2,permno,wsname,real_rho_path, resid_path_left);
    dlmwrite(strcat(analysis_directory,'spin_test_vertex_wise/means_GWC_', cells(i),'_pval.txt'), pval);
end

% Spin CT data
readleft = strcat(analysis_directory, '/results_parc/mean_all_CT_left_spin_test_DKT308_vertex.csv');
permno = 1000; % how many spins
wsname = sprintf(strcat(analysis_directory, '/spin_test_vertex_wise/means_CT_rot.mat'));
SpinPermuCIVET_left(readleft,permno,wsname)

cells = ["astro" "endo" "micro" "neuroex" "neuroin" "oligo" "opc"];
%cells = ["endo"];

% Compare CT vs cells
x_label = 'CT';
for i = 1:length(cells)
    disp(cells(i))
    y_label = cells(i);
    readleft1 = readleft;
    readleft2 = strcat(analysis_directory, 'AHBA_cell_type/cell_classes_AHBA_lh152parc_', cells(i),'_vertex.csv');
    real_rho_path = strcat(analysis_directory,'spin_test_vertex_wise/means_CT_', cells(i),'_rho.txt');
    resid_path_left = strcat(analysis_directory,'spin_test_vertex_wise/means_CT_', cells(i),'_resid_left.txt');
    pval=pvalvsNull_resid_left(x_label, y_label, readleft1,readleft2,permno,wsname,real_rho_path, resid_path_left);
    dlmwrite(strcat(analysis_directory,'spin_test_vertex_wise/means_CT_', cells(i),'_pval.txt'), pval);
end

% Spin t1t2_ratio_25 data
readleft = strcat(analysis_directory, '/results_parc/mean_all_t1t2_ratio_25_left_spin_test_DKT308_vertex.csv');
permno = 1000; % how many spins
wsname = sprintf(strcat(analysis_directory, '/spin_test_vertex_wise/means_t1t2_ratio_25_rot.mat'));
SpinPermuCIVET_left(readleft,permno,wsname)

cells = ["astro" "endo" "micro" "neuroex" "neuroin" "oligo" "opc"];
%cells = ["endo"];

% Compare t1t2_ratio_25 vs cells
x_label = 't1t2_ratio_25';
for i = 1:length(cells)
    disp(cells(i))
    y_label = cells(i);
    readleft1 = readleft;
    readleft2 = strcat(analysis_directory, 'AHBA_cell_type/cell_classes_AHBA_lh152parc_', cells(i),'_vertex.csv');
    real_rho_path = strcat(analysis_directory,'spin_test_vertex_wise/means_t1t2_ratio_25_', cells(i),'_rho.txt');
    resid_path_left = strcat(analysis_directory,'spin_test_vertex_wise/means_t1t2_ratio_25_', cells(i),'_resid_left.txt');
    pval=pvalvsNull_resid_left(x_label, y_label, readleft1,readleft2,permno,wsname,real_rho_path, resid_path_left);
    dlmwrite(strcat(analysis_directory,'spin_test_vertex_wise/means_t1t2_ratio_25_', cells(i),'_pval.txt'), pval);
end

% Spin t1t2_ratio_50 data
readleft = strcat(analysis_directory, '/results_parc/mean_all_t1t2_ratio_50_left_spin_test_DKT308_vertex.csv');
permno = 1000; % how many spins
wsname = sprintf(strcat(analysis_directory, '/spin_test_vertex_wise/means_t1t2_ratio_50_rot.mat'));
SpinPermuCIVET_left(readleft,permno,wsname)

cells = ["astro" "endo" "micro" "neuroex" "neuroin" "oligo" "opc"];
%cells = ["endo"];

% Compare t1t2_ratio_50 vs cells
x_label = 't1t2_ratio_50';
for i = 1:length(cells)
    disp(cells(i))
    y_label = cells(i);
    readleft1 = readleft;
    readleft2 = strcat(analysis_directory, 'AHBA_cell_type/cell_classes_AHBA_lh152parc_', cells(i),'_vertex.csv');
    real_rho_path = strcat(analysis_directory,'spin_test_vertex_wise/means_t1t2_ratio_50_', cells(i),'_rho.txt');
    resid_path_left = strcat(analysis_directory,'spin_test_vertex_wise/means_t1t2_ratio_50_', cells(i),'_resid_left.txt');
    pval=pvalvsNull_resid_left(x_label, y_label, readleft1,readleft2,permno,wsname,real_rho_path, resid_path_left);
    dlmwrite(strcat(analysis_directory,'spin_test_vertex_wise/means_t1t2_ratio_50_', cells(i),'_pval.txt'), pval);
end

% Spin t1t2_ratio_25_WM data
readleft = strcat(analysis_directory, '/results_parc/mean_all_t1t2_ratio_25_WM_left_spin_test_DKT308_vertex.csv');
permno = 1000; % how many spins
wsname = sprintf(strcat(analysis_directory, '/spin_test_vertex_wise/means_t1t2_ratio_25_WM_rot.mat'));
SpinPermuCIVET_left(readleft,permno,wsname)

cells = ["astro" "endo" "micro" "neuroex" "neuroin" "oligo" "opc"];
%cells = ["endo"];

% Compare t1t2_ratio_25_WM vs cells
x_label = 't1t2_ratio_25_WM';
for i = 1:length(cells)
    disp(cells(i))
    y_label = cells(i);
    readleft1 = readleft;
    readleft2 = strcat(analysis_directory, 'AHBA_cell_type/cell_classes_AHBA_lh152parc_', cells(i),'_vertex.csv');
    real_rho_path = strcat(analysis_directory,'spin_test_vertex_wise/means_t1t2_ratio_25_WM_', cells(i),'_rho.txt');
    resid_path_left = strcat(analysis_directory,'spin_test_vertex_wise/means_t1t2_ratio_25_WM_', cells(i),'_resid_left.txt');
    pval=pvalvsNull_resid_left(x_label, y_label, readleft1,readleft2,permno,wsname,real_rho_path, resid_path_left);
    dlmwrite(strcat(analysis_directory,'spin_test_vertex_wise/means_t1t2_ratio_25_WM_', cells(i),'_pval.txt'), pval);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Compare cell-type between themselves %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spin astro data
readleft = strcat(analysis_directory, 'AHBA_cell_type/cell_classes_AHBA_lh152parc_astro.txt');
permno = 1000; % how many spins
wsname = sprintf(strcat(analysis_directory, '/spin_test_vertex_wise/means_AHBA_astro_rot.mat'));
SpinPermuCIVET_left(readleft,permno,wsname)

% Compare astro with cells
cells = ["endo" "micro" "neuroex" "neuroin" "oligo" "opc"];
x_label = 'astro';
for i = 1:length(cells)
    disp(cells(i))
    y_label = cells(i);
    readleft1 = readleft;
    readleft2 = strcat(analysis_directory, 'AHBA_cell_type/cell_classes_AHBA_lh152parc_', cells(i),'_vertex.csv');
    real_rho_path = strcat(analysis_directory,'spin_test_vertex_wise/means_AHBA_astro_', cells(i),'_rho.txt');
    resid_path_left = strcat(analysis_directory,'spin_test_vertex_wise/means_AHBA_astro_', cells(i),'_resid_left.txt');
    pval=pvalvsNull_resid_left(x_label, y_label, readleft1,readleft2,permno,wsname,real_rho_path, resid_path_left);
    dlmwrite(strcat(analysis_directory,'spin_test_vertex_wise/means_AHBA_astro_', cells(i),'_pval.txt'), pval);
end

% Spin endo data
readleft = strcat(analysis_directory, 'AHBA_cell_type/cell_classes_AHBA_lh152parc_endo.txt');
permno = 1000; % how many spins
wsname = sprintf(strcat(analysis_directory, '/spin_test_vertex_wise/means_AHBA_endo_rot.mat'));
SpinPermuCIVET_left(readleft,permno,wsname)

% Compare endo with cells
cells = ["micro" "neuroex" "neuroin" "oligo" "opc"];
x_label = 'endo';
for i = 1:length(cells)
    disp(cells(i))
    y_label = cells(i);
    readleft1 = readleft;
    readleft2 = strcat(analysis_directory, 'AHBA_cell_type/cell_classes_AHBA_lh152parc_', cells(i),'_vertex.csv');
    real_rho_path = strcat(analysis_directory,'spin_test_vertex_wise/means_AHBA_endo_', cells(i),'_rho.txt');
    resid_path_left = strcat(analysis_directory,'spin_test_vertex_wise/means_AHBA_endo_', cells(i),'_resid_left.txt');
    pval=pvalvsNull_resid_left(x_label, y_label, readleft1,readleft2,permno,wsname,real_rho_path, resid_path_left);
    dlmwrite(strcat(analysis_directory,'spin_test_vertex_wise/means_AHBA_endo_', cells(i),'_pval.txt'), pval);
end

% Spin micro data
readleft = strcat(analysis_directory, 'AHBA_cell_type/cell_classes_AHBA_lh152parc_micro.txt');
permno = 1000; % how many spins
wsname = sprintf(strcat(analysis_directory, '/spin_test_vertex_wise/means_AHBA_micro_rot.mat'));
SpinPermuCIVET_left(readleft,permno,wsname)

% Compare micro with cells
cells = ["neuroex" "neuroin" "oligo" "opc"];
x_label = 'micro';
for i = 1:length(cells)
    disp(cells(i))
    y_label = cells(i);
    readleft1 = readleft;
    readleft2 = strcat(analysis_directory, 'AHBA_cell_type/cell_classes_AHBA_lh152parc_', cells(i),'_vertex.csv');
    real_rho_path = strcat(analysis_directory,'spin_test_vertex_wise/means_AHBA_micro_', cells(i),'_rho.txt');
    resid_path_left = strcat(analysis_directory,'spin_test_vertex_wise/means_AHBA_micro_', cells(i),'_resid_left.txt');
    pval=pvalvsNull_resid_left(x_label, y_label, readleft1,readleft2,permno,wsname,real_rho_path, resid_path_left);
    dlmwrite(strcat(analysis_directory,'spin_test_vertex_wise/means_AHBA_micro_', cells(i),'_pval.txt'), pval);
end

% Spin neuroex data
readleft = strcat(analysis_directory, 'AHBA_cell_type/cell_classes_AHBA_lh152parc_neuroex.txt');
permno = 1000; % how many spins
wsname = sprintf(strcat(analysis_directory, '/spin_test_vertex_wise/means_AHBA_neuroex_rot.mat'));
SpinPermuCIVET_left(readleft,permno,wsname)

% Compare neuroex with cells
cells = ["neuroin" "oligo" "opc"];
x_label = 'neuroex';
for i = 1:length(cells)
    disp(cells(i))
    y_label = cells(i);
    readleft1 = readleft;
    readleft2 = strcat(analysis_directory, 'AHBA_cell_type/cell_classes_AHBA_lh152parc_', cells(i),'_vertex.csv');
    real_rho_path = strcat(analysis_directory,'spin_test_vertex_wise/means_AHBA_neuroex_', cells(i),'_rho.txt');
    resid_path_left = strcat(analysis_directory,'spin_test_vertex_wise/means_AHBA_neuroex_', cells(i),'_resid_left.txt');
    pval=pvalvsNull_resid_left(x_label, y_label, readleft1,readleft2,permno,wsname,real_rho_path, resid_path_left);
    dlmwrite(strcat(analysis_directory,'spin_test_vertex_wise/means_AHBA_neuroex_', cells(i),'_pval.txt'), pval);
end

% Spin neuroin data
readleft = strcat(analysis_directory, 'AHBA_cell_type/cell_classes_AHBA_lh152parc_neuroin.txt');
permno = 1000; % how many spins
wsname = sprintf(strcat(analysis_directory, '/spin_test_vertex_wise/means_AHBA_neuroin_rot.mat'));
SpinPermuCIVET_left(readleft,permno,wsname)

% Compare neuroin with cells
cells = ["oligo" "opc"];
x_label = 'neuroin';
for i = 1:length(cells)
    disp(cells(i))
    y_label = cells(i);
    readleft1 = readleft;
    readleft2 = strcat(analysis_directory, 'AHBA_cell_type/cell_classes_AHBA_lh152parc_', cells(i),'_vertex.csv');
    real_rho_path = strcat(analysis_directory,'spin_test_vertex_wise/means_AHBA_neuroin_', cells(i),'_rho.txt');
    resid_path_left = strcat(analysis_directory,'spin_test_vertex_wise/means_AHBA_neuroin_', cells(i),'_resid_left.txt');
    pval=pvalvsNull_resid_left(x_label, y_label, readleft1,readleft2,permno,wsname,real_rho_path, resid_path_left);
    dlmwrite(strcat(analysis_directory,'spin_test_vertex_wise/means_AHBA_neuroin_', cells(i),'_pval.txt'), pval);
end

% Spin oligo data
readleft = strcat(analysis_directory, 'AHBA_cell_type/cell_classes_AHBA_lh152parc_oligo.txt');
permno = 1000; % how many spins
wsname = sprintf(strcat(analysis_directory, '/spin_test_vertex_wise/means_AHBA_oligo_rot.mat'));
SpinPermuCIVET_left(readleft,permno,wsname)

% Compare oligo with cells
cells = ["opc"];
x_label = 'oligo';
for i = 1:length(cells)
    disp(cells(i))
    y_label = cells(i);
    readleft1 = readleft;
    readleft2 = strcat(analysis_directory, 'AHBA_cell_type/cell_classes_AHBA_lh152parc_', cells(i),'_vertex.csv');
    real_rho_path = strcat(analysis_directory,'spin_test_vertex_wise/means_AHBA_oligo_', cells(i),'_rho.txt');
    resid_path_left = strcat(analysis_directory,'spin_test_vertex_wise/means_AHBA_oligo_', cells(i),'_resid_left.txt');
    pval=pvalvsNull_resid_left(x_label, y_label, readleft1,readleft2,permno,wsname,real_rho_path, resid_path_left);
    dlmwrite(strcat(analysis_directory,'spin_test_vertex_wise/means_AHBA_oligo_', cells(i),'_pval.txt'), pval);
end

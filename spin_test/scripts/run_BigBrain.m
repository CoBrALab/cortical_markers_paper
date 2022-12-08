
analysis_directory = '../../BigBrain/';

%%%%%%%%%%%%%%%%% BigBrain 25_GM

% Spin BigBrain_25_GM data
readleft = strcat(analysis_directory, 'inputs/BigBrain/BigBrain_GM_25_left_div100_20mm_40962_inv.txt');
readright = strcat(analysis_directory, 'inputs/BigBrain/BigBrain_GM_25_right_div100_20mm_40962_inv.txt');
permno = 1000; % how many spins
wsname = sprintf(strcat(analysis_directory, 'spin_test/BigBrain_GM_25_div100_20mm_40962_inv.mat'));
% SpinPermuCIVET(readleft,readright,permno,wsname)

% Compare BigBrain_25_GM to MRI markers
markers = ["BSC" "GWC" "CT" "t1t2_ratio_25_GM" "t1t2_ratio_50_GM" "t1t2_ratio_25_WM"];

x_label = 'GM BigBrain';
markers_label = ["BSC" "GWC" "CT" "GM T1w/T2w" "GM T1w/T2w (50%)" "SWM T1w/T2w"];

for i = 1:length(markers)
    disp(markers(i))
    y_label = markers_label(i);
    readleft2 = strcat(analysis_directory, 'inputs/descriptive_stats/mean_all_',markers(i),'_left_spin_test.csv');
    readright2 = strcat(analysis_directory, 'inputs/descriptive_stats/mean_all_',markers(i),'_right_spin_test.csv');
    real_rho_path = strcat(analysis_directory,'spin_test/BB25GM_vs_',markers(i),'_rho.txt');
    resid_path_left = strcat(analysis_directory,'spin_test/BB25GM_vs_',markers(i),'_resid_left.txt');
    resid_path_right = strcat(analysis_directory,'spin_test/BB25GM_vs_',markers(i),'_resid_right.txt');
    pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right);
    dlmwrite(strcat(analysis_directory,'spin_test/BB25GM_vs_',markers(i),'_pval.txt'), pval);
end

%%%%%%%%%%%%%%%%% BigBrain 50_GM

% Spin BigBrain_50_GM data
readleft = strcat(analysis_directory, 'inputs/BigBrain/BigBrain_GM_50_left_div100_20mm_40962_inv.txt');
readright = strcat(analysis_directory, 'inputs/BigBrain/BigBrain_GM_50_right_div100_20mm_40962_inv.txt');
permno = 1000; % how many spins
wsname = sprintf(strcat(analysis_directory, 'spin_test/BigBrain_GM_50_div100_20mm_40962_inv.mat'));
%SpinPermuCIVET(readleft,readright,permno,wsname)

% Compare BigBrain_50_GM to MRI markers
markers = ["BSC" "GWC" "CT" "t1t2_ratio_25_GM" "t1t2_ratio_50_GM" "t1t2_ratio_25_WM"];

x_label = 'GM BigBrain (50%)';
markers_label = ["BSC" "GWC" "CT" "GM T1w/T2w" "GM T1w/T2w (50%)" "SWM T1w/T2w"];

for i = 1:length(markers)
    disp(markers(i))
    y_label = markers_label(i);
    readleft2 = strcat(analysis_directory, 'inputs/descriptive_stats/mean_all_',markers(i),'_left_spin_test.csv');
    readright2 = strcat(analysis_directory, 'inputs/descriptive_stats/mean_all_',markers(i),'_right_spin_test.csv');
    real_rho_path = strcat(analysis_directory,'spin_test/BB50GM_vs_',markers(i),'_rho.txt');
    resid_path_left = strcat(analysis_directory,'spin_test/BB50GM_vs_',markers(i),'_resid_left.txt');
    resid_path_right = strcat(analysis_directory,'spin_test/BB50GM_vs_',markers(i),'_resid_right.txt');
    pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right);
    dlmwrite(strcat(analysis_directory,'spin_test/BB50GM_vs_',markers(i),'_pval.txt'), pval);
end

%%%%%%%%%%%%%%%%% BigBrain 25_WM

% Spin BigBrain_25_WM data
readleft = strcat(analysis_directory, 'inputs/BigBrain/BigBrain_WM_25_left_div100_20mm_40962_inv.txt');
readright = strcat(analysis_directory, 'inputs/BigBrain/BigBrain_WM_25_right_div100_20mm_40962_inv.txt');
permno = 1000; % how many spins
wsname = sprintf(strcat(analysis_directory, 'spin_test/BigBrain_WM_25_div100_20mm_40962_inv.mat'));
%SpinPermuCIVET(readleft,readright,permno,wsname)

% Compare BigBrain_25_WM to MRI markers
markers = ["BSC" "GWC" "CT" "t1t2_ratio_25_GM" "t1t2_ratio_50_GM" "t1t2_ratio_25_WM"];

x_label = 'SWM BigBrain';
markers_label = ["BSC" "GWC" "CT" "GM T1w/T2w" "GM T1w/T2w (50%)" "SWM T1w/T2w"];

for i = 1:length(markers)
    disp(markers(i))
    y_label = markers_label(i);
    readleft2 = strcat(analysis_directory, 'inputs/descriptive_stats/mean_all_',markers(i),'_left_spin_test.csv');
    readright2 = strcat(analysis_directory, 'inputs/descriptive_stats/mean_all_',markers(i),'_right_spin_test.csv');
    real_rho_path = strcat(analysis_directory,'spin_test/BB25WM_vs_',markers(i),'_rho.txt');
    resid_path_left = strcat(analysis_directory,'spin_test/BB25WM_vs_',markers(i),'_resid_left.txt');
    resid_path_right = strcat(analysis_directory,'spin_test/BB25WM_vs_',markers(i),'_resid_right.txt');
    pval=pvalvsNull_resid(x_label,y_label,readleft,readright,readleft2,readright2,permno,wsname,real_rho_path,resid_path_left,resid_path_right);
    dlmwrite(strcat(analysis_directory,'spin_test/BB25WM_vs_',markers(i),'_pval.txt'), pval);
end



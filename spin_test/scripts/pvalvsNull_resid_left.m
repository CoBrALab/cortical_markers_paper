function pval=pvalvsNull_resid_left(x_label,y_label,readleft1,readleft2,permno,wsname, real_rho_path, resid_path_left)
% Calculate the p-value of correlation between two surface maps based on
% the null distribution of spins of map 1
% FORMAT pvalvsNull(readleft1,readright1,readleft2,readright2,permno,wsname)
% readleft1     - the filename of the first left surface data to spin 
% readright1    - the filename of the first right surface data to spin 
% readleft2     - the filename of the second left surface data to spin 
% readright2    - the filename of the second right surface data to spin 
% permno       - the number of permutations used in SpinPermuFS/CIVET
% wsname       - the name of a workspace file output from SpinPermuFS/CIVET
% pval         - the calculated p-value
% Added 07/31/2020 (SMW): indicate vertices to exclude (e.g., medial wall)
% v_exclude_left  - left hemisphere vertices to exclude (indicate with 1)
% v_exclude_right - right hemisphere vertices to exclude (indicate with 1)
% Example   p=pvalvsNull('../data/depressionFSdataL.csv','../data/depressionFSdataR.csv','../data/ptsdFSdataL.csv','../data/ptsdFSdataR.csv',100,'../data/rotationFS.mat')
% will calculate the pvalue of correlation between prebuilt data, neurosynth map associated with 'depression',
% and 'PTSD' using the null distribution of depression maps spun 100 times
% from the SpinPermuFS.m
% Simiarly, it can be used for CIVET version as well.
% Aaron Alexander-Bloch & Siyuan Liu 
% pvalvsNull.m, 2018-04-22

%load the saved workspace from SpinPermu
load(wsname);

%read the data saved in csv and merge left and right surfaces into one
datal1=importdata(readleft1); % .data() part may or may not be needed
%datar1=importdata(readright1);
datal2=importdata(readleft2);
%datar2=importdata(readright2);

x_label = x_label;
y_label = y_label;

% Label medial wall vertices with NaN (07/31/2020):
%%%%%%%%%%%% MODIFIED MEDIAL WALL MASK FOR PARCELLATION %%%%%%%%%%%%%%%%%%%
v_exclude_left = logical(load('../data/CIVET_2.0_mask_left_short_inv_AHBA.txt')); % label of vertices in the medial wall is 1644825
%v_exclude_right = logical(load('../data/CIVET_2.0_mask_right_short_inv.txt'));
datal1(v_exclude_left) = NaN;
datal2(v_exclude_left) = NaN;
%datar1(v_exclude_right) = NaN;
%datar2(v_exclude_right) = NaN;

%data1=cat(1,datal1,datar1);
data1=datal1;
%data2=cat(1,datal2,datar2);
data2=datal2;

% Z-score input data
zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan')); % Function to Z-score surfaces with the NaNs from medial wall removal
Z_data1 = zscor_xnan(data1); %Z-scored values of first map
Z_data2 = zscor_xnan(data2); %Z-scored values of second map

%calculate the real Pearson's correlation between two interested maps
realrho=corr(data1,data2, 'rows','complete'); % 'rows','complete' to exclude NaN's
Z_realrho = corr(Z_data1,Z_data2, 'rows','complete'); %Sanity check to see if correlation with original vs z-scored data is the same (it should be)

% Get residuals from regression
[b,bint,resid] = regress(Z_data1,Z_data2); % Get residuals from regression of z-scored surfaces
resid_left = resid(1:40962,:); % Divide into left and right
%resid_right = resid(40963:81924,:);
%resid_left(v_exclude_left) = 0; %residual values on medial wall = 0
%resid_right(v_exclude_right) = 0;
dlmwrite(resid_path_left, resid_left); %write residuals as vertex files (.txt)
%dlmwrite(resid_path_right, resid_right);

% Scatter plot of vertices and regression line
model = fitlm(Z_data1,Z_data2); % Regression function on z-scored surfaces
figure('visible','off','DefaultAxesFontSize',18)

z_data2_above = NaN(40962,1);
z_data2_inside = NaN(40962,1);
z_data2_below = NaN(40962,1);

for i = 1:40962
   if Z_data2(i,1) > (table2array(model.Coefficients(1,1))+1) + Z_data1(i,1)*table2array(model.Coefficients(2,1))
       z_data2_above(i,1) = Z_data2(i,1);
   elseif Z_data2(i,1) < (table2array(model.Coefficients(1,1))-1) + Z_data1(i,1)*table2array(model.Coefficients(2,1))
       z_data2_below(i,1) = Z_data2(i,1);
   else
       z_data2_inside(i,1) = Z_data2(i,1);
   end
end

scatter(Z_data1, z_data2_above, 0.1, 'r'); %generate graph
hold on
scatter(Z_data1, z_data2_inside, 0.1, 'k'); %generate graph
scatter(Z_data1, z_data2_below, 0.1, 'b'); %generate graph

reg_line = table2array(model.Coefficients(1,1)) + Z_data1*table2array(model.Coefficients(2,1));
reg_line_top = (table2array(model.Coefficients(1,1))+1) + Z_data1*table2array(model.Coefficients(2,1));
reg_line_bot = (table2array(model.Coefficients(1,1))-1) + Z_data1*table2array(model.Coefficients(2,1));

plot(Z_data1, reg_line, 'LineWidth', 2, 'Color', 'k');
plot(Z_data1, reg_line_top, 'LineWidth', 2, 'Color', 'r');
plot(Z_data1, reg_line_bot, 'LineWidth', 2, 'Color', 'b');

% x_label = extractBefore(extractAfter(readleft1, 'mean_all_'), '_left_spin_test.csv');
% y_label = extractBefore(extractAfter(readleft2, 'mean_all_'), '_left_spin_test.csv');
xlabel(x_label,'FontSize', 25);
ylabel(y_label,'FontSize', 25);
title('');
legend('off')
graph_path = strcat(extractBefore(resid_path_left, '_resid'), '_graph.png'); %where to save graph
saveas(gcf, graph_path) %save graph as .png

% test the observed rho against null described by SpinPermu
nullrho=[];
for i=1:permno
%tempdata=cat(2,bigrotl(i,:),bigrotr(i,:))';
tempdata=cat(2,bigrotl(i,:))';
nullrho=cat(1,nullrho,corr(tempdata,data2, 'rows','complete')); % 'rows','complete' to exclude NaN's
end
%assuming sign is preserved, calculate the probability that the observed
%correlation coeffcient is above the null distribution
pval=length(find(abs(nullrho)>abs(realrho)))/permno; % added abs() 07/31/2020pval

% create & save histogram of randomly generated correlations (can be valuable downstream):
figure('visible','off','DefaultAxesFontSize',18)
h=figure;
hist(nullrho,permno);
xlabel(sprintf("Pearson R"))
title(sprintf('Pearson R= %d;  p= %d',realrho,pval))
graph_path = strcat(extractBefore(resid_path_left, '_resid'), '_hist.png');
saveas(h, graph_path)
%dlmwrite(strcat(output_dir,'correspondence_stat_null.txt'), nullrho);

coeffs = [realrho; Z_realrho; b; table2array(model.Coefficients(2,1))]; %all coefficients should be the same (correlation of original surfaces, correlation of z-scored surfaces, regression on z-scored surfaces for residuals, regression on z-scored surfaces for vertex plot)
dlmwrite(real_rho_path, coeffs);


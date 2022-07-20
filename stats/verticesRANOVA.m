% t-tests for number of reliable vertices in different conditions
% Written by C. Yilmaz, Boyaci Lab, 2019
% MATLAB R2016b
% ----------------------------------------------------------------------- %
clear all; close all;
EstimateFolder = '/Users/cemreyilmaz/Documents/fear-pRF/data/estimates/';
StatFolder = '/Users/cemreyilmaz/Documents/fear-pRF/results/';
% import data
cd(EstimateFolder)
data = readtable('Vertices_Eccentricity_vs_Sigma.csv');
within = table([1 2 3]','VariableNames',{'Conditions'});
rm = fitrm(data,'Scrambled-Emotional~Hemispheres','WithinDesign',within); 
ranovatbl = ranova(rm);
% ranavotbl cannot be saved as it is. So, I rearranged it.
ranovatbl = table(ranovatbl);
ranovatbl=[ranovatbl.Properties.RowNames,ranovatbl.ranovatbl(1,:);ranovatbl.ranovatbl(2,:);...
    ranovatbl.ranovatbl(3,:)];
tbl = multcompare(rm,'Hemispheres','Alpha', 0.01, ...
    'ComparisonType', 'bonferroni');
% save the tables as csv
writetable(ranovatbl,[StatFolder, 'ranova_vertices.csv'],'WriteRowNames',true)
writetable(tbl,[StatFolder, 'posthoc_multcomp_vertices.csv'],'WriteRowNames',true)

% t-tests for number of reliable vertices in different conditions
% Written by C. Yilmaz, Boyaci Lab, 2019
% MATLAB R2016b
% ----------------------------------------------------------------------- %
% effect size toolbox
% https://www.mathworks.com/matlabcentral/fileexchange/32398-hhentschke-measures-of-effect-size-toolbox?w.mathworks.com
addpath(genpath('/Users/cemreyilmaz/Documents/fear-pRF/Codes/statAnalysis/hhentschke_measures-of-effect-size-toolbox'))
% ----------------------------------------------------------------------- %
clear all; close all;
EstimateFolder = '/Users/cemreyilmaz/Documents/fear-pRF/data/estimates/';
StatFolder = '/Users/cemreyilmaz/Documents/fear-pRF/results/';
% import data
data = readtable([EstimateFolder 'Vertices_Eccentricity_vs_Sigma.csv']);
within = table([1 2 3]','VariableNames',{'Conditions'});
rm = fitrm(data,'Scrambled-Emotional~Hemispheres','WithinDesign',within); 
stats = mes1way([data.Scrambled,data.Neutral,data.Emotional],'eta2');
p = table({'isDep';'confLevel';'eta2'}, ...
    [stats.isDep;stats.confLevel;stats.eta2]);
ranovatbl = ranova(rm);
% ranavotbl cannot be saved as it is. So, I rearranged it.
ranovatbl = table(ranovatbl);
ranovatbl=[ranovatbl.Properties.RowNames,ranovatbl.ranovatbl(1,:);ranovatbl.ranovatbl(2,:);...
    ranovatbl.ranovatbl(3,:)];
ranovatbl = [ranovatbl,p];
tbl1 = multcompare(rm,'Hemispheres','Alpha', 0.05, ...
    'ComparisonType', 'bonferroni');
tbl2 = multcompare(rm,'Conditions','Alpha', 0.05, ...
    'ComparisonType', 'bonferroni');
stats12 = mes(data.Scrambled,data.Neutral,'hedgesg','nBoot',10000);
stats13 = mes(data.Scrambled,data.Emotional,'hedgesg','nBoot',10000);
stats23 = mes(data.Neutral,data.Emotional,'hedgesg','nBoot',10000);
cohensd = table([stats12.hedgesg;stats13.hedgesg;stats12.hedgesg;...
    stats23.hedgesg;stats13.hedgesg;stats23.hedgesg],'VariableNames',{'cohens-d'});
tbl2 = [tbl2,cohensd];
% save the tables as csv
writetable(ranovatbl,[StatFolder, 'ranova_vertices.csv'],'WriteRowNames',true)
writetable(tbl1,[StatFolder, 'posthoc_hemispheres_vertices.csv'],'WriteRowNames',true)
writetable(tbl2,[StatFolder, 'posthoc_conditions_vertices.csv'],'WriteRowNames',true)

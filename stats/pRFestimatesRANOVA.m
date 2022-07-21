% Repeated Measures ANOVA of the pRF estimations
% Written by C. Yilmaz, Boyaci Lab, 2019
% MATLAB R2016b
% This code uses the outputs of pRFestimates.m
% Uncomment rmANOVA or post-hoc for analyses. Find them into the code.
% You can save the results as csv by uncommenting the last two parts of code.
% ----------------------------------------------------------------------- %
clear all; close all;
Val = 'Sigma'; % Sigma: pRF size / Beta: BOLD amplitude
StatFolder = '/Users/cemreyilmaz/Documents/fear-pRF/results';
EstFolder = '/Users/cemreyilmaz/Documents/fear-pRF/data/estimates';
cd(EstFolder)
var = 'Hemispheres'; % post-hoc variable: ROI / Hemispheres
% ----------------------------------------------------------------------- %
% effect size toolbox
% https://www.mathworks.com/matlabcentral/fileexchange/32398-hhentschke-measures-of-effect-size-toolbox?w.mathworks.com
addpath(genpath('/Users/cemreyilmaz/Documents/fear-pRF/Codes/statAnalysis/hhentschke_measures-of-effect-size-toolbox'))
% ----------------------------------------------------------------------- %
% pRFcenter location 
Data = readtable('pRFcenter.csv');
within = table([1 2 3]','VariableNames',{'Conditions'});
rm_pRF = fitrm(Data,'Scrambled-Emotional~Subjects+Hemispheres*ROI',...
    'WithinDesign',within); 
ranovatbl_pRF = ranova(rm_pRF);
tbl_pRF = multcompare(rm_pRF,'Conditions','By',var,'Alpha', 0.05, ...
    'ComparisonType', 'bonferroni');
stats_pRF = mes1way([Data.Scrambled,Data.Neutral,Data.Emotional],'eta2');
% Value, e.g. Sigma
Data = readtable([Val '.csv']);
within = table([1 2 3]','VariableNames',{'Conditions'});
rm_val = fitrm(Data,'Scrambled-Emotional~Subjects+Hemispheres*ROI',...
    'WithinDesign',within);
ranovatbl_val = ranova(rm_val);
tbl_val = multcompare(rm_val,'Conditions','By',var,'Alpha', 0.05, ...
    'ComparisonType', 'bonferroni');
stats_val = mes1way([Data.Scrambled,Data.Neutral,Data.Emotional],'eta2');
% ----------------------------------------------------------------------- %
% Save results for rmANOVA
% ranavotbl cannot be saved as it is. So, I rearranged it.
p = table({'isDep';'nBoot';'confLevel';'eta2';'eta2Ci-min';'eta2Ci-max'}, ...
    [stats_pRF.isDep;stats_pRF.nBoot;stats_pRF.confLevel;stats_pRF.eta2;...
    stats_pRF.eta2Ci(1);stats_pRF.eta2Ci(2)]);
ranovatbl_pRF = table(ranovatbl_pRF);
ranovatbl_pRF=[ranovatbl_pRF.Properties.RowNames,
    ranovatbl_pRF.ranovatbl_pRF(1,:);ranovatbl_pRF.ranovatbl_pRF(2,:);...
    ranovatbl_pRF.ranovatbl_pRF(3,:);ranovatbl_pRF.ranovatbl_pRF(4,:);...
    ranovatbl_pRF.ranovatbl_pRF(5,:);ranovatbl_pRF.ranovatbl_pRF(6,:)];
ranovatbl_pRF = [ranovatbl_pRF,p];
% ----------------------------------------------------------------------- %
v = table({'isDep';'nBoot';'confLevel';'eta2';'eta2Ci-min';'eta2Ci-max'}, ...
    [stats_val.isDep;stats_val.nBoot;stats_val.confLevel;stats_val.eta2;...
    stats_val.eta2Ci(1);stats_val.eta2Ci(2)]);
ranovatbl_val = table(ranovatbl_val);
ranovatbl_val=[ranovatbl_val.Properties.RowNames,
    ranovatbl_val.ranovatbl_val(1,:);ranovatbl_val.ranovatbl_val(2,:);...
    ranovatbl_val.ranovatbl_val(3,:);ranovatbl_val.ranovatbl_val(4,:);...
    ranovatbl_val.ranovatbl_val(5,:);ranovatbl_val.ranovatbl_val(6,:)];
ranovatbl_val = [ranovatbl_val,v];
cd(StatFolder)
writetable(ranovatbl_pRF, 'rmANOVA_pRFcenter.csv','WriteRowNames',true)
writetable(ranovatbl_val, ['rmANOVA_' Val '.csv'],'WriteRowNames',true)
% ----------------------------------------------------------------------- %
% Save results for multiple comparison (post-hoc)
cd(StatFolder)
writetable(tbl_val,['Posthoc_pRFcenter_' var '.csv'],'WriteRowNames',true)
writetable(tbl_val,['Posthoc_' Val '_' var '.csv'],'WriteRowNames',true)

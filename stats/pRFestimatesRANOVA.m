% Repeated Measures ANOVA of the pRF estimations
% Written by C. Yilmaz, Boyaci Lab, 2019
% MATLAB R2016b
% This code uses the outputs of pRFestimates.m
% Uncomment rmANOVA or post-hoc for analyses. Find them into the code.
% You can save the results as csv by uncommenting the last two parts of code.
% ----------------------------------------------------------------------- %
clear all; close all;
Val = 'Beta'; % Sigma: pRF size / Beta: BOLD amplitude
StatFolder = '/Users/cemreyilmaz/Documents/fear-pRF/results';
EstFolder = '/Users/cemreyilmaz/Documents/fear-pRF/data/estimates';
cd(EstFolder)
var = 'Hemispheres'; % post-hoc variable: ROI / Hemispheres
% ----------------------------------------------------------------------- %
% pRFcenter location 
Data = readtable('pRFcenter.csv');
within = table([1 2 3]','VariableNames',{'Conditions'});
rm_pRF = fitrm(Data,'Scrambled-Emotional~Subjects+Hemispheres*ROI',...
    'WithinDesign',within); 
ranovatbl_pRF = ranova(rm_pRF);
tbl_pRF = multcompare(rm_pRF,var,'Alpha', 0.05, ...
    'ComparisonType', 'bonferroni');

% Value, e.g. Sigma
Data = readtable([Val '.csv']);
within = table([1 2 3]','VariableNames',{'Conditions'});
rm_val = fitrm(Data,'Scrambled-Emotional~Subjects+Hemispheres*ROI',...
    'WithinDesign',within);
ranovatbl_val = ranova(rm_val);
tbl_val = multcompare(rm_val,var,'Alpha', 0.05, ...
    'ComparisonType', 'bonferroni');
% ----------------------------------------------------------------------- %
% Save results for rmANOVA
% ranavotbl cannot be saved as it is. So, I rearranged it.
ranovatbl_pRF = table(ranovatbl_pRF);
ranovatbl_pRF=[ranovatbl_pRF.Properties.RowNames,
    ranovatbl_pRF.ranovatbl_pRF(1,:);
    ranovatbl_pRF.ranovatbl_pRF(2,:);...
    ranovatbl_pRF.ranovatbl_pRF(3,:)];
ranovatbl_val = table(ranovatbl_val);
ranovatbl_val=[ranovatbl_val.Properties.RowNames,
    ranovatbl_val.ranovatbl_val(1,:);
    ranovatbl_val.ranovatbl_val(2,:);...
    ranovatbl_val.ranovatbl_val(3,:)];
cd(StatFolder)
writetable(ranovatbl_pRF, 'rmANOVA_pRFcenter.csv')
writetable(ranovatbl_val, ['rmANOVA_' Val '.csv'])
% ----------------------------------------------------------------------- %
% Save results for multiple comparison (post-hoc)
cd(StatFolder)
writetable(tbl_val,['Posthoc_pRFcenter_' var '.csv'])
writetable(tbl_val,['Posthoc_' Val '_' var '.csv'])

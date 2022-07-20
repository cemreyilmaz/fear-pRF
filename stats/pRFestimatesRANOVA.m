% Repeated Measures ANOVA of the pRF estimations
% Written by C. Yilmaz, Boyaci Lab, 2019
% MATLAB R2016b
% This code uses the outputs of pRFestimates.m
% Uncomment rmANOVA or post-hoc for analyses. Find them into the code.
% You can save the results as csv by uncommenting the last two parts of code.
% ----------------------------------------------------------------------- %
clear all; close all;
Val = 'Sigma'; % Sigma: pRF size, Beta: BOLD amplitude
StatFolder = '/Users/cemreyilmaz/Documents/fear-pRF/results';
EstFolder = '/Users/cemreyilmaz/Documents/fear-pRF/data/estimates';
cd(EstFolder)
var = 'ROI'; % post-hoc variable
% ----------------------------------------------------------------------- %
% pRFcenter location 
Data = readtable('pRFcenter.csv');
within = table([1 2 3]','VariableNames',{'Conditions'});
disp('fitting rm: pRFcenter')
rm_pRF = fitrm(Data,'Scrambled-Emotional~Subjects+Hemispheres*ROI',...
    'WithinDesign',within); 
% uncomment if you are performing rmANOVA
disp('doing ranova: pRFcenter')
ranovatbl_pRF = ranova(rm_pRF);
disp('done: pRFcenter_SE')
% uncomment if you are performing post-hoc
disp('doing post-hoc: pRFcenter')
tbl_pRF = multcompare(rm_pRF,var);
disp('done: pRFcenter')

% Value, e.g. Sigma
Data = readtable([Val '.csv']);
within = table([1 2 3]','VariableNames',{'Conditions'});
disp('fitting rm: sigma_SE')
rm_val = fitrm(Data,'Scrambled-Emotional~Subjects+Hemispheres*ROI',...
    'WithinDesign',within);
% uncomment if you are performing rmANOVA
disp(['doing ranova: ' Val])
ranovatbl_val = ranova(rm_val);
disp(['done: ' Val])
% uncomment if you are performing post-hoc
disp(['doing post-hoc: ' Val])
tbl_val = multcompare(rm_val,var);
disp(['done: ' Val])
% ----------------------------------------------------------------------- %
% Save results for rmANOVA
cd(StatFolder)
Rows = {'Conditions'; 'Subj:Cond';'Hemi:Cond';'ROI:Cond';...
    'Hemi:ROI:Cond';'Error'};

SumSq = table2array(ranovatbl_pRF(:,1));
DF = table2array(ranovatbl_pRF(:,2));
MeanSq = table2array(ranovatbl_pRF(:,3));
F = table2array(ranovatbl_pRF(:,4));
pValue = table2array(ranovatbl_pRF(:,5));
pValueGG = table2array(ranovatbl_pRF(:,6));
pValueHF = table2array(ranovatbl_pRF(:,7));
pValueLB = table2array(ranovatbl_pRF(:,8));
T_pRFcenter = table(Rows,SumSq,DF,MeanSq,F,pValue,pValueGG,pValueHF,pValueLB);
writetable(T_pRFcenter, 'rmANOVA_pRFcenter.csv')

SumSq = table2array(ranovatbl_val(:,1));
DF = table2array(ranovatbl_val(:,2));
MeanSq = table2array(ranovatbl_val(:,3));
F = table2array(ranovatbl_val(:,4));
pValue = table2array(ranovatbl_val(:,5));
pValueGG = table2array(ranovatbl_val(:,6));
pValueHF = table2array(ranovatbl_val(:,7));
pValueLB = table2array(ranovatbl_val(:,8));
T_val = table(Rows,SumSq,DF,MeanSq,F,pValue,pValueGG,pValueHF,pValueLB);
writetable(T_val, ['rmANOVA_' Val '.csv'])
% ----------------------------------------------------------------------- %
% Save results for multiple comparison (post-hoc)
cd(StatFolder)
Difference = table2array([tbl_pRF(1,3); tbl_val(1,3)]);
StdError = table2array([tbl_pRF(1,4); tbl_val(1,4)]);
pValue = table2array([tbl_pRF(1,5); tbl_val(1,5)]);
Lower = table2array([tbl_pRF(1,6); tbl_val(1,6)]);
Upper = table2array([tbl_pRF(1,7); tbl_val(1,7)]);
Rows = {'(L-R)_pRF'; ['(L-R)_' Val]};
T = table(Rows,Difference,StdError,pValue,Lower,Upper);
writetable(T,['Posthoc_pRFcenter_' Val '_' var '.csv'])